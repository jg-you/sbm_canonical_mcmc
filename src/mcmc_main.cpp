/* ~~~~~~~~~~~~~~~~ Notes ~~~~~~~~~~~~~~~~ */
// A few base assumption go into this program:
// 
// - We assume that node identifiers are zero indexed contiguous integers.
// - We assume that block memberships are zero indexed contiguous integers.
// - We assume that the SBM is of the undirected and simple variant.

/* ~~~~~~~~~~~~~~~~ Includes ~~~~~~~~~~~~~~~~ */
// STL
#include <iostream>
#include <chrono>
#include <vector>
#include <utility>
#include <random>
#include <string>
// Boost
#include <boost/program_options.hpp>
// Program headers
#include "types.h"
#include "blockmodel.h"
#include "output_functions.h"
#include "metropolis_hasting.h"
#include "graph_utilities.h"
#include "config.h"

namespace po = boost::program_options;


int main(int argc, char const *argv[]) {
    /* ~~~~~ Program options ~~~~~~~*/
    std::string edge_list_path;
    float_vec_t probabilities;
    uint_vec_t n;
    unsigned int burn_in;
    unsigned int sampling_steps;
    unsigned int sampling_frequency;
    bool randomize = false;
    bool use_ppm = false;
    bool use_single_vertex = false;
    bool maximize = false;
    std::string cooling_schedule;
    float_vec_t cooling_schedule_kwargs(2,0);
    unsigned int seed = 0;

    po::options_description description("Options");
    description.add_options()
    ("edge_list_path,e", po::value<std::string>(&edge_list_path), 
        "Path to edge list file.")
    ("probabilities,P", po::value<float_vec_t>(&probabilities)->multitoken(), 
        "In normal mode (SBM): probability matrix in row major order. In PPM mode: p_in followed by p_out.")
    ("n,n", po::value<uint_vec_t>(&n)->multitoken(), 
        "Block sizes vector.\n")
    ("burn_in,b", po::value<unsigned int>(&burn_in)->default_value(1000),
        "Burn-in time.")
    ("sampling_steps,t", po::value<unsigned int>(&sampling_steps)->default_value(1000),
        "Number of sampling steps in marginalize mode. Length of the simulated annealing process.")
    ("sampling_frequency,f", po::value<unsigned int>(&sampling_frequency)->default_value(10),
        "Number of step between each sample in marginalize mode. Unused in likelihood maximization mode.")
    ("randomize,r",
        "Randomize initial block state.")
    ("use_ppm,u",
        "Use PPM transition ratios (defaults to SBM).")
    ("use_single_vertex,s",
        "Use single vertex proposal distribution (defaults to vertices swap).")
    ("maximize,m",
        "Maximize likelihood instead of marginalizing.")
    ("cooling_schedule,c", po::value<std::string>(&cooling_schedule)->default_value("exponential"),
        "Cooling schedule for the simulated annealing algorithm. Options are exponential, linear, logarithmic and constant.")
    ("cooling_schedule_kwargs,a", po::value<float_vec_t>(&cooling_schedule_kwargs)->multitoken(),
        "Additional arguments for the cooling schedule provided as a list of floats. "\
        "Depends on the choice of schedule:\n"\
         "Exponential: T_0 (init. temperature > 0)\n"\
         "             alpha (in ]0,1[).\n"\
         "Linear: T_0 (init. temperature > 0)\n"\
         "        eta (rate of decline).\n"\
         "Logarithmic: c (rate of decline)\n"\
         "             d (delay > 1)\n"\
         "Constant: T (temperature > 0)")
    ("seed,d", po::value<unsigned int>(&seed),
        "Seed of the pseudo random number generator (Mersenne-twister 19937). A random seed is used if seed is not specified.")
    ("help,h", "Produce this help message.")
    ;
    po::variables_map var_map;
    po::store(po::parse_command_line(argc,argv,description), var_map);
    po::notify(var_map);
    if (var_map.count("help") > 0 || argc == 1) {
	#if OUTPUT_HISTORY == 0
		std::cout << "MCMC algorithms for the SBM (final output only)\n";
       	#else
		std::cout << "MCMC algorithms for the SBM (output intermediate states)\n";
	#endif
	std::cout << "Usage:\n"
                  << "  "+std::string(argv[0])+" [--option_1=value] [--option_s2=value] ...\n";
        std::cout << description;
        return 0;
    }
    if (var_map.count("edge_list_path") == 0) {
        std::cout << "edge_list_path is required (-e flag)\n";
        return 1;
    }
    if (var_map.count("n") == 0) {
        std::cout << "n is required (-n flag)\n";
        return 1;
    }
    if (var_map.count("randomize") > 0) {
        randomize = true;
    }
    if (var_map.count("use_ppm") > 0) {
        use_ppm = true;
    }
    if (var_map.count("use_single_vertex") > 0) {
        use_single_vertex = true;
    }
    if (var_map.count("maximize") > 0) {
        maximize = true;
        if (var_map.count("cooling_schedule_kwargs") == 0)
        {
          // defaults
          if (cooling_schedule == "exponential")
          {
            cooling_schedule_kwargs[0] = 1;
            cooling_schedule_kwargs[1] = 0.99;
          } 
          if (cooling_schedule == "linear")
          {
            cooling_schedule_kwargs[0] = sampling_steps + 1;
            cooling_schedule_kwargs[1] = 1;
          }
          if (cooling_schedule == "logarithmic")
          {
            cooling_schedule_kwargs[0] = 1;
            cooling_schedule_kwargs[1] = 1;
          }
          if (cooling_schedule == "constant")
          {
            cooling_schedule_kwargs[0] = 1;
          }
        }
        else
        {
          // kwards not defaulted, must check.
          if (cooling_schedule == "exponential")
          {
            if (cooling_schedule_kwargs[0] <= 0)
            {
              std::cerr << "Invalid cooling schedule argument for linear schedule: T_0 must be grater than 0.\n";
              std::cerr << "Passed value: T_0=" << cooling_schedule_kwargs[0] << "\n";
              return 1;
            }
            if (cooling_schedule_kwargs[1] <= 0 || cooling_schedule_kwargs[1] >= 1)
            {
              std::cerr << "Invalid cooling schedule argument for exponential schedule: alpha must be in ]0,1[.\n";
              std::cerr << "Passed value: alpha=" << cooling_schedule_kwargs[1] << "\n";
              return 1;
            }
          }
          else if (cooling_schedule == "linear")
          {
            if (cooling_schedule_kwargs[0] <= 0)
            {
              std::cerr << "Invalid cooling schedule argument for linear schedule: T_0 must be grater than 0.\n";
              std::cerr << "Passed value: T_0=" << cooling_schedule_kwargs[0] << "\n";
              return 1;
            }
            if (cooling_schedule_kwargs[1] <= 0 || cooling_schedule_kwargs[1] > cooling_schedule_kwargs[0])
            {
              std::cerr << "Invalid cooling schedule argument for linear schedule: eta must be in ]0, T_0].\n";
              std::cerr << "Passed value: T_0=" << cooling_schedule_kwargs[0]
                        << ", eta=" << cooling_schedule_kwargs[1] << "\n";
              return 1;
            }
            if (cooling_schedule_kwargs[1] * sampling_steps > cooling_schedule_kwargs[0])
            {
              std::cerr << "Invalid cooling schedule argument for linear schedule: eta * sampling_steps must be smaller or equal to T_0.\n";
              std::cerr << "Passed value: eta*sampling_steps=" << cooling_schedule_kwargs[1] * sampling_steps
                        << ", T_0=" << cooling_schedule_kwargs[0] << "\n";
              return 1;
            }
          }
          else if (cooling_schedule == "logarithmic")
          {
            if (cooling_schedule_kwargs[0] <= 0)
            {
              std::cerr << "Invalid cooling schedule argument for logarithmic schedule: c must be greater than 0.\n";
              std::cerr << "Passed value: c=" << cooling_schedule_kwargs[0] << "\n";
              return 1;
            }
            if (cooling_schedule_kwargs[1] <= 0)
            {
              std::cerr << "Invalid cooling schedule argument for logarithmic schedule: d must be greater than 0.\n";
              std::cerr << "Passed value: d=" << cooling_schedule_kwargs[1] << "\n";
              return 1;
            }
          }
          else if (cooling_schedule == "constant")
          {
            if (cooling_schedule_kwargs[0] <= 0)
            {
              std::cerr << "Invalid cooling schedule argument for constant schedule: temperature must be greater than 0.\n";
              std::cerr << "Passed value: T=" << cooling_schedule_kwargs[0] << "\n";
              return 1;
            }
          }
          else
          {
            std::cerr << "Invalid cooling schedule. Options are exponential, linear, logarithmic.\n";
            return 1;
          }
        }
    }
    if (var_map.count("seed") == 0) {
        // seeding based on the clock
        seed = (unsigned int) std::chrono::high_resolution_clock::now().time_since_epoch().count();
    }
    /* ~~~~~ Setup objects ~~~~~~~*/
    std::mt19937 engine(seed);
    // number of blocks
    unsigned int g = n.size();
    // number of vertices
    unsigned int N = 0;
    for (unsigned int i = 0; i < g; ++i) {
      N += n[i];
    }
    // Graph structure
    edge_list_t edge_list;
    load_edge_list(edge_list, edge_list_path);
    adj_list_t adj_list = edge_to_adj(edge_list, N);
    edge_list.clear();
    // memberships from block sizes
    uint_vec_t memberships_init;
    {
      unsigned int accu = 0;
      for (auto it = n.begin(); it != n.end(); ++it) {
          accu += *it;
      }
      memberships_init.resize(accu, 0);
      unsigned shift = 0;
      for (unsigned int r = 0; r < n.size(); ++r) {
          for (unsigned int i = 0; i < n[r]; ++i) {
              memberships_init[shift + i] = r;
          }
          shift += n[r];
      }
    }
    // blockmodel
    blockmodel_t blockmodel(memberships_init, g, adj_list.size(), &adj_list);
    memberships_init.clear();
    if (randomize) {
        blockmodel.shuffle(engine);
    }
    // probabilities
    float_mat_t p(g, float_vec_t(g, 0));
    if (!use_ppm)
    {
      for (unsigned int r = 0; r < g; ++r)
      {
        for (unsigned int s = 0; s < g; ++s)
        {
          p[r][s] = probabilities[s + r*g];
        }
      }
    }
    else
    {
      for (unsigned int r = 0; r < g; ++r)
      {
        p[r][r] = probabilities[0];
        for (unsigned int s = r + 1; s < g; ++s)
        {
          p[r][s] = probabilities[1];
          p[s][r] = probabilities[1];
        }
      }
    }
    // Bind proper Metropolis-Hasting algorithm
    std::shared_ptr<metropolis_hasting> algorithm;
    if (!use_ppm && use_single_vertex) {
        algorithm = std::make_shared<mh_single_vertex_sbm>();
    }
    else if (use_ppm && use_single_vertex) {
        algorithm = std::make_shared<mh_single_vertex_ppm>();    
    }
    else if (!use_ppm &&  !use_single_vertex) {
        algorithm = std::make_shared<mh_vertices_swap_sbm>();     
    }
    else if (use_ppm && !use_single_vertex) {
        algorithm = std::make_shared<mh_vertices_swap_ppm>();     
    }

    /* ~~~~~ Logging ~~~~~~~*/
    #if LOGGING == 1
    std::clog << "edge_list_path: " << edge_list_path << "\n";
    std::clog << "probabilities:\n";
    output_mat<float_mat_t>(p, std::clog);
    std::clog << "sizes (g=" << n.size() << "): ";
    for (auto it = n.begin(); it != n.end(); ++it)
        std::clog << *it << " ";
    std::clog << "\n";
    std::clog << "burn_in: " << burn_in << "\n";
    std::clog << "sampling_steps: " << sampling_steps << "\n";
    std::clog << "sampling_frequency: " << sampling_frequency << "\n";
    if (maximize) {std::clog << "maximize: true\n";}
    else {std::clog << "maximize: false\n";}
    if (use_ppm) {std::clog << "use_ppm: true\n";}
    else {std::clog << "use_ppm: false\n";}
    if (use_single_vertex) {std::clog << "use_single_vertex: true\n";}
    else {std::clog << "use_single_vertex: false\n";}
    if (randomize) {std::clog << "randomize: true\n";}
    else {std::clog << "randomize: false\n";}
    if (maximize)
    {
      std::clog << "cooling_schedule: " << cooling_schedule << "\n";
      std::clog << "cooling_schedule_kwargs: ";
      output_vec<float_vec_t>(cooling_schedule_kwargs);
    }
    std::clog << "seed: " << seed << "\n";
    #endif

    /* ~~~~~ Actual algorithm ~~~~~~~*/
    double rate = 0;
    uint_mat_t marginal(adj_list.size(), uint_vec_t(g, 0));
    if (maximize)
    {
      if (cooling_schedule == "exponential")
      {
        algorithm->anneal(blockmodel, p, &exponential_schedule, cooling_schedule_kwargs, sampling_steps, engine);
      }
      if (cooling_schedule == "linear")
      {
        algorithm->anneal(blockmodel, p, &linear_schedule, cooling_schedule_kwargs, sampling_steps, engine);
      }
      if (cooling_schedule == "logarithmic")
      {
        algorithm->anneal(blockmodel, p, &logarithmic_schedule, cooling_schedule_kwargs, sampling_steps, engine);
      }
      if (cooling_schedule == "constant")
      {
        algorithm->anneal(blockmodel, p, &constant_schedule, cooling_schedule_kwargs, sampling_steps, engine);
      }
      output_vec<uint_vec_t>(blockmodel.get_memberships(), std::cout);
    }
    else  // marginalize
    {  
      rate = algorithm->marginalize(blockmodel, marginal, p, burn_in, sampling_frequency, sampling_steps, engine);
      uint_vec_t memberships(blockmodel.get_N(), 0);
      for (unsigned int i = 0; i < blockmodel.get_N(); ++i)
      {
        unsigned int max = 0;
        for (unsigned int r = 0; r < g; ++r)
        {
          if (marginal[i][r] > max)
          {
            memberships[i] = r;
            max = marginal[i][r];
          }
        }
      }
      output_vec<uint_vec_t>(memberships, std::cout);
      std::clog << "acceptance ratio " <<  rate  <<  "\n";
    }
    return 0;
}
