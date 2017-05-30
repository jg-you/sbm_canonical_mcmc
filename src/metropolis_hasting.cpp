#include "metropolis_hasting.h"

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Non class methods
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Implemented from
// http://www.fys.ku.dk/~andresen/BAhome/ownpapers/permanents/annealSched.pdf
double exponential_schedule(unsigned int t, float_vec_t cooling_schedule_kwargs)
{
  // kwargs is the speed of the exponential cooling.
  return cooling_schedule_kwargs[0] * std::pow(cooling_schedule_kwargs[1], t);
}
double linear_schedule(unsigned int t, float_vec_t cooling_schedule_kwargs)
{
  // kwargs are the initial temperature and a rate of linear cooling.
  return cooling_schedule_kwargs[0] - cooling_schedule_kwargs[1] * t;
}
double logarithmic_schedule(unsigned int t, float_vec_t cooling_schedule_kwargs)
{
  // kwargs are the rate of lienar cooling and a delay (typically 1).
  return cooling_schedule_kwargs[0] / std::log(t + cooling_schedule_kwargs[1]);
}
double constant_schedule(unsigned int t, float_vec_t cooling_schedule_kwargs)
{
  // kwargs are the rate of lienar cooling and a delay (typically 1).
  return cooling_schedule_kwargs[0];
}



// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// metropolis_hasting class
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bool metropolis_hasting::step(blockmodel_t& blockmodel,
                              const float_mat_t& p,
                              double temperature,
                              std::mt19937 & engine)
{
  std::vector<mcmc_move_t> moves = sample_proposal_distribution(blockmodel, engine);
  double a = std::pow(transition_ratio(blockmodel, p, moves), 1 / temperature);
  if (random_real(engine) < a || a > 1)
  {
    blockmodel.apply_mcmc_moves(moves);
    return true;
  }
  return false;
}
double metropolis_hasting::marginalize(blockmodel_t& blockmodel,
                                       uint_mat_t& marginal_distribution,
                                       const float_mat_t& p,
                                       unsigned int burn_in_time,
                                       unsigned int sampling_frequency,
                                       unsigned int num_samples,
                                       std::mt19937& engine)
{
  unsigned int accetped_steps = 0;
  // Burn-in period
  for (unsigned int t = 0; t < burn_in_time; ++t)
  {
    step(blockmodel, p, 1.0, engine);
  }
  // Sampling
  for (unsigned int t = 0; t < sampling_frequency * num_samples; ++t)
  {
    if (t % sampling_frequency == 0)
    {
      // Sample the blockmodel
      uint_vec_t memberships = blockmodel.get_memberships();
      #if OUTPUT_HISTORY == 1 // compile time output
      output_vec<uint_vec_t>(memberships, std::cout);
      #endif
      for (unsigned int i = 0; i < blockmodel.get_N(); ++i)
      {
        marginal_distribution[i][memberships[i]] += 1;
      }
    }
    if (step(blockmodel, p, 1.0, engine))
    {
      ++accetped_steps;
    }
  }
  return (double) accetped_steps / ((double) sampling_frequency * num_samples);
}
void metropolis_hasting::anneal(blockmodel_t& blockmodel,
                                  const float_mat_t& p,
                                  double (*cooling_schedule)(unsigned int, float_vec_t),
                                  float_vec_t cooling_schedule_kwargs,
                                  unsigned int duration,
                                  std::mt19937& engine)
{
  for (unsigned int t = 0; t < duration; ++t)
  {
    #if OUTPUT_HISTORY == 1  // compile time output
    output_vec<uint_vec_t>(blockmodel.get_memberships(), std::cout);
    #endif
    step(blockmodel, p, cooling_schedule(t, cooling_schedule_kwargs), engine);
  }
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// virtual functions implementation
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/* Implementation for the single vertex change (SBM) */
std::vector<mcmc_move_t> mh_single_vertex_sbm::sample_proposal_distribution(blockmodel_t& blockmodel, std::mt19937& engine)
{
    return blockmodel.single_vertex_change(engine);
}
double mh_single_vertex_sbm::transition_ratio(const blockmodel_t& blockmodel, const float_mat_t& p, const std::vector<mcmc_move_t> moves)
{
  int_vec_t ki = blockmodel.get_k(moves[0].vertex);
  int_vec_t n = blockmodel.get_size_vector();
  unsigned int r = moves[0].source;
  unsigned int s = moves[0].target;
  // Get the part of the probability associated to  block r and s.
  double a = std::pow((1 - p[s][s]) / (1 - p[r][s]), n[s] - ki[s]) *
             std::pow((1 - p[r][s]) / (1 - p[r][r]), n[r] - ki[r] - 1) *
             std::pow(p[r][s] / p[r][r], ki[r]) *
             std::pow(p[s][s] / p[r][s], ki[s]);
  // compute the rest
  if (n.size() > 2)
  {
    for (unsigned int l = 0; l < n.size(); ++l)
    {
      if (l != r && l != s)
      {
        a *= std::pow((1 - p[s][l]) / (1 - p[r][l]), n[l] - ki[l]) *
             std::pow(p[s][l] / p[r][l], ki[l]);
      }
    }
  }
  return a;
}
/* Implementation for the single vertex change (PPM) */
std::vector<mcmc_move_t> mh_single_vertex_ppm::sample_proposal_distribution(blockmodel_t& blockmodel, std::mt19937& engine)
{
  return blockmodel.single_vertex_change(engine);
}
double mh_single_vertex_ppm::transition_ratio(const blockmodel_t& blockmodel, const float_mat_t& p, const std::vector<mcmc_move_t> moves)
{
  int_vec_t ki = blockmodel.get_k(moves[0].vertex);
  int_vec_t n = blockmodel.get_size_vector();
  unsigned int r = moves[0].source;
  unsigned int s = moves[0].target;
  double a = std::pow((1 - p[0][0]) / (1 - p[0][1]), n[s] - ki[s] - n[r] + ki[r] + 1) *
             std::pow(p[0][0] / p[0][1], ki[s] - ki[r]);
  return a;
}
/* Implementation for the vertices swap (SBM) */
std::vector<mcmc_move_t> mh_vertices_swap_sbm::sample_proposal_distribution(blockmodel_t& blockmodel, std::mt19937& engine)
{
  return blockmodel.vertices_swap(engine);
}
double mh_vertices_swap_sbm::transition_ratio(const blockmodel_t& blockmodel, const float_mat_t& p, const std::vector<mcmc_move_t> moves)
{
  int_vec_t ki = blockmodel.get_k(moves[0].vertex);
  int_vec_t kj = blockmodel.get_k(moves[1].vertex);
  unsigned int r = moves[0].source;
  unsigned int s = moves[1].source;
  int a_xy = 0;
  if (blockmodel.are_connected(moves[0].vertex, moves[1].vertex))
  {
    a_xy = 1;
  }
  // Get the part of the probability associated to  block r and s.
  double a = std::pow((p[r][s] / p[r][r]) * (1 - p[r][r]) / (1 - p[r][s]), ki[r] - kj[r] + a_xy) *
             std::pow((p[s][s] / p[r][s]) * (1 - p[r][s]) / (1 - p[s][s]), ki[s] - kj[s] - a_xy);
  // compute the rest
  if (ki.size() > 2)
  {
    for (unsigned int l = 0; l < ki.size(); ++l)
    {
      if (l != r && l != s)
      {
        a *= std::pow((p[s][l] / p[r][l]) * (1 - p[r][l]) / (1 - p[s][l]), ki[l] - kj[l]);
      }
    }
  }
  return a;
}
/* Implementation for the vertices swap (PPM) */
std::vector<mcmc_move_t> mh_vertices_swap_ppm::sample_proposal_distribution(blockmodel_t& blockmodel, std::mt19937& engine)
{
    return blockmodel.vertices_swap(engine);
}
double mh_vertices_swap_ppm::transition_ratio(const blockmodel_t& blockmodel, const float_mat_t& p, const std::vector<mcmc_move_t> moves)
{
    int_vec_t ki = blockmodel.get_k(moves[0].vertex);
    int_vec_t kj = blockmodel.get_k(moves[1].vertex);
    unsigned int r = moves[0].source;
    unsigned int s = moves[1].source;
    int a_xy = 0;
    if (blockmodel.are_connected(moves[0].vertex, moves[1].vertex)) {
        a_xy = 1;
    }
    // Get the part of the probability associated to  block r and s.
    double a = std::pow((p[0][1] / p[0][0]) * (1 - p[0][0]) / (1 - p[0][1]), ki[r] - kj[r] + a_xy) *
               std::pow((p[0][0] / p[0][1]) * (1 - p[0][1]) / (1 - p[0][0]), ki[s] - kj[s] - a_xy);
    return a;
}
