#ifndef METROPOLIS_HASTING_H
#define METROPOLIS_HASTING_H

#include <cmath>
#include <vector>
#include <iostream>
#include "types.h"
#include "blockmodel.h"
#include "output_functions.h"

/* Cooling schedules */
double exponential_schedule(unsigned int t, float_vec_t cooling_schedule_kwargs);
double linear_schedule(unsigned int t, float_vec_t cooling_schedule_kwargs);
double logarithmic_schedule(unsigned int t, float_vec_t cooling_schedule_kwargs);
double constant_schedule(unsigned int t, float_vec_t cooling_schedule_kwargs);

class metropolis_hasting
{
protected:
  std::uniform_real_distribution<> random_real;
public:
  // Ctor
  metropolis_hasting() : random_real(0,1) {;}

  // Virtual methods
  virtual std::vector<mcmc_move_t> sample_proposal_distribution(blockmodel_t& blockmodel, std::mt19937& engine)
    {return std::vector<mcmc_move_t>();}  // bogus virtual implementation
  virtual double transition_ratio(const blockmodel_t& blockmodel, const float_mat_t & p, const std::vector<mcmc_move_t> moves)
    {return 0;}  // bogus virtual implementation

  // Common methods
  bool step(blockmodel_t& blockmodel,
            const float_mat_t & p,
            double temperature,
            std::mt19937 & engine);
  double marginalize(blockmodel_t& blockmodel,
                     uint_mat_t & marginal_distribution,
                     const float_mat_t & p,
                     unsigned int burn_in_time,
                     unsigned int sampling_frequency,
                     unsigned int num_samples,
                     std::mt19937& engine);
  double anneal(blockmodel_t& blockmodel,
                const float_mat_t & p,
                double (*cooling_schedule)(unsigned int, float_vec_t),
                float_vec_t cooling_schedule_kwargs,
                unsigned int duration,
                std::mt19937& engine);
};

/* Inherited classes with specific definitions */
class mh_single_vertex_sbm : public metropolis_hasting
{
public:
  std::vector<mcmc_move_t> sample_proposal_distribution(blockmodel_t& blockmodel, std::mt19937& engine);
  double transition_ratio(const blockmodel_t& blockmodel, const float_mat_t & p, const std::vector<mcmc_move_t> moves);
};

class mh_single_vertex_ppm : public metropolis_hasting
{
public:
  std::vector<mcmc_move_t> sample_proposal_distribution(blockmodel_t& blockmodel, std::mt19937& engine);
  double transition_ratio(const blockmodel_t& blockmodel, const float_mat_t & p, const std::vector<mcmc_move_t> moves);
};

class mh_vertices_swap_sbm : public metropolis_hasting
{
public:
  std::vector<mcmc_move_t> sample_proposal_distribution(blockmodel_t& blockmodel, std::mt19937& engine);
  double transition_ratio(const blockmodel_t& blockmodel, const float_mat_t & p, const std::vector<mcmc_move_t> moves);
};

class mh_vertices_swap_ppm : public metropolis_hasting
{
public:
  std::vector<mcmc_move_t> sample_proposal_distribution(blockmodel_t& blockmodel, std::mt19937& engine);
  double transition_ratio(const blockmodel_t& blockmodel, const float_mat_t & p, const std::vector<mcmc_move_t> moves);
};

#endif // METROPOLIS_HASTING_H