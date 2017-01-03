#ifndef BLOCKMODEL_H
#define BLOCKMODEL_H

#include <random>
#include <utility>
#include <algorithm> // std::shuffle
#include <vector>
#include "types.h"

class blockmodel_t {
public:
  blockmodel_t(const uint_vec_t & memberships, unsigned int g, unsigned int N, adj_list_t * adj_list_ptr);

  std::vector<mcmc_move_t> single_vertex_change(std::mt19937& engine);
  std::vector<mcmc_move_t> vertices_swap(std::mt19937& engine);

  int_vec_t get_k(unsigned int vertex) const;
  bool are_connected(unsigned int vertex_a, unsigned int vertex_b) const;
  int_vec_t get_size_vector() const;
  uint_vec_t get_memberships() const;
  uint_mat_t get_m() const;
  unsigned int get_N() const;
  unsigned int get_g() const;

  void apply_mcmc_moves(std::vector<mcmc_move_t> moves);

  void shuffle(std::mt19937& engine);

private:
    /// State variable
  adj_list_t * adj_list_ptr_;
  int_mat_t k_;
  int_vec_t n_;
  uint_vec_t memberships_;
    /// Internal distribution. Generator must be passed as a service
  std::uniform_int_distribution<> random_block_;
  std::uniform_int_distribution<> random_node_;
    /// Private methods
    /* Compute the degree matrix from scratch. */
  void compute_k();
};

#endif // BLOCKMODEL_H