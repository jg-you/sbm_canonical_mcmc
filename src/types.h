#ifndef TYPES_H
#define TYPES_H

#include <vector>
#include <set>
#include <utility>

typedef std::pair<unsigned int, unsigned int> edge_t;
typedef std::vector<edge_t> edge_list_t;
typedef std::set<unsigned int> neighbourhood_t;
typedef std::vector<neighbourhood_t> adj_list_t;

typedef struct mcmc_move_t
{
  unsigned int vertex;
  unsigned int source;
  unsigned int target;
} mcmc_move_t;


typedef std::vector<unsigned int> uint_vec_t;
typedef std::vector<int> int_vec_t;
typedef std::vector<float> float_vec_t;
typedef std::vector< std::vector<unsigned int> > uint_mat_t;
typedef std::vector< std::vector<int> > int_mat_t;
typedef std::vector< std::vector<float> > float_mat_t;

#endif // TYPES_H