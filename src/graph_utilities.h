#ifndef GRAPH_UTILITIES_H
#define GRAPH_UTILITIES_H

#include <string>
#include <fstream>
#include <sstream>
#include "types.h"

/* Load an edge list. Result passed by reference. Returns true on success. */
bool load_edge_list(edge_list_t & edge_list, const std::string edge_list_path);
/* Convert adjacency list to edge list. Result passed by reference. */
adj_list_t edge_to_adj(const edge_list_t & edge_list, unsigned int num_vertices=0);

#endif // GRAPH_UTILITIES_H