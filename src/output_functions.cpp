#include "output_functions.h"

void output_edge_list(const edge_list_t & edge_list, std::ostream & stream)
{
  for (auto edge = edge_list.begin(); edge != edge_list.end(); ++edge)
  {
    stream << edge->first << " " << edge->second << "\n"; 
  }
  return;
}

void output_adj_list(const adj_list_t & adj_list, std::ostream & stream)
{
  for (unsigned int n = 0; n < adj_list.size(); ++n)
  {
    stream << n << " : ";
    for (auto neighbour = adj_list[n].begin(); neighbour != adj_list[n].end(); ++neighbour)
    {
      stream << *neighbour << " ";
    }
    stream << "\n";
  }
  return;
}