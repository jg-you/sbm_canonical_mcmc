#include "blockmodel.h"


blockmodel_t::blockmodel_t(const uint_vec_t & memberships, unsigned int g, unsigned int N, adj_list_t * adj_list_ptr) :
random_block_(0, g - 1),
random_node_(0, N - 1)
{
  memberships_ = memberships;
  adj_list_ptr_ = adj_list_ptr;
  n_.resize(g, 0);
  for (unsigned int j = 0; j < memberships.size(); ++j)
  {
    ++n_[memberships[j]];
  }
  compute_k();
}


std::vector<mcmc_move_t> blockmodel_t::single_vertex_change(std::mt19937& engine)
{
  std::vector<mcmc_move_t> moves(1);
  moves[0].vertex = random_node_(engine);
  moves[0].source = memberships_[moves[0].vertex];
  moves[0].target = random_block_(engine);
  return moves;
}
std::vector<mcmc_move_t> blockmodel_t::vertices_swap(std::mt19937& engine)
{
  std::vector<mcmc_move_t> moves(2);
  moves[0].vertex = random_node_(engine);
  moves[1].vertex = random_node_(engine);
  moves[0].source = memberships_[moves[0].vertex];
  moves[0].target = memberships_[moves[1].vertex];
  moves[1].source = memberships_[moves[1].vertex];
  moves[1].target = memberships_[moves[0].vertex];
  return moves;
}

int_vec_t blockmodel_t::get_k(unsigned int vertex) const {return k_[vertex];}
bool blockmodel_t::are_connected(unsigned int vertex_a, unsigned int vertex_b) const
{
  if (adj_list_ptr_->at(vertex_a).find(vertex_b) != adj_list_ptr_->at(vertex_a).end())
  {
    return true;
  }
  return false;
}
int_vec_t blockmodel_t::get_size_vector() const {return n_;}
uint_vec_t blockmodel_t::get_memberships() const {return memberships_;}
uint_mat_t blockmodel_t::get_m() const
{
  uint_mat_t m(get_g(), uint_vec_t(get_g(), 0));
  for (auto vertex = 0; vertex < adj_list_ptr_->size(); ++vertex)
  {
    for (auto neighbour = adj_list_ptr_->at(vertex).begin(); neighbour != adj_list_ptr_->at(vertex).end(); ++neighbour)
    {
      ++m[memberships_[vertex]][memberships_[*neighbour]];
    }
  }
  for (unsigned int r = 0; r < get_g(); ++r)
  {
    for (unsigned int s = 0; s < get_g(); ++s)
    {
      m[r][s] /= 2;  // edges are counted twice (the adj_list is symmetric)
      m[r][s] = m[s][r];  // symmetrize m matrix.
    }
  }
  return m;
}
unsigned int blockmodel_t::get_N() const {return memberships_.size();}
unsigned int blockmodel_t::get_g() const {return n_.size();}

void blockmodel_t::apply_mcmc_moves(std::vector<mcmc_move_t> moves)
{
  for (unsigned int i = 0; i < moves.size(); ++i)
  {
        // Change block degrees and block sizes
    for (auto neighbour = adj_list_ptr_->at(moves[i].vertex).begin();
     neighbour != adj_list_ptr_->at(moves[i].vertex).end();
     ++neighbour)
    {
      --k_[*neighbour][moves[i].source];
      ++k_[*neighbour][moves[i].target];
    }
    --n_[moves[i].source];
    ++n_[moves[i].target];
        // Set new memberships
    memberships_[moves[i].vertex] = moves[i].target;
  }
}


void blockmodel_t::shuffle(std::mt19937& engine)
{
  std::shuffle(memberships_.begin(), memberships_.end(), engine);
  compute_k();
}

void blockmodel_t::compute_k()
{
  k_.clear();
  k_.resize(adj_list_ptr_->size());
  for (unsigned int i = 0; i < adj_list_ptr_->size(); ++i)
  {
    k_[i].resize(this->n_.size(), 0);
    for (auto nb = adj_list_ptr_->at(i).begin(); nb != adj_list_ptr_->at(i).end(); ++nb)
    {
      ++k_[i][memberships_[*nb]];
    }
  }
}
