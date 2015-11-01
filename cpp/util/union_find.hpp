#ifndef UNION_FIND_HPP
#define UNION_FIND_HPP

#include <vector>

namespace union_find{
struct Node{
  int parent;
  int weight;
  int id;
  Node():parent(-1), weight(1), id(-1){}
  bool isroot() const{ return parent == -1;}
};

typedef std::vector<Node> Nodes;

inline bool isroot(Nodes const& nodes, int index)
{
  return nodes[index].isroot();
}

inline int root_index(Nodes const& nodes, int index)
{
  while(!isroot(nodes, index)) index = nodes[index].parent;
  return index;
}

inline Node root(Nodes const& nodes, int index)
{
  return nodes[root_index(nodes, index)];
}

inline int cluster_id(Nodes const& nodes, int index)
{
  return nodes[index].id;
}

inline int unify(Nodes& nodes, int i0, int i1)
{
  int r0 = root_index(nodes, i0);
  int r1 = root_index(nodes, i1);
  if(r0 == r1) return r0;

  if(nodes[r0].weight < nodes[r1].weight){
    std::swap(r0, r1);
  }
  nodes[r1].parent = r0;
  if(nodes[r1].weight == nodes[r0].weight){
    ++nodes[r0].weight;
  }
  return r0;
}

inline int clusterize(Nodes& nodes)
{
  const int nnodes = nodes.size();
  int nc = 0;
  for(int i=0; i<nnodes; ++i){
    if(isroot(nodes, i)){
      nodes[i].id = nc;
      ++nc;
    }
  }
  for(int i=0; i<nnodes; ++i){
    nodes[i].id = root(nodes,i).id;
  }
  return nc;
}

}

#endif
