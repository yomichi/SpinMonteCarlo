#ifndef UNION_FIND_HPP
#define UNION_FIND_HPP

#include <vector>

namespace union_find{
struct Node{
  int parent;
  int size;
  int id;
  Node():parent(-1), size(1), id(-1){}
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
  return root(nodes, index).id;
}
inline int cluster_size(Nodes const& nodes, int index)
{
  return root(nodes, index).size;
}

inline int unify(Nodes& nodes, int i0, int i1)
{
  int r0 = root_index(nodes, i0);
  int r1 = root_index(nodes, i1);
  if(r0 == r1) return r0;

  nodes[r1].parent = r0;
  nodes[r0].size += nodes[r1].size;
  return r0;
}

}

#endif
