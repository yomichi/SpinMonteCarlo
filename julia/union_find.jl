export UnionFind, is_root, root, unify, new_node, num_nodes, num_clusters, reset

type UnionFind
  parents :: Vector{Int}
  nnodes :: Int
  nclusters :: Int
  UnionFind() = new(zeros(Int,0),0,0)
  UnionFind(n::Int) = new(reshape(1:n,n),n,n)
end

num_nodes(u::UnionFind) = u.nnodes
num_clusters(u::UnionFind) = u.nclusters

is_root(u::UnionFind, n::Int) = u.parents[n] == n

function root_and_weight(u::UnionFind, n::Int)
  r,w = n,0
  while !is_root(u,r)
    r = u.parents[r]
    w += 1
  end
  u.parents[n] = r
  return r,w
end

root(u::UnionFind, n::Int) = root_and_weight(u,n)[1]

function unify(u::UnionFind, n1::Int, n2::Int)
  r1,w1 = root_and_weight(u,n1)
  r2,w2 = root_and_weight(u,n2)
  if r1 != r2
    u.nclusters -= 1
    if w1<w2
      u.parents[r1] = r2
    else
      u.parents[r2] = r1
    end
  else
    return r1
  end
end

function new_node(u::UnionFind) 
  push!(u.parents,length(u.parents)+1)
  u.nclusters += 1
  u.nnodes += 1
end

