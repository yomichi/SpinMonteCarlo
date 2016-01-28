export Lattice, dim, size, num_sites, num_bonds, neighbors, source, target
export chain_lattice, square_lattice

type Lattice
  dim :: Int
  size :: Vector{Int}
  nsites :: Int
  nbonds :: Int
  neighbors :: Matrix{Int}
  source :: Vector{Int}
  target :: Vector{Int}
end

dim(lat::Lattice) = lat.dim
size(lat::Lattice, dim::Int) = lat.size[dim]
num_sites(lat::Lattice) = lat.nsites
num_bonds(lat::Lattice) = lat.nbonds
neighbors(lat::Lattice, site::Int) = lat.neighbors[:,site]
source(lat::Lattice, bond::Int) = lat.source[bond]
target(lat::Lattice, bond::Int) = lat.target[bond]

function chain_lattice(L::Int)
  neighbors = zeros(Int,2,L)
  source = zeros(Int,L)
  target = zeros(Int,L)
  for s in 1:L
    neighbors[1,s] = mod1(s+1,L)
    neighbors[2,s] = mod1(s-1,L)
    source[s] = s
    target[s] = neighbors[1,s]
  end
  Lattice(1,[L],L,L,neighbors,source,target)
end

square_lattice(L::Int) = square_lattice(L,L)
function square_lattice(L::Int, W::Int)
  s2xy(s::Int) = mod(s-1,L),div(s-1,L)
  xy2s(x::Int, y::Int) = mod(y,W)*L + mod(x,L) + 1
  nsites = L*W
  nbonds = 2nsites
  neighbors = zeros(Int,4,nsites)
  source = zeros(Int,nbonds)
  target = zeros(Int,nbonds)
  for s in 1:nsites
    x,y = s2xy(s)
    neighbors[1,s] = xy2s(x+1,y)
    neighbors[2,s] = xy2s(x,y+1)
    neighbors[3,s] = xy2s(x-1,y)
    neighbors[4,s] = xy2s(x,y-1)
    source[2s-1:2s] = s
    target[2s-1:2s] = neighbors[1:2,s]
  end
  Lattice(2,[L,W],nsites,nbonds,neighbors,source,target)
end
