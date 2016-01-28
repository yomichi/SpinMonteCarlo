include("observable.jl")
include("lattice.jl")
include("union_find.jl")

type Cluster
  spin :: Int
  size :: Int
end

function update!(spins::AbstractArray{Int}, p::Real)
  uf = UnionFind(nsites)
  for bond in 1:nbonds
    s1,s2 = source(lat, bond), target(lat, bond)
    if spins[s1] == spins[s2] && rand() < p
      unify(uf, s1,s2)
    end
  end
  nc = clusterize!(uf)
  clusters = Vector{Cluster}(nc)
  for cl in clusters
    cl.size = 0
    cl.spin = rand([1,-1])
  end
  for site in 1:nsites
    id = cluster_id(uf, site)
    cl = clusters[id]
    spins[site] = cl.spin
    cl.size += 1
  end
  return clusters
end

function ising_SW(lat::Lattice, T::Float64, Sweeps::Int, Thermalization::Int)
  nsites = num_sites(lat)
  nbonds = num_bonds(lat)
  spins = ones(Int,nsites)

  p = 1.0-exp(-2.0/T)

  obs = ObservableSet()
  add(obs,"Magnetization")
  add(obs,"Magnetization^2")
  add(obs,"Magnetization^4")
  add(obs,"Energy")
  add(obs,"Energy^2")

  function measure(clusters)
    m2 = 0.0
    m4 = 0.0

    for i in 1:length(clusters)
      c2 = square(clusters[i].size)
      m4 += square(mi)
      m4 += 6.0 * c2 * m2
      m2 += mi
    end

    add(obs["Magnetization"], 0.0)
    add(obs["Magnetization^2"], m2)
    add(obs["Magnetization^4"], m4)
    #add(obs["Energy"], en)
    #add(obs["Energy^2"], square(en))
  end

  for mcs in 1:Thermalization
    update()
  end
  for mcs in 1:Sweeps
    measure(update())
  end

  return obs
end

L = 16
lat = square_lattice(L)
Sweeps = 8192
Thermalization = Sweeps>>3

for T in 2.0:0.1:2.4
  observables = ising_SW(lat, T, Sweeps, Thermalization)
  for obs in observables
    open("$(obs[1])_$L.dat", "a") do os
      write(os, "$T $(mean(obs[2])) $(error(obs[2]))\n")
    end
  end
end

