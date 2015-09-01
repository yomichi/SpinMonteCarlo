require("observable.jl")
require("lattice.jl")

using UnionFinds

type cluster
  spin :: Int
  size :: Int
end

function ising_SW(lat::Lattice, T::Float64, Sweeps::Int, Thermalization::Int)
  nsites = num_sites(lat)
  nbonds = num_bonds(lat)
  spins = ones(Int,nsites)

  p = 1.0-exp(-2.0/T)

  function update()
    uf = UnionFind(nsites)
    for bond in 1:nbonds
      s1,s2 = source(lat, bond), target(lat, bond)
      if spins[s1] == spins[s2] && rand() < p
        unify(uf, s1,s2)
      end
    end
    clusters = Dict{Int,cluster}()
    for site in 1:nsites
      r = root(uf,site)
      if !has(clusters,r)
        clusters[r] = cluster(2(rand(0:1))-1, 0)
      end
      cl = clusters[r]
      spins[site] = cl.spin
      cl.size += 1
    end
    return (Float64)[ square(cl[2].size/nsites) for cl in clusters]
  end

  obs = ObservableSet()
  add(obs,"Magnetization")
  add(obs,"Magnetization^2")
  add(obs,"Magnetization^4")
  add(obs,"Energy")
  add(obs,"Energy^2")

  function measure(ar::Vector)
    m2 = 0.0
    m4 = 0.0

    for i in 1:length(ar)
      mi = ar[i]
      m2 += mi
      m4 += square(mi)
      for j in 1:i-1
        m4 += 6.0*mi*ar[j]
      end
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

