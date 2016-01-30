include("lattice.jl")
include("union_find.jl")

using MCObservables

square(x) = x*x

type Cluster
  spin :: Int
  size :: Int
end

type Ising_SW
  lat :: Lattice
  T :: Float64
  beta :: Float64
  p :: Float64
  spins :: Vector{Int}
  clusters :: Vector{Cluster}

  function Ising_SW(lat::Lattice, T::Real)
    model = new()
    model.lat = lat
    model.T = T
    model.beta = 1.0/T
    model.p = -expm1(-2model.beta)
    model.spins = rand([1,-1], num_sites(lat))
    return model
  end
end

function update!(model::Ising_SW)
  uf = UnionFind(nsites)
  @inbounds for bond in 1:nbonds
    s1,s2 = source(lat, bond), target(lat, bond)
    if spins[s1] == spins[s2] && rand() < p
      unify(uf, s1,s2)
    end
  end
  nc = clusterize!(uf)
  clusters = Cluster[ Cluster(rand([1,-1]), 0) for i in 1:nc]
  @inbounds for site in 1:nsites
    id = cluster_id(uf, site)
    cl = clusters[id]
    spins[site] = cl.spin
    cl.size += 1
  end
  model.clusters = clusters
  return nothing
end

function measure!(obs::MCObservable, model::Ising_SW)
  N = num_sites(model.lat)
  m2 = 0.0
  m4 = 0.0
  for cl in model.clusters
    c2 = square(cl.size/N)
    m4 += square(c2)
    m4 += 6.0 * m2 * c2
    m2 += c2
  end
  obs["Magnetization"] << 0.0
  obs["Magnetization^2"] << m2
  obs["Magnetization^4"] << m4
end

function ising_SW(lat::Lattice; T::Real=1.0, Sweeps::Int=8192, Thermalization::Int=Sweeps>>3)
  model = Ising_SW(lat, T)

  obs = BinningObservableSet()
  makeMCObservable(obs, "Magnetization")
  makeMCObservable(obs, "Magnetization^2")
  makeMCObservable(obs, "Magnetization^4")

  for mcs in 1:Thermalization
    update!(model)
  end
  for mcs in 1:Sweeps
    update!(model)
    measure!(obs, model)
  end

  jk_obj = jackknife(obs)
  jk_obj["Binder"] = jk_obj["Magnetization^4"] / (jk_obj["Magnetization^2"]^2)

  return obs
end

