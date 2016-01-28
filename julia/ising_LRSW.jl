include("union_find.jl")

using Distributions
using MCObservables

square(x) = x*x

type LR_Ising
  L :: Integer
  T :: Float64
  beta :: Float64
  sigma :: Float64
  spins :: Vector{Int}
  Js :: Vector{Real}
  bond_dist :: Distributions.Categorical
  poisson_dist :: Distributions.Poisson

  function LR_Ising(L::Integer, T::Real, sigma::Real)
    beta = 1.0/T

    Js = zeros(L)
    invLs = 1.0/(L^sigma)
    @inbounds for i in 1:L-1
      r = i/L
      Js[i] = (zeta(sigma, r) + zeta(sigma, 1-r)) * invLs
    end
    Js[end] = 2zeta(sigma) * invLs
    lambda_tot = sum(Js)
    ps = Js ./ lambda_tot
    bond_dist = Distributions.Categorical(ps)
    lambda_tot *= 2beta
    poisson_dist = Distributions.Poisson(lambda_tot)

    spins = ones(Int,L)

    return new(L, T, beta, sigma, spins, Js, bond_dist, poisson_dist)
  end
end

type Cluster
  spin :: Int
  size :: Int
end

function update!(model::LR_Ising)
  L = model.L
  uf = UnionFind(L)
  nbonds = rand(model.poisson_dist)
  @inbounds for ibond in 1:nbonds
    s1 = rand(1:L)
    s2 = mod1(s1+rand(model.bond_dist), L)
    if model.spins[s1] == model.spins[s2]
      unify!(uf, s1,s2)
    end
  end
  nc = clusterize!(uf)
  clusters = Vector{Cluster}(nc)
  @inbounds for ci in 1:nc
    clusters[ci] = Cluster(0, rand([1,-1]))
  end
  @inbounds for site in 1:L
    id = cluster_id(uf, site)
    cl = clusters[id]
    model.spins[site] = cl.spin
    cl.size += 1
  end
  return clusters
end

function measure!(obs, clusters, L)
  m2 = 0.0
  m4 = 0.0

  @inbounds for i in 1:length(clusters)
    c2 = square(clusters[i].size/L)
    m4 += square(c2)
    m4 += 6.0 * c2 * m2
    m2 += c2
  end

  obs["Magnetization"] << 0.0
  obs["Magnetization^2"] << m2
  obs["Magnetization^4"] << m4
end


function ising_LRSW(L::Integer, T::Float64, sigma::Real, Sweeps::Int, Thermalization::Int)

  model = LR_Ising(L, T, sigma)

  obs = BinningObservableSet()
  makeMCObservable(obs, "Magnetization")
  makeMCObservable(obs, "Magnetization^2")
  makeMCObservable(obs, "Magnetization^4")
  makeMCObservable(obs, "Energy")
  makeMCObservable(obs, "Energy^2")

  @inbounds for mcs in 1:Thermalization
    update!(model)
  end
  @inbounds for mcs in 1:Sweeps
    clusters = update!(model)
    measure!(obs, clusters, L)
  end

  jk_obs = jackknife(obs)

  jk_obs["Binder"] = jk_obs["Magnetization^4"]/(jk_obs["Magnetization^2"]^2)

  return jk_obs
end
