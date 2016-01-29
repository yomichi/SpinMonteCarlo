include("union_find.jl")

using Distributions
using MCObservables

square(x) = x*x

abstract LongRangeIsingChain

type LRIC_FT <: LongRangeIsingChain
  L :: Int
  T :: Float64
  beta :: Float64
  sigma :: Float64
  spins :: Vector{Int}
  Js :: Vector{Float64}
  bond_dist :: Distributions.AliasTable
  poisson_dist :: Distributions.Poisson

  function LRIC_FT(L::Integer, T::Real, sigma::Real)
    beta = 1.0/T
    Js = zeros(L-1)
    if sigma > 1.0
      invLs = 1.0/(L^sigma)
      @inbounds for i in 1:L-1
        r = i/L
        Js[i] = zeta(sigma, r) * invLs
      end
    elseif sigma == 0.0
      Js[1:end] = 0.5/L
    else
      error("sigma should be zero (Infinite range) or more than one")
    end
    lambda_tot = sum(Js)
    ps = Js ./ lambda_tot
    bond_dist = sampler(Categorical(ps))
    lambda_tot *= 2beta * L
    poisson_dist = Distributions.Poisson(lambda_tot)
    spins = rand([1,-1],L)
    return new(L, T, beta, sigma, spins, Js, bond_dist, poisson_dist)
  end
end

type LRIC_naiveSW <: LongRangeIsingChain
  L :: Int
  T :: Float64
  beta :: Float64
  sigma :: Float64
  spins :: Vector{Int}
  Js :: Vector{Float64}
  ps :: Vector{Float64}

  function LRIC_naiveSW(L::Integer, T::Real, sigma::Real)
    beta = 1.0/T

    Js = zeros(L-1)
    if sigma > 1.0
      invLs = 1.0/(L^sigma)
      @inbounds for i in 1:L-1
        r = i/L
        Js[i] = (zeta(sigma, r) + zeta(sigma, 1-r)) * invLs
      end
    elseif sigma == 0.0
      Js[1:end] = 1.0/L
    else
      error("sigma should be zero (Infinite range) or more than one")
    end

    ps = Float64[ -expm1(-2beta*J) for J in Js]

    spins = rand([1,-1],L)
    return new(L, T, beta, sigma, spins, Js, ps)
  end
end


type Cluster
  spin :: Int
  size :: Int
end

function uf_bond(model::LRIC_FT)
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
  return uf
end
function uf_bond(model::LRIC_naiveSW)
  L = model.L
  L2 = div(L,2)
  uf = UnionFind(L)
  @inbounds for s1 in 1:L
    for s2 in (s1+1):L
      r = min(s2-s1, L-(s2-s1))
      if model.spins[s1] == model.spins[s2] && rand() < model.ps[r]
        unify!(uf, s1,s2)
      end
    end
  end
  return uf
end
function update!(model::LongRangeIsingChain)
  L = model.L
  uf = uf_bond(model)
  nc = clusterize!(uf)
  clusters = Vector{Cluster}(nc)
  @inbounds for ci in 1:nc
    clusters[ci] = Cluster(rand([1,-1]), 0)
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

  for cl in clusters
    c2 = square(cl.size/L)
    m4 += square(c2)
    m4 += 6.0 * c2 * m2
    m2 += c2
  end

  obs["Magnetization"] << 0.0
  obs["Magnetization^2"] << m2
  obs["Magnetization^4"] << m4
end

function ising_LRSW{Model<:LongRangeIsingChain}(::Type{Model};
  L::Integer=100, T::Float64=1.0, sigma::Real=2.0, Sweeps::Int=8192, Thermalization::Int=Sweeps>>3)

  model = Model(L, T, sigma)

  obs = BinningObservableSet()
  makeMCObservable(obs, "Magnetization")
  makeMCObservable(obs, "Magnetization^2")
  makeMCObservable(obs, "Magnetization^4")
  makeMCObservable(obs, "Energy")
  makeMCObservable(obs, "Energy^2")
  makeMCObservable(obs, "Time")
  makeMCObservable(obs, "Speed")

  @inbounds for mcs in 1:Thermalization
    update!(model)
  end
  @inbounds for mcs in 1:Sweeps
    tic()
    clusters = update!(model)
    measure!(obs, clusters, L)
    t = toq()
    obs["Time"] << t
    obs["Speed"] << 1.0/t
  end

  jk_obs = jackknife(obs)

  jk_obs["Binder"] = jk_obs["Magnetization^4"]/(jk_obs["Magnetization^2"]*jk_obs["Magnetization^2"])

  return jk_obs
end
