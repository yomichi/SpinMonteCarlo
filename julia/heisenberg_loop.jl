include("lattice.jl")
include("union_find.jl")

square(x) = x*x

using DataStructures
using MCObservables

@enum(OperatorType, ot_diag, ot_offdiag, ot_dummy)

type LocalOperator
  bond :: Int
  time :: Float64
  typ :: OperatorType
  dID :: Int
  uID :: Int
end
LocalOperator(b,t) = LocalOperator(b,t,ot_diag,-1,-1)

isdummy(op :: LocalOperator) = op.typ == ot_dummy
istail(ops :: Cons{LocalOperator}) = isa(ops.tail, Nil)
isdummy(ops :: Cons{LocalOperator}) = isdummy(ops.head)

type Cluster
  mag :: Float64
  length :: Float64
  toflip :: Bool
end
Cluster(tf::Bool) = Cluster(0.0, 0.0, tf)

type AFHLoop
  lat :: Lattice
  T :: Float64
  beta :: Float64
  inv_op_density :: Float64
  spins :: Vector{Int}
  operators :: Cons{LocalOperator}
  clusters :: Vector{Cluster}
  function AFHLoop(lat::Lattice, T::Real)
    nsites = num_sites(lat)
    nbonds = num_bonds(lat)
    beta = 1.0/T
    inv_op_density = 2T/num_bonds(lat)
    spins = ones(Int, nsites)
    operators = list(LocalOperator(0, 0.0, ot_dummy, 0, 0),LocalOperator(0, 1.0, ot_dummy, 0,0))
    return new(lat, T, beta, inv_op_density, spins, operators)
  end
end

function update!(model::AFHLoop)
  t = randexp()*model.inv_op_density
  ops = model.operators
  @inbounds while t < 1.0 || !istail(ops.tail)
    if istail(ops) || t < ops.tail.head.time
      ## insert an operator

      b = rand(1:num_bonds(model.lat))
      s1, s2 = source(model.lat, b), target(model.lat, b)
      if model.spins[s1] != model.spins[s2]
        lop = LocalOperator(b, t)
        ops.tail = cons(lop, ops.tail)
        t += randexp()*model.inv_op_density
      else
        t += randexp()*model.inv_op_density
        continue
      end
    elseif ops.tail.head.typ == ot_diag
      ## remove an operator
      ops.tail = ops.tail.tail
      continue
    end

    next_op = ops.tail.head
    s1, s2 = source(model.lat, next_op.bond), target(model.lat, next_op.bond)
    if next_op.typ == ot_offdiag
      model.spins[s1] *= -1
      model.spins[s2] *= -1
    end
    ops = ops.tail
  end

  N = num_sites(model.lat)
  currents = collect(1:N)
  uf = UnionFind(N)
  @inbounds for op in model.operators
    if isdummy(op)
      continue
    end
    s1, s2 = source(model.lat, op.bond), target(model.lat, op.bond)
    op.dID = unify!(uf, currents[s1], currents[s2])
    op.uID = currents[s1] = currents[s2] = add_node!(uf)
  end
  @inbounds for s in 1:N
    unify!(uf, s, currents[s])
  end

  nc = clusterize!(uf)
  clusters = Cluster[Cluster(tf) for tf in rand(Bool, nc)]
  @inbounds for op in model.operators
    if isdummy(op)
      continue
    end
    uid = cluster_id(uf, op.uID)
    ucl = clusters[uid]
    did = cluster_id(uf, op.dID)
    dcl = clusters[did]
    if ucl.toflip != dcl.toflip
      op.typ = ifelse(op.typ == ot_diag , ot_offdiag , ot_diag)
    end
    ucl.length -= 2.0*op.time
    dcl.length += 2.0*op.time
  end
  @inbounds for site in 1:N
    id = cluster_id(uf, site)
    cl = clusters[id]
    model.spins[site] *= ifelse(cl.toflip, -1,1)
    cl.length += 1.0
    cl.mag += model.spins[site]
  end
  model.clusters = clusters
end

function measure!(obs, model::AFHLoop)
  m2 = 0.0
  l2 = 0.0
  for cl in model.clusters
    m2 += square(cl.mag)
    l2 += square(cl.length)
  end
  m2 *= 0.25*model.beta/num_sites(model.lat)
  l2 *= 0.25*model.beta/num_sites(model.lat)

  obs["Uniform Susceptibility"] << m2
  obs["Staggered Susceptibility"] << l2
end

function looper(lat::Lattice; T::Real=1.0, Sweeps::Integer=8192, Thermalizations::Integer=Sweeps>>3)

  model = AFHLoop(lat, T)

  obs = BinningObservableSet()
  makeMCObservable(obs, "Uniform Susceptibility")
  makeMCObservable(obs, "Staggered Susceptibility")

  for mcs in 1:Thermalizations
    update!(model)
  end
  for mcs in 1:Sweeps
    update!(model)
    measure!(obs, model)
  end
  return obs
end
