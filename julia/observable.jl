import Base:  mean, reset, error
export AbstractObservable, Observable, ObservableSet
export add!, reset, mean, variance, error, name

square(x) = x*x

abstract AbstractObservable

type Observable <: AbstractObservable
  xs :: Vector{Float64}
  num :: Int
  sum :: Float64
  sum2 :: Float64
end

Observable() = Observable(Float64[], 0, 0.0, 0.0)
function add!(o::Observable, val::Real)
  push!(o.xs, val)
  o.num += 1
  o.sum += val
  o.sum2 += val*val
  return o
end
function reset(o::Observable)
  o.xs = Float64[]
  o.num = 0
  o.sum = o.sum2 = 0.0
  return o
end
mean(o::Observable) = o.sum / o.num
variance(o::Observable) = o.num*(o.sum2/o.num - square(mean(o)))/(o.num-1)
error(o::Observable) = sqrt(variance(o)/o.num)
name(o::Observable) = obs.name


type Jackknife <: AbstractObservable
  xs :: Vector{Float64}
  num :: Int
  sum :: Float64
  sum2 :: Float64
end

function Jackknife(o::Observable)
  n1 = 1.0/(o.num-1)
  xs = fill(o.sum*n1, o.num)
  LinAlg.BLAS.axpy!(-n1, o.xs, xs)
  Jackknife(xs, o.num, sum(xs), sumabs2(xs))
end

typealias ObservableSet Dict{AbstractString, Observable}

function add!(obs::ObservableSet, name::AbstractString) 
  if !haskey(obs,name)
    obs[name] = Observable()
  end
  obs
end

