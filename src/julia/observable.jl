export Observable,ObservableSet
export add,reset,mean,variance,error,name

type Observable
  name :: String
  num :: Int
  sum :: Float64
  sum2 :: Float64
end

Observable() = Observable("")
Observable(name::String) = Observable(name, 0, 0.0, 0.0)
function add(o::Observable, val::Float64)
  o.num += 1
  o.sum += val
  o.sum2 += val*val
  return o
end
function reset(o::Observable)
  o.num = 0
  o.sum = o.sum2 = 0.0
  return o
end
mean(o::Observable) = o.sum / o.num
variance(o::Observable) = o.num*(o.sum2/o.num - square(mean(o)))/(o.num-1)
error(o::Observable) = sqrt(variance(o)/o.num)
name(o::Observable) = obs.name

typealias ObservableSet Dict{String, Observable}

function add(obs::ObservableSet, name::String) 
  if !has(obs,name)
    obs[name] = Observable(name)
  end
  obs
end

