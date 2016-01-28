include("SpinMonteCarlo.jl")
using SpinMonteCarlo

square(x) = x*x

type Ising
  lat :: Lattice
  spins :: Vector{Int}
  T :: Float64
  beta :: Float64
  h :: Float64
  mag :: Float64
  ene :: Float64

  function Ising(L, T, h, initial::Symbol=:up)
    @show T, h
    lat = square_lattice(L)
    nsites = num_sites(lat)
    nbonds = num_bonds(lat)
    spins = (initial == :up ? ones(Int,nsites)
                           : (initial == :down ? -ones(Int, nsites)
                                               : rand([1,-1], nsites)))

    mag = sum(spins)
    ene = -mag*h
    for site in 1:nsites
      s = spins[site]
      for site2 in neighbors(lat, site)
        s2 = spins[site2]
        ene -= (s == s2 ? 1 : -1)
      end
    end
    new(lat, spins, T, 1.0/T, h, mag, ene)
  end
end

function update!(model::Ising)
  nsites = num_sites(model.lat)
  for site in 1:nsites
    de = 2model.spins[site] * (sum([model.spins[s] for s in neighbors(model.lat,site)]) + model.h)
    p = exp(-model.beta*de)
    if rand()<p
      model.mag -= 2model.spins[site]
      model.ene += de
      model.spins[site] *= -1
    end
  end
end

function measure!(obs, model)
  nsites = num_sites(model.lat)
  nbonds = num_bonds(model.lat)
  m = model.mag/nsites
  en = model.ene/nbonds
  add!(obs["Magnetization"], m)
  add!(obs["Magnetization^2"], square(m))
  add!(obs["Magnetization^4"], square(square(m)))
  add!(obs["Energy"], en)
  add!(obs["Energy^2"], square(en))
end

function ising_localupdate(L::Int, T::Float64, h::Float64, Sweeps::Int, Thermalization::Int; initial::Symbol = :up)

  model = Ising(L, T, h, initial)

  obs = ObservableSet()
  add!(obs,"Magnetization")
  add!(obs,"Magnetization^2")
  add!(obs,"Magnetization^4")
  add!(obs,"Energy")
  add!(obs,"Energy^2")

  for mcs in 1:Thermalization
    update!(model)
  end
  for mcs in 1:Sweeps
    update!(model)
    measure!(obs, model)
  end

  return obs
end

