using SpinMonteCarlo

function ising_localupdate(lat::Lattice, T::Float64, Sweeps::Int, Thermalization::Int)
  nsites = num_sites(lat)
  nbonds = num_bonds(lat)
  spins = ones(Int,nsites)
  mag = nsites
  ene = -nbonds

  beta = 1.0/T

  function update()
    for site in 1:nsites
      de = 2spins[site] * sum([spins[s] for s in neighbors(lat,site)])
      p = exp(-beta*de)
      if rand()<p
        mag -= 2spins[site]
        ene += de
        spins[site] *= -1
      end
    end
  end

  obs = ObservableSet()
  add(obs,"Magnetization")
  add(obs,"Magnetization^2")
  add(obs,"Magnetization^4")
  add(obs,"Energy")
  add(obs,"Energy^2")

  function measure()
    m = convert(Float64,  mag/nsites)
    en = convert(Float64, ene/nbonds)
    add(obs["Magnetization"], m)
    add(obs["Magnetization^2"], square(m))
    add(obs["Magnetization^4"], square(square(m)))
    add(obs["Energy"], en)
    add(obs["Energy^2"], square(en))
  end

  for mcs in 1:Thermalization
    update()
  end
  for mcs in 1:Sweeps
    update()
    measure()
  end

  return obs
end

L = 16
lat = square_lattice(L)
Sweeps = 8192
Thermalization = Sweeps>>3

for T in 2.0:0.1:2.4
  observables = ising_localupdate(lat, T, Sweeps, Thermalization)
  for obs in observables
    open("$(obs[1])_$L.dat", "a") do os
      write(os, "$T $(mean(obs[2])) $(error(obs[2]))\n")
    end
  end
end

