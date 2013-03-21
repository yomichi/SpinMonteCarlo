using SpinMonteCarlo

type LocalOperator
  bond :: Int
  t :: Float64
  ty :: Int
  dID :: Int
  uID :: Int
end

function looper(lat::Lattice, T::Float64, Sweeps::Int64, Thermalizations::Int64)
  nsites = num_site(lat)
  nbonds = num_bond(lat)
  beta = 1.0/T
  lambda = 0.5beta*nbonds

  function update()
    currents = reshape(1:nsites, nsites)
  end
  function measure()
  end
end
