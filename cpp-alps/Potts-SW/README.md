# Potts-SW
Swendsen-Wang algorithm for classical Ising / Potts model implemented with ALPS/parapack library.

# Requirement
- ALPS v2

# Build
Type `ALPS_HOME=/where/ALPS/is/built ./make.sh`.
Executable file `potts` will be built in `build` directory.

# Run

    # copy template files
    cp -r template 2D && cd 2D
    
    # edit input parameters
    vi init_param.py
    
    # generate input files
    python init_param.py
    
    # run
    mpiexec -np 4 ../build/potts --mpi params.in.xml
    
    # extract results
    python extract.py
    
    # several output files have been generated
    ls *.dat

# Input parameters
- `MODEL`
    - "Potts" or "Ising" model
    - "Potts"
        - `E = -J \sum_{<ij>} \delta_{S_i, S_j}`
            - `S_i = 1, 2, \dots, q` is a spin at site i
            - `\delta` is Kroneker delta
            - `\sum_{<ij>}` runs all nearest neighbor pairs
    - "Ising"
        - `E = -J \sum_{<ij>} S_i S_j`
        - `S_i = +/- 1`
- `q`
    - The number of states of one spin
- `J`
    - Coupling constant (only ferromagnet system (J>0) can be dealt with)
    - By default, `J = 1`
- `LATTICE`
    - Lattice name ('chain lattice', 'square lattice', 'simple cubic lattice', ...)
    - please consult [ALPS web page](http://alps.comp-phys.org/mediawiki/index.php/Tutorials:LatticeHOWTO).
- `L`, `W`, `H`
    - Lattice length, width, and height.
    - `W` and `H` are equal to `L` and `W` by default.
- `T`
    - Temperature
- `SWEEPS`
    - Monte Carlo steps
    - By default, `SWEEPS = 65536`
- `THERMALIZATION`
    - Thermalization steps
    - By default, `THERMALIZATION = SWEEPS >> 3`

# Observables
Except for `|Magnetization|` and `Connected Susceptibility`,
improved estimator is implemented.

- squared magnetization per site `<M^2>`
    - Observable name :  `Magnetization^2`
    - result file (`extract.py`) : `mag2-L*.dat`
- biquadratic magnetization per site `<M^4>`
    - Observable name : `Magnetization^4`
    - result file (`extract.py`) : `mag4-L*.dat`
- absolute magnetization per site `<|M|>`
    - Observable name : `|Magnetization|`
    - result file (`extract.py`) : `amag-L*.dat`
- Binder ratio defined as `<M^4>/<M^2>^2`
    - Observable name : `Binder Ratio`
    - result file (`extract.py`) : `binder-L*.dat`
- Magnetic susceptibility `N <M^2> / T`
    - Obaservable name : `Susceptibility`
    - This _survives_ in ordered state
    - result file (`extract.py`) : `sus-L*.dat`
- Magnetic susceptibility `N (<M^2> - <|M|>^2 ) /T`
    - Observable name : `Connected Susceptibility`
    - This _vanishes_ in ordered state
    - result file (`extract.py`) : `csus-L*.dat`
- total energy (_not_ per site or bond)
    - Observable name : `Energy`
    - result file (`extract.py`) : `ene-L*.dat`
- specific heat
    - Observable name : `Specific Heat`
    - result file (`extract.py`) : `spec-L*.dat`

# License
This program is distributed under the Boost Software License.
See [LICENSE_1_0.txt](LICENSE_1_0.txt) for detail.
