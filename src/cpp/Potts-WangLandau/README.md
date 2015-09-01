# Potts-WL
WangLandau algorithm for classical Potts model on square lattices implemented in C++.

# Requirement
- Boost ProgramOption Library

# Build
Type `./make.sh`.
Executable file `potts` will be built in `build` directory.

# Run

    # show help message
    ./potts --help

    options:
      -h [ --help ]                         show this message
      -q [ --q ] arg (=2)                   number of state of a spin
      -L [ --L ] arg (=10)                  length of lattice
      -f [ --flatness ] arg (=0.80000000000000004)
                                            target flatness of histogram
      -F [ --final ] arg (=1e-08)           minimum updating factor
      -i [ --interval ] arg (=32)           interval between checks for flatness
      -v [ --verbose ]                      show histogram verbose info

    # q=2 Potts model (Ising model) 
    ./potts -q 2

    # L=16 square lattice
    ./potts -L 16

    # show verbose informations such as smallest entry and mean of histogram
    ./potts -v

    # after each WL iteration, log of DoS and histogram are dumped into `hist-#.dat`
    ls *.dat
    
# License
This program is distributed under the Boost Software License.
See [LICENSE_1_0.txt](LICENSE_1_0.txt) for detail.
