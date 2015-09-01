#include "SW.hpp"
#include <alps/parapack/parapack.h>

PARAPACK_SET_VERSION("Swendsen-Wang algorithm");
PARAPACK_REGISTER_WORKER(potts::Worker, "cluster");
PARAPACK_REGISTER_EVALUATOR(potts::Evaluator, "cluster");

int main(int argc, char** argv) { return alps::parapack::start(argc, argv); }
