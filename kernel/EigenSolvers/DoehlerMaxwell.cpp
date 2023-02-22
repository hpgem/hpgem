#include <numeric>
#include <iostream>

#include <petscsys.h>
#include "DoehlerMaxwell.h"
#include "Logger.h"

namespace hpgem {

namespace EigenSolvers {

DoehlerMaxwellSolver::DoehlerMaxwellSolver() {
  std::cout << "Initialize DoehlerMaxwellSolver!!!!!" << std::endl;
}


}  // namespace EigenSolvers

}  // namespace hpgem
