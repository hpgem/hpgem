#ifndef HPGEM_ABSTRACTEIGENVALUESOLVER_H
#define HPGEM_ABSTRACTEIGENVALUESOLVER_H

#include <memory>

#include "AbstractEigenvalueResult.h"
#include "EigenvalueProblem.h"

/// Solver for the EigenvalueProblem
template <std::size_t DIM>
class AbstractEigenvalueSolver {
   public:
    virtual ~AbstractEigenvalueSolver() = default;
    /// Solve the EigenvalueProblem
    virtual std::unique_ptr<AbstractEigenvalueResult<DIM>> solve(
        const EigenvalueProblem<DIM>& input) = 0;
};

#endif  // HPGEM_ABSTRACTEIGENVALUESOLVER_H
