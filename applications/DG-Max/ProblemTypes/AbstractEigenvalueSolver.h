#ifndef HPGEM_ABSTRACTEIGENVALUESOLVER_H
#define HPGEM_ABSTRACTEIGENVALUESOLVER_H

#include <memory>

#include "EigenValueProblem.h"
#include "AbstractEigenvalueResult.h"

template <std::size_t DIM>
class AbstractEigenvalueSolver {
   public:
    virtual ~AbstractEigenvalueSolver() = default;
    /// Solve an eigenvalue problem.
    virtual std::unique_ptr<AbstractEigenvalueResult<DIM>> solve(
        const EigenValueProblem<DIM>& input) = 0;
};

#endif  // HPGEM_ABSTRACTEIGENVALUESOLVER_H
