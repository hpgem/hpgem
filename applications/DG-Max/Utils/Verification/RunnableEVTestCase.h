#ifndef HPGEM_RUNNABLEEVTESTCASE_H
#define HPGEM_RUNNABLEEVTESTCASE_H

#include <memory>

#include "ProblemTypes/AbstractEigenvalueResult.h"
#include "EVTestCase.h"

namespace DGMax {

template <std::size_t DIM>
class RunnableEVTestCase {
   public:
    virtual ~RunnableEVTestCase() = default;

    EVConvergenceResult runWithResults(bool breakOnError);

   protected:
    /// Number of mesh levels
    virtual std::size_t getNumberOfLevels() const = 0;
    /// Absolute tolerance in checking the results
    virtual double getTolerance() const = 0;
    /// Expected result (null if none)
    virtual const EVConvergenceResult* getExpected() const = 0;
    /// Run the actual algorithm on a single level.
    virtual std::unique_ptr<AbstractEigenvalueResult<DIM>> runInternal(
        std::size_t level) = 0;

   private:
    bool compareWithExpected(std::size_t level,
                             const AbstractEigenvalueResult<DIM>& result) const;
};

}  // namespace DGMax

#endif  // HPGEM_RUNNABLEEVTESTCASE_H
