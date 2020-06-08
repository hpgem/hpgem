#ifndef HPGEM_DIVDGMAXEVTESTCASE_H
#define HPGEM_DIVDGMAXEVTESTCASE_H

#include <memory>

#include "Algorithms/DivDGMaxEigenvalue.h"
#include "Utils/Verification/RunnableEVTestCase.h"
#include "Utils/Verification/EVTestCase.h"

namespace DGMax {

template <std::size_t DIM>
class DivDGMaxEVTestCase : public RunnableEVTestCase<DIM> {
   public:
    DivDGMaxEVTestCase(EVTestCase<DIM> testCase,
                       std::vector<std::string> meshFileNames, double tolerance,
                       std::size_t order,
                       typename DivDGMaxDiscretization<DIM>::Stab stab,
                       EVConvergenceResult* expected)
        : testCase_(testCase),
          meshFileNames_(std::move(meshFileNames)),
          tolerance_(tolerance),
          order_(order),
          stab_(stab),
          expected_(expected){};

   protected:
    std::size_t getNumberOfLevels() const override {
        return meshFileNames_.size();
    }
    double getTolerance() const override { return tolerance_; }
    const EVConvergenceResult* getExpected() const override {
        return expected_;
    }
    /// Run the actual algorithm on a single level.
    std::unique_ptr<AbstractEigenvalueResult<DIM>> runInternal(
        std::size_t level) override;

   private:
    EVTestCase<DIM> testCase_;
    std::vector<std::string> meshFileNames_;
    double tolerance_;
    std::size_t order_;
    typename DivDGMaxDiscretization<DIM>::Stab stab_;
    // TODO: This is usually static data
    EVConvergenceResult* expected_;
};

}  // namespace DGMax

#endif  // HPGEM_DIVDGMAXEVTESTCASE_H
