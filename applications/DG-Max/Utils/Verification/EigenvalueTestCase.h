#ifndef HPGEM_EIGENVALUETESTCASE_H
#define HPGEM_EIGENVALUETESTCASE_H

#include <utility>
#include <vector>

#include "LinearAlgebra/SmallVector.h"
#include "EigenvalueResult.h"

namespace DGMax {

template <std::size_t DIM>
class EigenvalueTestCase {
   public:
    EigenvalueTestCase(const LinearAlgebra::SmallVector<DIM>& kpoint,
                       std::vector<std::string>  meshes,
                       size_t structureId, size_t numberOfEigenvalues)
        : kpoint_(kpoint),
          structureId_(structureId),
          numberOfEigenvalues_(numberOfEigenvalues) {}

    const LinearAlgebra::SmallVector<DIM>& getKPoint() const {
        return kpoint_;
    };

    std::size_t getStructureId() const {
        return structureId_;
    }

    std::size_t getNumberOfEigenvalues() const {
        return numberOfEigenvalues_;
    }

   private:
    LinearAlgebra::SmallVector<DIM> kpoint_;
    std::size_t structureId_;
    std::size_t numberOfEigenvalues_;
};

}  // namespace DGMax

#endif  // HPGEM_EIGENVALUETESTCASE_H
