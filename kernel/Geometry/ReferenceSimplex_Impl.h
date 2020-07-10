#ifndef HPGEM_KERNEL_REFERENCESIMPLEX_IMPL_H
#define HPGEM_KERNEL_REFERENCESIMPLEX_IMPL_H

#include "ReferenceSimplex.h"
#include "PointReference.h"

#include "../LinearAlgebra/SmallVector.h"

namespace Geometry {
// Helper functions for the constructor
template <std::size_t DIM>
LinearAlgebra::SmallVector<DIM> createCenter() {
    LinearAlgebra::SmallVector<DIM> result;
    result.set(1.0 / (DIM + 1));
    return result;
}

template <std::size_t DIM>
ReferenceSimplex<DIM>::ReferenceSimplex(const ReferenceGeometryType& geoT,
                                        std::string name)
    : ReferenceGeometry(geoT, name), center_(createCenter<DIM>()) {
    LinearAlgebra::SmallVector<DIM> zero;
    zero.set(0);
    for (std::size_t i = 0; i <= DIM; ++i) {
        this->points_[i].setCoordinates(zero);
        if (i > 0) this->points_[i].setCoordinate(i - 1, 1);
    }
}

template <std::size_t DIM>
LinearAlgebra::SmallVector<DIM + 1>
    ReferenceSimplex<DIM>::baryCentricCoordinates(
        const Geometry::PointReference<DIM>& p) const {
    LinearAlgebra::SmallVector<DIM + 1> result;
    result[0] = 1;
    for (std::size_t i = 0; i < DIM; ++i) {
        result[0] -= p[i];
        result[i + 1] = p[i];
    }
    return result;
}

}  // namespace Geometry

#endif  // HPGEM_KERNEL_REFERENCESIMPLEX_IMPL_H
