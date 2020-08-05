#include "PhysicalGeometryBase.h"

namespace hpgem {

namespace Geometry {
using namespace LinearAlgebra;

template <std::size_t DIM>
double computeDiameterDim(const PhysicalGeometryBase& geom) {
    PointPhysical<DIM> p, diam;
    double result = 0;
    for (std::size_t i = 0; i < geom.getNumberOfNodes() - 1; ++i) {
        p = geom.getLocalNodeCoordinates(i);
        for (std::size_t j = i + 1; j < geom.getNumberOfNodes(); ++j) {
            diam = p - geom.getLocalNodeCoordinates(j);
            // Delay taking the square root, as that is a relatively expensive
            // operation.
            result = std::max(result, diam.getCoordinates().l2NormSquared());
        }
    }

    return std::sqrt(result);
}

double PhysicalGeometryBase::computeDiameter() const {
    switch (referenceGeometry_->getDimension()) {
        case 1:
            return computeDiameterDim<1>(*this);
        case 2:
            return computeDiameterDim<2>(*this);
        case 3:
            return computeDiameterDim<3>(*this);
        case 4:
            return computeDiameterDim<4>(*this);
        default:
            logger.assert_always(
                false,
                "Diameter computation not implemented in this dimension");
            return -1;
    }
}
}  // namespace Geometry
}  // namespace hpgem
