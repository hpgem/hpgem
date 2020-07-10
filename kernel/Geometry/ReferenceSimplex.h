#ifndef HPGEM_KERNEL_REFERENCESIMPLEX_H
#define HPGEM_KERNEL_REFERENCESIMPLEX_H

#include "../Logger.h"

#include "ReferenceGeometry.h"

namespace hpgem {

namespace LinearAlgebra {
template <std::size_t DIM>
class SmallVector;
}

namespace Geometry {
/// Implementation of a reference simplex (triangle, tetrahedron)
///
/// The standard reference simplex is formed by the origin and each of the
/// unit vectors, numbered in increasing order of the coordinate (e.g. O,X,Y)
///
/// Note that the reference line is based on the hypercube.
/// \tparam DIM The dimension of the simplex.
template <std::size_t DIM>
class ReferenceSimplex : public ReferenceGeometry {
   public:
    ReferenceSimplex(const ReferenceGeometryType& geoT, std::string name);

    /// Compute the bary centric coordinates of a point
    /// \param p The point
    /// \return The barycentric coordinates
    LinearAlgebra::SmallVector<DIM + 1> baryCentricCoordinates(
        const PointReference<DIM>& p) const;

    template <std::size_t DIM1>
    friend std::ostream& operator<<(std::ostream& os,
                                    const ReferenceSimplex<DIM1>& simplex);

    const PointReferenceBase& getCenter() const final { return center_; }

    std::size_t getNumberOfNodes() const final { return DIM + 1; }

    const PointReferenceBase& getReferenceNodeCoordinate(
        const std::size_t& localIndex) const final {
        logger.assert_debug(localIndex >= 0 && localIndex < getNumberOfNodes(),
                            "Local index % should be in the range [0, %]",
                            localIndex, getNumberOfNodes());
        return points_[localIndex];
    }

    const PointReference<DIM>& getPoints() { return points_; }

   private:
    const PointReference<DIM> center_;
    // Leave non const to allow actual initialization
    PointReference<DIM> points_[DIM + 1];
};

template <std::size_t DIM>
std::ostream& operator<<(std::ostream& os,
                         const ReferenceSimplex<DIM>& simplex) {
    os << simplex.name_ << " =( ";
    for (std::size_t i = 0; i < simplex.getNumberOfNodes(); ++i) {
        os << simplex.points_[i] << "\t";
    }
    os << ")" << std::endl;
    return os;
}
}  // namespace Geometry

}  // namespace hpgem

#endif  // HPGEM_KERNEL_REFERENCESIMPLEX_H
