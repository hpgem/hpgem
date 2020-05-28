/*
 * File:   VTKElementOrdering.h
 * Author: brinkf
 *
 * Created on November 26, 2014, 12:14 PM
 */

#ifndef VTKELEMENTORDERING_HPP
#define VTKELEMENTORDERING_HPP

#include <array>
#include <typeinfo>

namespace Output {
///\brief given a local node index in the node ordering VTK uses, return the
/// local node index in hpGEM numbering
inline std::size_t tohpGEMOrdering(std::size_t VTKIndex,
                                   const Geometry::ReferenceGeometry* shape) {
    logger.assert_debug(VTKIndex < shape->getNumberOfNodes(),
                        "A % only has % indices", shape->getName(),
                        shape->getNumberOfNodes());
    if (typeid(*shape) == typeid(Geometry::ReferenceSquare)) {
        return std::array<std::size_t, 4>{{0, 1, 3, 2}}[VTKIndex];
    }
    if (typeid(*shape) == typeid(Geometry::ReferenceCube)) {
        return std::array<std::size_t, 8>{{0, 1, 3, 2, 4, 5, 7, 6}}[VTKIndex];
    }
    if (typeid(*shape) == typeid(Geometry::ReferencePyramid)) {
        return std::array<std::size_t, 5>{{3, 4, 2, 1, 0}}[VTKIndex];
    }
    // all other cases use the same numbering
    return VTKIndex;
}

}  // namespace Output

#endif /* VTKELEMENTORDERING_HPP */
