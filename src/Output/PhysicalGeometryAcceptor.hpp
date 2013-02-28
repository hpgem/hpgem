/*
 * PhysicalGeometryAcceptor.hpp
 *
 *  Created on: Feb 16, 2013
 *      Author: nicorivas
 */

#ifndef PHYSICALGEOMETRYACCEPTOR_HPP_
#define PHYSICALGEOMETRYACCEPTOR_HPP_

#include "Geometry/PhysicalTetrahedron.hpp"
#include "Geometry/PhysicalPyramid.hpp"
#include "Geometry/PhysicalTriangularPrism.hpp"
#include "Geometry/PhysicalHexahedron.hpp"
#include "Geometry/PhysicalQuadrilateral.hpp"
#include "Geometry/PhysicalTriangle.hpp"
#include "Geometry/PhysicalLine.hpp"
//class Geometry::PhysicalQuadrilateral;

namespace Output
{
    /// TODO: Implement other geometries.
    class PhysicalGeometryAcceptor
    {

    public:

//        virtual void acceptHyperCubeGeometry(const Geometry::PhysicalHypercube&) = 0;
        virtual void acceptTetrahedronGeometry(const Geometry::PhysicalTetrahedron*) = 0;
        virtual void acceptPyramidGeometry(const Geometry::PhysicalPyramid*) = 0;
        virtual void acceptTriangularPrismGeometry(const Geometry::PhysicalTriangularPrism*) = 0;
        virtual void acceptHexahedronGeometry(const Geometry::PhysicalHexahedron*) = 0;
        virtual void acceptQuadrilateralGeometry(const Geometry::PhysicalQuadrilateral*) = 0;
        virtual void acceptTriangleGeometry(const Geometry::PhysicalTriangle*) = 0;
        virtual void acceptLineGeometry(const Geometry::PhysicalLine*) = 0;
        virtual ~PhysicalGeometryAcceptor() {};
    };
}
#endif
