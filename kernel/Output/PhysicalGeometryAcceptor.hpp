/*
 This file forms part of hpGEM. This package has been developed over a number of years by various people at the University of Twente and a full list of contributors can be found at
 http://hpgem.org/about-the-code/team
 
 This code is distributed using BSD 3-Clause License. A copy of which can found below.
 
 
 Copyright (c) 2014, University of Twente
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 
 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
