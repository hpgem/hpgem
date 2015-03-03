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


#include "TecplotPhysicalGeometryIterator.hpp"
#include "Geometry/PhysicalGeometry.hpp"
#include "Geometry/PointPhysical.hpp"

#include "Geometry/PhysicalHexahedron.hpp"
#include "Geometry/PhysicalLine.hpp"
#include "Geometry/PhysicalOctachoron.hpp"
#include "Geometry/PhysicalPyramid.hpp"
#include "Geometry/PhysicalQuadrilateral.hpp"
#include "Geometry/PhysicalTetrahedron.hpp"
#include "Geometry/PhysicalTriangle.hpp"
#include "Geometry/PhysicalTriangularPrism.hpp"

#include <typeinfo>

namespace Output
{
    TecplotPhysicalGeometryIterator::TecplotPhysicalGeometryIterator()
    {
        // sort in the sequences in BACKWARD order!
        // see tecplot manual chapter 4.3

        // 4D FE hyper cube: ORDER IS GUESS
        /*
        hypercubeNodes.push_back(14);
        hypercubeNodes.push_back(15);
        hypercubeNodes.push_back(13);
        hypercubeNodes.push_back(12);
        hypercubeNodes.push_back(10);
        hypercubeNodes.push_back(11);
        hypercubeNodes.push_back(9);
        hypercubeNodes.push_back(8);

        hypercubeNodes.push_back(6);
        hypercubeNodes.push_back(7);
        hypercubeNodes.push_back(5);
        hypercubeNodes.push_back(4);
        hypercubeNodes.push_back(2);
        hypercubeNodes.push_back(3);
        hypercubeNodes.push_back(1);
        hypercubeNodes.push_back(0);

        // 3D FE volumes: (all as bricks)
        */
        hexahedronNodes.push_back(6);
        hexahedronNodes.push_back(7);
        hexahedronNodes.push_back(5);
        hexahedronNodes.push_back(4);
        hexahedronNodes.push_back(2);
        hexahedronNodes.push_back(3);
        hexahedronNodes.push_back(1);
        hexahedronNodes.push_back(0);

        triangularPrismNodes.push_back(5);
        triangularPrismNodes.push_back(5);
        triangularPrismNodes.push_back(4);
        triangularPrismNodes.push_back(3);
        triangularPrismNodes.push_back(2);
        triangularPrismNodes.push_back(2);
        triangularPrismNodes.push_back(1);
        triangularPrismNodes.push_back(0);

        pyramidNodes.push_back(0);
        pyramidNodes.push_back(0);
        pyramidNodes.push_back(0);
        pyramidNodes.push_back(0);
        pyramidNodes.push_back(3);
        pyramidNodes.push_back(4);
        pyramidNodes.push_back(2);
        pyramidNodes.push_back(1);

        tetrahedronNodes.push_back(3);
        tetrahedronNodes.push_back(3);
        tetrahedronNodes.push_back(3);
        tetrahedronNodes.push_back(3);
        tetrahedronNodes.push_back(2);
        tetrahedronNodes.push_back(2);
        tetrahedronNodes.push_back(1);
        tetrahedronNodes.push_back(0);

        // 2D FE-surfaces:
        quadrilateralNodes.push_back(2);
        quadrilateralNodes.push_back(3);
        quadrilateralNodes.push_back(1);
        quadrilateralNodes.push_back(0);

        triangleNodes.push_back(2);
        triangleNodes.push_back(1);
        triangleNodes.push_back(0);
        triangleNodes.push_back(0); // double so that we have four points

        // for the FE-QUADRILATERAL element type

        // 1D data - is this supported by tecplot?
        lineNodes.push_back(1);
        lineNodes.push_back(0);
    }

    void TecplotPhysicalGeometryIterator::acceptG(const Geometry::PhysicalGeometry* geo)
    {
        if(typeid(geo)==typeid(Geometry::PhysicalLine*))
        {
            acceptLineGeometry(dynamic_cast<const Geometry::PhysicalLine*>(geo));
        }
        else if(typeid(geo)==typeid(Geometry::PhysicalTriangle*))
        {
            acceptTriangleGeometry(dynamic_cast<const Geometry::PhysicalTriangle*>(geo));
        }
        else if(typeid(geo)==typeid(Geometry::PhysicalQuadrilateral*))
        {
            acceptQuadrilateralGeometry(dynamic_cast<const Geometry::PhysicalQuadrilateral*>(geo));
        }
        else if(typeid(geo)==typeid(Geometry::PhysicalHexahedron*))
        {
            acceptHexahedronGeometry(dynamic_cast<const Geometry::PhysicalHexahedron*>(geo));
        }
        else if(typeid(geo)==typeid(Geometry::PhysicalTetrahedron*))
        {
            acceptTetrahedronGeometry(dynamic_cast<const Geometry::PhysicalTetrahedron*>(geo));
        }
        else if(typeid(geo)==typeid(Geometry::PhysicalTriangularPrism*))
        {
            acceptTriangularPrismGeometry(dynamic_cast<const Geometry::PhysicalTriangularPrism*>(geo));
        }
        else if(typeid(geo)==typeid(Geometry::PhysicalPyramid*))
        {
            acceptPyramidGeometry(dynamic_cast<const Geometry::PhysicalPyramid*>(geo));
        }
        else if(typeid(geo)==typeid(Geometry::PhysicalOctachoron*))
        {
            throw "not implemented";
            //acceptOctachoronGeometry(dynamic_cast<const Geometry::PhysicalOctachoron*>(geo));
        }
    }

    /*
    void TecplotPhysicalGeometryIterator::acceptHyperCubeGeometry(const Geometry::PhysicalHypercube&)
    {
        currentSequencePtr = &hypercubeNodes;
        currentNode = hypercubeNodes.size() - 1;
    }
    */

    
    void TecplotPhysicalGeometryIterator::acceptHexahedronGeometry(const Geometry::PhysicalHexahedron*)
    {
        currentSequencePtr = &hexahedronNodes;
        currentNode = hexahedronNodes.size() - 1;
    }

    void TecplotPhysicalGeometryIterator::acceptTriangularPrismGeometry(const Geometry::PhysicalTriangularPrism*)
    {
        currentSequencePtr = &triangularPrismNodes;
        currentNode = triangularPrismNodes.size() - 1;
    }

    void TecplotPhysicalGeometryIterator::acceptPyramidGeometry(const Geometry::PhysicalPyramid*)
    {
        currentSequencePtr = &pyramidNodes;
        currentNode = pyramidNodes.size() - 1;
    }

    void TecplotPhysicalGeometryIterator::acceptTetrahedronGeometry(const Geometry::PhysicalTetrahedron*)
    {
        currentSequencePtr = &tetrahedronNodes;
        currentNode = tetrahedronNodes.size() - 1;
    }

    void TecplotPhysicalGeometryIterator::acceptQuadrilateralGeometry(const Geometry::PhysicalQuadrilateral*)
    {
        currentSequencePtr = &quadrilateralNodes;
        currentNode = quadrilateralNodes.size() - 1; // inverse transverse
    }

    void TecplotPhysicalGeometryIterator::acceptTriangleGeometry(const Geometry::PhysicalTriangle*)
    {
        currentSequencePtr = &triangleNodes;
        currentNode = triangleNodes.size() - 1;
    }

    void TecplotPhysicalGeometryIterator::acceptLineGeometry(const Geometry::PhysicalLine*)
    {
        currentSequencePtr = &lineNodes;
        currentNode = lineNodes.size() - 1;
    }

    bool TecplotPhysicalGeometryIterator::more() const
    {
        // this test is why InternalIndexType has to be int !!!
        return (currentNode >= 0);
    }

    std::size_t TecplotPhysicalGeometryIterator::getNodeNr()
    {
        return currentSequencePtr->operator[](currentNode--);
    }
}
