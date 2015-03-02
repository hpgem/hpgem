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

#ifndef TECPLOTPHYSICALGEOMETRYITERATOR_HH
#define TECPLOTPHYSICALGEOMETRYITERATOR_HH

#include <vector>
#include "PhysicalGeometryAcceptor.hpp"

namespace Geometry {
	class PhysicalGeometry;
}

namespace Output
{

    /// TODO: Implement other geometries. For that we need the physical geometries.
    class TecplotPhysicalGeometryIterator: public PhysicalGeometryAcceptor
    {

    public:

        /// return reference to singleton instance:
        static TecplotPhysicalGeometryIterator& Instance()
        {
            static TecplotPhysicalGeometryIterator theInstance;
            return theInstance;
        }

		void acceptG(const Geometry::PhysicalGeometry* geo);

        /// \brief Choose node sequence for the hypercube.
        //virtual void acceptHyperCubeGeometry(const Geometry::PhysicalHypercube&);

        /// \brief Choose node sequence for the hexahedron.
        virtual void acceptHexahedronGeometry(const Geometry::PhysicalHexahedron*);

        /// \brief Choose node sequence for the prism.
        virtual void acceptTriangularPrismGeometry(const Geometry::PhysicalTriangularPrism*);

        /// \brief Choose node sequence for the pyramid.
        virtual void acceptPyramidGeometry(const Geometry::PhysicalPyramid*);

        /// \brief Choose node sequence for the tetrahedron.
        virtual void acceptTetrahedronGeometry(const Geometry::PhysicalTetrahedron*);

        /// \brief Choose node sequence for quadrilateral.
        virtual void acceptQuadrilateralGeometry(const Geometry::PhysicalQuadrilateral*);

        /// \brief Choose node sequence for triangle.
        virtual void acceptTriangleGeometry(const Geometry::PhysicalTriangle*);

        /// \brief Choose node sequence for line.
        virtual void acceptLineGeometry(const Geometry::PhysicalLine*);

        /// \brief Check whether all nodes of current shape have been used.
        bool more() const;

        /// \brief Return the current node number.
        std::size_t getNodeNr();

    private:
        TecplotPhysicalGeometryIterator();
        TecplotPhysicalGeometryIterator(const TecplotPhysicalGeometryIterator&)=delete;
        TecplotPhysicalGeometryIterator& operator=(const TecplotPhysicalGeometryIterator&)=delete;
        virtual ~TecplotPhysicalGeometryIterator() { }

        typedef int InternalIndexType;
        typedef std::vector<std::size_t> VectorOfNodeIndexes;
        VectorOfNodeIndexes hypercubeNodes;
        VectorOfNodeIndexes hexahedronNodes;
        VectorOfNodeIndexes cubeNodes;
        VectorOfNodeIndexes triangularPrismNodes;
        VectorOfNodeIndexes pyramidNodes;
        VectorOfNodeIndexes tetrahedronNodes;
        VectorOfNodeIndexes quadrilateralNodes;
        VectorOfNodeIndexes triangleNodes;
        VectorOfNodeIndexes lineNodes;

        VectorOfNodeIndexes* currentSequencePtr;
        InternalIndexType currentNode;
    };

		}
#endif
