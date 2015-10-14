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
#include "PhysicalGeometryAcceptor.h"

namespace Geometry
{
    class PhysicalGeometryBase;
}

namespace Output
{
    
    /// \todo: Implement other geometries.
    class TecplotPhysicalGeometryIterator : public PhysicalGeometryAcceptor
    {
        
    public:
        
        /// return reference to singleton instance:
        static TecplotPhysicalGeometryIterator& Instance()
        {
            static TecplotPhysicalGeometryIterator theInstance;
            return theInstance;
        }
        
        void acceptG(const Geometry::PhysicalGeometryBase* geo);
        
        /// \brief Choose node sequence for the hexahedron.
        void acceptHexahedronGeometry(const Geometry::PhysicalGeometryBase*) override final;

        /// \brief Choose node sequence for the prism.
        void acceptTriangularPrismGeometry(const Geometry::PhysicalGeometryBase*) override final;

        /// \brief Choose node sequence for the pyramid.
        void acceptPyramidGeometry(const Geometry::PhysicalGeometryBase*) override final;

        /// \brief Choose node sequence for the tetrahedron.
        void acceptTetrahedronGeometry(const Geometry::PhysicalGeometryBase*) override final;

        /// \brief Choose node sequence for quadrilateral.
        void acceptQuadrilateralGeometry(const Geometry::PhysicalGeometryBase*) override final;

        /// \brief Choose node sequence for triangle.
        void acceptTriangleGeometry(const Geometry::PhysicalGeometryBase*) override final;

        /// \brief Choose node sequence for line.
        void acceptLineGeometry(const Geometry::PhysicalGeometryBase*) override final;

        /// \brief Check whether all nodes of current shape have been used.
        bool more() const;

        /// \brief Return the current node number.
        std::size_t getNodeNumber();

    private:
        TecplotPhysicalGeometryIterator();
        TecplotPhysicalGeometryIterator(const TecplotPhysicalGeometryIterator&) = delete;
        TecplotPhysicalGeometryIterator& operator=(const TecplotPhysicalGeometryIterator&) = delete;
        
        using InternalIndexType = int;
        using VectorOfNodeIndexes = std::vector<std::size_t>;
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
