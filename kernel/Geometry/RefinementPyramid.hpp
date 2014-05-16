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

#ifndef RefinementPyramid_hpp
#define RefinementPyramid_hpp

#include "Geometry/RefinementGeometry.hpp"

namespace Geometry
{
    class RefinementPyramid : public Geometry::RefinementGeometry
    {
    public:
        typedef unsigned int                    DimT;
        typedef Point                        PointT;
        typedef PointPhysical                PointPhysicalT;
        typedef PointReference               PointReferenceT;
        typedef PhysicalGeometry             PhysicalGeometryT;
        typedef ReferenceGeometry            ReferenceGeometryT;
        typedef std::vector<PointPhysicalT>     VectorOfPointPhysicalsT;
        typedef std::vector<unsigned int>       VectorOfIndicesT;

        /// Constructors.
        RefinementPyramid(const ReferenceGeometryT* const referenceGeometry,
                          const PhysicalGeometryT* const physicalGeometry)
                        : referenceGeometry_(referenceGeometry), physicalGeometry_(physicalGeometry)
        {
          std::cout << "RefinementPyramid(referenceGeometry, physicalGeometry)\n";
        }

        RefinementPyramid(const RefinementPyramid& other)
                        : referenceGeometry_(other.referenceGeometry_), 
                          physicalGeometry_(other.physicalGeometry_)
        {
          std::cout << "RefinementPyramid(other)\n";
        }

        virtual ~RefinementPyramid()
        {}

        /// For debugging and checkpointing: a human-readable name.
        virtual std::string getName() const
        { 
          return "RefinementPyramid";
        }

        //---------------------- Refinement definitions -----------------------------------------
        
        /// Number of new nodes due to the refinement that should be added to
        /// the vector of localNodeIndices
        virtual DimT nrOfNewNodes(int refineType) const { return 0; }
        
        /// New physical nodes due to refinement to the nodes vector
        /// \param newPoints On input, this vector contains the element's physical nodes.  
        /// On exit, the new physical nodes are added.
        virtual void getAllNodes(int refineType, VectorOfPointPhysicalsT& nodes) const {}
        
        /// Number of sub-elements due to the refinement
        virtual DimT nrOfSubElements(int refineType) const { return 0; }

        /// Assembly nodes for sub-element
        virtual void subElementLocalNodeIndices(int refineType, DimT iSubElement, VectorOfIndicesT& LocalNodeIdx) const {}
        
        /// Local indices pairs of sub-elements connected by a sub-Internal Face
        virtual void adjacentSubElementsPairs(int refineType,
                        VectorOfIndicesT& elemIdx1, VectorOfIndicesT& localFaceIdx1,
                        VectorOfIndicesT& elemIdx2, VectorOfIndicesT& localFaceIdx2) const {}

        /// Number of sub-elements on a parent's face.
        virtual DimT nrOfSubElementsOnFace(int refineType, DimT faLocalIndex) const { return 0; }

        /// Get sub-elements' local index on a parent's face.
        virtual void subElementsOnFace(int refineType, DimT faLocalIndex, VectorOfIndicesT& localSubElemIdx) const {}
        
        /// Get sub-face's local face number of on a parent's face.
        virtual DimT getLocalSubFaceNr(int refineType, DimT localFaceNr, DimT subElementIdx) const { return 0; }

    private:
        RefinementPyramid() : referenceGeometry_(NULL), physicalGeometry_(NULL) {}

        /// The physicalGeometry object contains pointers to the actual physical nodes, and
        /// a container of global node indexes.
        const PhysicalGeometryT*        physicalGeometry_;

        /// The corresponding referenceGeometry object
        const ReferenceGeometryT* const referenceGeometry_;
    };
}
#endif
