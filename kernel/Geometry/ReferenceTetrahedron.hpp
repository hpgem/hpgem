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
//

#ifndef ____ReferenceTetrahedron__
#define ____ReferenceTetrahedron__

#include <iostream>

#include "ReferenceGeometry.hpp"
#include <vector>

namespace Geometry
{
    /* The ordering of the vertex and faces in a tetrahedron:
     *
     * 3 o
     *   |\
     *   |  \2
     *   |  / \
     *   |/     \
     * 0 o--------o 1
     *
     */
    class ReferenceTetrahedron: public ReferenceGeometry
    {
    public:
        typedef ReferenceGeometry ReferenceTetrahedronT;

        using ReferenceTetrahedronT::PointReferenceT;
        using ReferenceTetrahedronT::IndexT;
        using ReferenceTetrahedronT::String;
        using ReferenceTetrahedronT::const_iterator;

    public:
        static ReferenceTetrahedron& Instance()
        {
            static ReferenceTetrahedron theInstance;
            return theInstance;
        }

    private:

        ReferenceTetrahedron();

        ReferenceTetrahedron(const ReferenceTetrahedron& copy);
        
    public:

        //! (see ReferenceGeometry.hpp)
        bool            isInternalPoint(const PointReferenceT& point) const;
        
        //! (see ReferenceGeometry.hpp)
        void            getCenter(PointReferenceT& point) const;
        
        //! (see ReferenceGeometry.hpp)
        void            getNode(const IndexT& i, PointReferenceT& point) const;
        
        //! (see ReferenceGeometry.hpp)
        String          getName() const {return "ReferenceTetrahedron";}
        
        //! Given a face index, and an index of the node position relative to the face,
        //! return the local index of the node.
        std::size_t             getLocalNodeIndex(std::size_t face, std::size_t node) const 
        {
            return localNodeIndexes_[face][node];
        }

        /// Output routine.
        friend std::ostream& operator<<(std::ostream& os, const ReferenceTetrahedron& point);
        
        // ================================== Codimension 0 ========================================

        //! (see MappingCodimensions.hpp)
        std::size_t                                      getCodim0MappingIndex(const ListOfIndexesT&, const ListOfIndexesT&) const;

        //! (see MappingCodimensions.hpp)
        const MappingReferenceToReference* getCodim0MappingPtr(const IndexT) const;

        using MappingCodimensions::getCodim0MappingPtr;

        // ================================== Codimension 1 ========================================

        //! (see MappingCodimensions.hpp)
        std::size_t                             getNrOfCodim1Entities() const {return 4;}

        //! (see MappingCodimensions.hpp)
        void                                     getCodim1EntityLocalIndices(const IndexT, ListOfIndexesT& faceNodesLocal) const;

        //! (see MappingCodimensions.hpp)
        const MappingReferenceToReference* getCodim1MappingPtr(const IndexT) const;

        //! (see MappingCodimensions.hpp)
        const ReferenceGeometry*              getCodim1ReferenceGeometry(const IndexT) const;

        // ================================== Codimension 2 ========================================

        //! (see MappingCodimensions.hpp)
        std::size_t                             getNrOfCodim2Entities() const {return 6;};

        //! (see MappingCodimensions.hpp)
        void                                     getCodim2EntityLocalIndices(const IndexT, ListOfIndexesT& faceNodesLocal) const;

        //! (see MappingCodimensions.hpp)
        const MappingReferenceToReference* getCodim2MappingPtr(const IndexT) const;

        //! (see MappingCodimensions.hpp)
        const ReferenceGeometry*              getCodim2ReferenceGeometry(const IndexT) const;

        // ================================== Codimension 3 ========================================

        //! (see MappingCodimensions.hpp)
        std::size_t                             getNrOfCodim3Entities() const {return 4;};

        //! (see MappingCodimensions.hpp)
        void                                     getCodim3EntityLocalIndices(const IndexT, ListOfIndexesT& faceNodesLocal) const;
            
        // =============================== Refinement mappings =====================================
        
        //! Transform a reference point using refinement mapping
        void refinementTransform(int refineType, std::size_t subElementIdx,
                      const PointReferenceT& p, PointReferenceT& pMap) const {}

        //! Transformation matrix of this refinement when located on the LEFT side
        void getRefinementMappingMatrixL(int refineType, std::size_t subElementIdx,
                    LinearAlgebra::Matrix& Q) const {}

        //! Transformation matrix of this refinement when located on the RIGHT side
        void getRefinementMappingMatrixR(int refineType, std::size_t subElementIdx,
                    LinearAlgebra::Matrix& Q) const {}

        //! Refinement mapping on codim1 for a given refinement on codim0
        //! Note: this should also applied on other dimensions
        void getCodim1RefinementMappingMatrixL(int refineType, DimT subElementIdx,
                                DimT faLocalIndex, LinearAlgebra::Matrix& Q) const {}

        //! Refinement mapping on codim1 for a given refinement on codim0
        //! Note: this should also applied on other dimensions
        void getCodim1RefinementMappingMatrixR(int refineType, DimT subElementIdx,
                                DimT faLocalIndex, LinearAlgebra::Matrix& Q) const {}

    private:
        
        //! Local node indexes contains the numbering of the vertex of the shape, ordered by faces.
        //! See top comment for the corresponding numbering.
        static std::size_t                                          localNodeIndexes_[4][3];
        static std::size_t                                          localNodesOnEdge_[6][2];
        
        //! Codimension 1 mappings, from a square to a tetrahedron face. TODO: Where is this used? clarify here.
        const MappingReferenceToReference*            mappingsTriangleToTetrahedron_[4];
        const MappingReferenceToReference*            mappingsTetrahedronToTetrahedron_[1];

        //! Pointer to the Codimension 1 reference geometry.
        ReferenceGeometry* const                         referenceGeometryCodim1Ptr_;
        ReferenceGeometry* const                         referenceGeometryCodim2Ptr_;
        
        //! List of valid quadrature rules for this reference geometry
        std::vector<QuadratureRules::GaussQuadratureRule*> lstGaussQuadratureRules_;
    };
};

#endif /* defined(____ReferenceTetrahedron__) */
