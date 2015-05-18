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

#ifndef REFERENCECUBE_HH
#define REFERENCECUBE_HH

#include "ReferenceGeometry.h"
#include "PointReference.h"
#include <vector>

namespace Geometry
{
    /*! The ordering of the vertex and faces in a cube:
     *
     *     6o---------o7
     *     /|        /|
     *    / |       / |
     *  4o---------o5 |
     *   | 2o------|--o3
     *   | /       | /
     *   |/        |/
     *  0o---------o1
     *
     */
    class ReferenceCube : public ReferenceGeometry
    {
    public:
        using typename ReferenceGeometry::VectorOfReferencePointsT;

        using ReferenceGeometry::String;
        
        static ReferenceCube& Instance()
        {
            static ReferenceCube theInstance;
            return theInstance;
        }

        ReferenceCube(const ReferenceCube& other) = delete;
        
        //! (see ReferenceGeometry.h)
        bool isInternalPoint(const PointReference& p) const override final;
        
        //! Output routine
        friend std::ostream& operator<<(std::ostream& os, const ReferenceCube& cube);

        // ================================== Codimension 0 ========================================
        
        //! (see MappingCodimensions.h)
        std::size_t getCodim0MappingIndex(const ListOfIndexesT&, const ListOfIndexesT&) const override final;

        //! (see MappingCodimensions.h)
        const MappingReferenceToReference* getCodim0MappingPtr(const std::size_t) const override final;

        using MappingCodimensions::getCodim0MappingPtr;

        // ================================== Codimension 1 ========================================
        
        //! (see MappingCodimensions.h)
        std::size_t getNrOfCodim1Entities() const override final
        {
            return 6;
        }
        
        //! (see MappingCodimensions.h)
        std::vector<std::size_t> getCodim1EntityLocalIndices(const std::size_t) const override final;

        //! (see MappingCodimensions.h)
        const MappingReferenceToReference* getCodim1MappingPtr(const std::size_t) const override final;

        //! (see MappingCodimensions.h)
        const ReferenceGeometry* getCodim1ReferenceGeometry(const std::size_t) const override final;

        // ================================== Codimension 2 ========================================
        
        //! (see MappingCodimensions.h)
        std::size_t getNrOfCodim2Entities() const override final
        {
            return 12;
        }
        
        //! (see MappingCodimensions.h)
        std::vector<std::size_t> getCodim2EntityLocalIndices(const std::size_t) const override final;

        //! (see MappingCodimensions.h)
        const MappingReferenceToReference* getCodim2MappingPtr(const std::size_t) const override final;

        //! (see MappingCodimensions.h)
        const ReferenceGeometry* getCodim2ReferenceGeometry(const std::size_t) const override final;

        // ================================== Codimension 3 ========================================
        
        //! (see MappingCodimensions.h)
        std::size_t getNrOfCodim3Entities() const override final
        {
            return 8;
        }

        //! (see MappingCodimensions.h)
        std::vector<std::size_t> getCodim3EntityLocalIndices(const std::size_t node) const override final
        {
            return std::vector<std::size_t>(1, node);
        }
        
        // =============================== Refinement mappings =====================================
        
        //! Transform a reference point using refinement mapping
        void refinementTransform(int refineType, std::size_t subElementIdx, const PointReference& p, PointReference& pMap) const;

        //! Transformation matrix of this refinement when located on the LEFT side
        void getRefinementMappingMatrixL(int refineType, std::size_t subElementIdx, LinearAlgebra::MiddleSizeMatrix& Q) const;

        //! Transformation matrix of this refinement when located on the RIGHT side
        void getRefinementMappingMatrixR(int refineType, std::size_t subElementIdx, LinearAlgebra::MiddleSizeMatrix& Q) const;

        //! Refinement mapping on codim1 for a given refinement on codim0
        void getCodim1RefinementMappingMatrixL(int refineType, std::size_t subElementIdx, std::size_t faLocalIndex, LinearAlgebra::MiddleSizeMatrix& Q) const;

        //! Refinement mapping on codim1 for a given refinement on codim0
        void getCodim1RefinementMappingMatrixR(int refineType, std::size_t subElementIdx, std::size_t faLocalIndex, LinearAlgebra::MiddleSizeMatrix& Q) const;

    private:
        
        ReferenceCube();
        
        //! Local node indexes contains the numbering of the vertex of the shape, ordered by faces.
        //! See top comment for the corresponding numbering.
        static std::size_t localNodeIndexes_[6][4];

        //! The nodes on edge contains the local index of the two nodes in every edge.
        static std::size_t localNodesOnEdge_[12][2]; //!< 12 edges with 2 nodes
        
        //! Codimension 1 mappings, from a line to a square. TODO: Where is this used? clarify here.
        const MappingReferenceToReference* mappingsSquareToCube_[6];

        //! Codimension 0 mappings, from a square to a square. TODO: Where is this used? clarifiy here.
        const MappingReferenceToReference* mappingsCubeToCube_[8];

        //! Pointer to the Codimension 1 reference geometry, in this case, to ReferenceSquare.
        ReferenceGeometry* const referenceGeometryCodim1Ptr_;

        //! Pointer to the Codimension 2 reference geometry, in this case, to ReferenceLine.
        ReferenceGeometry* const referenceGeometryCodim2Ptr_;

        //! List of valid quadrature rules for this reference geometry
        std::vector<QuadratureRules::GaussQuadratureRule*> lstGaussQuadratureRules_;
    };
}
#endif
