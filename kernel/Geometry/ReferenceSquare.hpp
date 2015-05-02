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
#ifndef ____ReferenceSquare__
#define ____ReferenceSquare__

#include "ReferenceGeometry.hpp"

#include <vector>
#include <iostream>

namespace Geometry
{

    /* Behold the reference square:
     *
     * (-1,+1) 2---3---3 (+1,+1)
     *         |       |
     *         1       2
     *         |       |
     * (-1,-1) 0---0---1 (+1,-1)
     *
     */
    class ReferenceSquare : public ReferenceGeometry
    {
    public:
        using ReferenceGeometryT = ReferenceGeometry ;
        //         typedef QuadratureRules::GaussQuadratureRule<TwoD>  GaussQuadratureRuleT;
        //         typedef std::vector<QuadratureRules::GaussQuadratureRule<TwoD>*>            ListOfGaussQuadratureRulePtrT;
        using ReferenceGeometryT::IndexT;
        using ReferenceGeometryT::PointReferenceT;
        using ReferenceGeometryT::String;
        using ReferenceGeometryT::ListOfIndexesT;

    public:

        static ReferenceSquare& Instance()
        {
            static ReferenceSquare theInstance;
            return theInstance;
        }

    private:

        ReferenceSquare();

        ReferenceSquare(const ReferenceSquare& copy);

    public:

        //! (see ReferenceGeometry.hpp)
        bool isInternalPoint(const PointReferenceT& point) const;

        //! (see ReferenceGeometry.hpp)
        void getCenter(PointReferenceT& point) const;

        //! (see ReferenceGeometry.hpp)
        void getNode(const IndexT& i, PointReferenceT& point) const;

        //! (see ReferenceGeometry.hpp)

        String getName() const
        {
            return "ReferenceSquare";
        }

        //! Given a face index, and an index of the node position relative to the face,
        //! return the local index of the node.

        std::size_t getLocalNodeIndex(std::size_t face, std::size_t node) const
        {
            return localNodeIndexes_[face][node];
        }

        //! (see ReferenceGeometry.hpp) //! (see ReferenceGeometry.hpp) duplicating the referenceGeometry.getNumberofNodes(), thus commented out
        //IndexT          getId() const {return 4;}

        //! Output routine.
        friend std::ostream& operator<<(std::ostream& os, const ReferenceSquare& point);

        // ================================== Codimension 0 ========================================

        //! (see MappingCodimensions.hpp)
        std::size_t getCodim0MappingIndex(const ListOfIndexesT&, const ListOfIndexesT&) const;

        //! (see MappingCodimensions.hpp)
        const MappingReferenceToReference* getCodim0MappingPtr(const IndexT) const;

        using MappingCodimensions::getCodim0MappingPtr;

        // ================================== Codimension 1 ========================================

        //! (see MappingCodimensions.hpp)

        std::size_t getNrOfCodim1Entities() const
        {
            return 4;
        }

        //! (see MappingCodimensions.hpp)
        void getCodim1EntityLocalIndices(const IndexT, ListOfIndexesT& faceNodesLocal) const;

        //! (see MappingCodimensions.hpp)
        const MappingReferenceToReference* getCodim1MappingPtr(const IndexT) const;

        //! (see MappingCodimensions.hpp)
        const ReferenceGeometry* getCodim1ReferenceGeometry(const IndexT) const;

        // ================================== Codimension 2 ========================================

        //! (see MappingCodimensions.hpp)

        std::size_t getNrOfCodim2Entities() const
        {
            return 4;
        };

        void getCodim2EntityLocalIndices(const std::size_t vertex, std::vector<std::size_t>& vertexNodesLocal) const
        {
            vertexNodesLocal[0] = vertex;
            return;
        }

        const ReferenceGeometry* getCodim2ReferenceGeometry(const std::size_t) const;

        // =============================== Refinement mappings =====================================

        //! Transform a reference point using refinement mapping
        void refinementTransform(int refineType, std::size_t subElementIdx,
                const PointReferenceT& p, PointReferenceT& pMap) const;

        //! Transformation matrix of this refinement when located on the LEFT side
        void getRefinementMappingMatrixL(int refineType, std::size_t subElementIdx,
                LinearAlgebra::Matrix& Q) const;

        //! Transformation matrix of this refinement when located on the RIGHT side
        void getRefinementMappingMatrixR(int refineType, std::size_t subElementIdx,
                LinearAlgebra::Matrix& Q) const;

        //! Refinement mapping on codim1 for a given refinement on codim0
        //! Note: this should also applied on other dimensions
        void getCodim1RefinementMappingMatrixL(int refineType, DimT subElementIdx,
                DimT faLocalIndex, LinearAlgebra::Matrix& Q) const;

        //! Refinement mapping on codim1 for a given refinement on codim0
        //! Note: this should also applied on other dimensions
        void getCodim1RefinementMappingMatrixR(int refineType, DimT subElementIdx,
                DimT faLocalIndex, LinearAlgebra::Matrix& Q) const;

    private:
        //! Local node indexes contains the numbering of the vertex of the shape, ordered by faces.
        //! See top comment for the corresponding numbering.
        static std::size_t localNodeIndexes_[4][2];

        //! Codimension 1 mappings, from a line to a square. TODO: Where is this used? clarify here.
        const MappingReferenceToReference* mappingsLineToSquare_[4];

        //! Codimension 0 mappings, from a square to a square. TODO: Where is this used? clarifiy here.
        const MappingReferenceToReference* mappingsSquareToSquare_[8];

        //! Pointer to the Codimension 1 reference geometry, in this case, to ReferenceLine.
        ReferenceGeometry * const referenceGeometryCodim1Ptr_;

        //! List of valid quadrature rules for this reference geometry
        std::vector<QuadratureRules::GaussQuadratureRule*> lstGaussQuadratureRules_;
    };

}
;
#endif
