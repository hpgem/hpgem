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

#ifndef ____ReferenceTriangle__
#define ____ReferenceTriangle__

#include "ReferenceGeometry.h"

#include <iostream>
#include <vector>

namespace Geometry
{
    
    /* The ordering of the vertex and faces in a triangle:
     *
     *   (0,1) 2
     *         | \
     *         1   2
     *         |     \
     *   (0,0) 0---0---1 (1,0)
     *
     */
    class ReferenceTriangle : public ReferenceGeometry
    {
    public:
        using ReferenceGeometry::String;
        using ReferenceGeometry::ListOfIndexesT;
        
        static ReferenceTriangle& Instance()
        {
            static ReferenceTriangle theInstance;
            return theInstance;
        }

        ReferenceTriangle(const ReferenceTriangle& copy) = delete;

        /// /see (see ReferenceGeometry.h)
        bool isInternalPoint(const PointReference<2>& point) const override final;
        
        /// Output routine.
        friend std::ostream& operator<<(std::ostream& os, const ReferenceTriangle& point);

        const PointReferenceBase& getCenter() const override final
        {
            return *center_;
        }

        std::size_t getNumberOfNodes() const override final
        {
            return 3;
        }

        const PointReferenceBase& getReferenceNodeCoordinate(const std::size_t& i) const override final
        {
            logger.assert(i < getNumberOfNodes(), "Asked for node %, but there are only % nodes", i, getNumberOfNodes());
            return *points_[i];
        }

        // ================================== Codimension 0 ========================================
        
        //! (see MappingCodimensions.h)
        std::size_t getCodim0MappingIndex(const ListOfIndexesT&, const ListOfIndexesT&) const override final;

        //! (see MappingCodimensions.h)
        const MappingReferenceToReference<0>* getCodim0MappingPtr(const std::size_t) const override final;

        using MappingCodimensions::getCodim0MappingPtr;

        // ================================== Codimension 1 ========================================
        
        //! (see MappingCodimensions.h)
        std::size_t getNrOfCodim1Entities() const override final
        {
            return 3;
        }
        
        //! (see MappingCodimensions.h)
        std::vector<std::size_t> getCodim1EntityLocalIndices(const std::size_t) const override final;

        //! (see MappingCodimensions.h)
        const MappingReferenceToReference<1>* getCodim1MappingPtr(const std::size_t) const override final;

        //! (see MappingCodimensions.h)
        const ReferenceGeometry* getCodim1ReferenceGeometry(const std::size_t) const override final;

        // ================================== Codimension 2 ========================================
        
        //! (see MappingCodimensions.h)
        std::size_t getNrOfCodim2Entities() const override final
        {
            return 3;
        }

        std::vector<std::size_t> getCodim2EntityLocalIndices(const std::size_t node) const override final
        {
            return std::vector<std::size_t>(1, node);
        }
        
        const ReferenceGeometry* getCodim2ReferenceGeometry(const std::size_t) const override final;

        // =============================== Refinement mappings =====================================
        
        //! Transform a reference point using refinement mapping
        void refinementTransform(int refineType, std::size_t subElementIdx, const PointReference<2>& p, PointReference<2>& pMap) const override final
        {
        }
        
        //! Transformation matrix of this refinement when located on the LEFT side
        void getRefinementMappingMatrixL(int refineType, std::size_t subElementIdx, LinearAlgebra::MiddleSizeMatrix& Q) const override final
        {
        }
        
        //! Transformation matrix of this refinement when located on the RIGHT side
        void getRefinementMappingMatrixR(int refineType, std::size_t subElementIdx, LinearAlgebra::MiddleSizeMatrix& Q) const override final
        {
        }
        
        //! Refinement mapping on codim1 for a given refinement on codim0
        //! Note: this should also applied on other dimensions
        void getCodim1RefinementMappingMatrixL(int refineType, std::size_t subElementIdx, std::size_t faLocalIndex, LinearAlgebra::MiddleSizeMatrix& Q) const override final
        {
        }
        
        //! Refinement mapping on codim1 for a given refinement on codim0
        //! Note: this should also applied on other dimensions
        void getCodim1RefinementMappingMatrixR(int refineType, std::size_t subElementIdx, std::size_t faLocalIndex, LinearAlgebra::MiddleSizeMatrix& Q) const override final
        {
        }
        
    private:
        
        ReferenceTriangle();
        
        //! Local node indexes contains the numbering of the vertex of the shape, ordered by faces.
        //! See top comment for the corresponding numbering.
        static std::size_t localNodeIndexes_[3][2];

        //! Codimension 1 mappings, from a line to a triangle. TODO: Where is this used? clarify here.
        const MappingReferenceToReference<1>* mappingsLineToTriangle_[3];

        //! Codimension 0 mappings, from a triangle to a triangle. TODO: Where is this used? clarifiy here.
        const MappingReferenceToReference<0>* mappingsTriangleToTriangle_[6];

        //! Pointer to the Codimension 1 reference geometry, in this case, to ReferenceLine.
        ReferenceGeometry * const referenceGeometryCodim1Ptr_;

        //! List of valid quadrature rules for this reference geometry
        std::vector<QuadratureRules::GaussQuadratureRule*> lstGaussQuadratureRules_;

        std::vector<const PointReference<2>* > points_;

        const PointReference<2>* center_;
    };

}
#endif /* defined(____ReferenceTriangle__) */
