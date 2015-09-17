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
#ifndef ____ReferencePyramid__
#define ____ReferencePyramid__

#include "ReferenceGeometry.h"
#include <vector>
// created for the shape globally
namespace Geometry
{
    /* The ordering of the vertex and faces in a pyramid (top view; 0 is above the other nodes):
     *
     * 3 o--------o 4
     *   |\     / |
     *   |  \ /   |
     *   | 0 o    |
     *   |  / \   |
     *   |/     \ |
     * 1 o--------o 2
     *
     * Note that the axes are oriented differently than for the other reference geometries
     */
    class ReferencePyramid : public ReferenceGeometry
    {
    public:
        using ReferenceGeometry::String;

        static ReferencePyramid& Instance()
        {
            static ReferencePyramid theInstance;
            return theInstance;
        }
        
        ReferencePyramid(const ReferencePyramid& copy) = delete;

    public:
        
        //! (see ReferenceGeometry.h)
        bool isInternalPoint(const PointReference<3>& point) const override final;
        
        /// Output routine.
        friend std::ostream& operator<<(std::ostream& os, const ReferencePyramid& point);

        const PointReferenceBase& getCenter() const override final
        {
            return *center_;
        }

        std::size_t getNumberOfNodes() const override final
        {
            return 5;
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
            return 5;
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
            return 8;
        }
        
        //! (see MappingCodimensions.h)
        std::vector<std::size_t> getCodim2EntityLocalIndices(const std::size_t) const override final;

        //! (see MappingCodimensions.h)
        const MappingReferenceToReference<2>* getCodim2MappingPtr(const std::size_t) const override final;

        //! (see MappingCodimensions.h)
        const ReferenceGeometry* getCodim2ReferenceGeometry(const std::size_t) const override final;

        // ================================== Codimension 3 ========================================
        
        //! (see MappingCodimensions.h)
        std::size_t getNrOfCodim3Entities() const override final
        {
            return 5;
        }
        
        //! (see MappingCodimensions.h)
        std::vector<std::size_t> getCodim3EntityLocalIndices(const std::size_t) const override final;

        // =============================== Refinement mappings =====================================
        
        //! Transform a reference point using refinement mapping
        void refinementTransform(int refineType, std::size_t subElementIdx, const PointReference<3>& p, PointReference<3>& pMap) const override final
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
        
        ReferencePyramid();
        
        //! Local node indexes contains the numbering of the vertex of the shape, ordered by faces.
        //! See top comment for the corresponding numbering.
        static std::size_t localNodeIndexes_[5][4];

        //! The nodes on edge contains the local index of the two nodes in every edge.
        static std::size_t localNodesOnEdge_[8][2];

        //! Pointer to the Codimension 1 reference geometry: square.
        ReferenceGeometry* const referenceGeometryCodim1SquarePtr_;

        //! Pointer to the Codimension 1 reference geometry: triangle.
        ReferenceGeometry* const referenceGeometryCodim1TrianglePtr_;

        //! Pointer to the Codimension 2 reference geometry, in this case, to ReferenceLine.
        ReferenceGeometry* const referenceGeometryCodim2Ptr_;

        //! Codimension 1 mappings, from a face to a pyramid. (used to map a coordinate to a pyramid from one of its faces)
        const MappingReferenceToReference<1>* mappingsFaceToPyramid_[5];

        //! Codimension 0 mappings, from a Pyramid to a pyramid. Only needed when an element has a pyramid as a face.
        //const MappingReferenceToReference<0>* mappingsPyramidToPyramid_[1];

        //! List of valid quadrature rules for this reference geometry
        std::vector<QuadratureRules::GaussQuadratureRule*> lstGaussQuadratureRules_;

        std::vector<const PointReference<3>* > points_;

        const PointReference<3>* center_;
    };
}
#endif
