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

#ifndef ____ReferencePoint__
#define ____ReferencePoint__

#include <iostream>

#include "ReferenceGeometry.h"
#include <vector>

namespace Geometry
{
    /// \class ReferencePoint
    /// \brief Reference geometry of dimensions 0.
    /// \details
    /// Most of the codimension functions are not well defined for dimension 0, so be careful.
    class ReferencePoint : public ReferenceGeometry
    {
    public:
        using ReferenceGeometryT = ReferenceGeometry;

        using ReferenceGeometryT::PointReferenceT;
        using ReferenceGeometryT::VectorOfReferencePointsT;
        using ReferenceGeometryT::IndexT;
        using ListOfIndexesT = std::vector<IndexT>;

    public:
        static ReferencePoint& Instance()
        {
            static ReferencePoint theInstance;
            return theInstance;
        }
        
    private:
        
        ReferencePoint();

        ReferencePoint(const ReferencePoint&);

    public:
        
        /// \brief Return true. TODO: Why?
        bool isInternalPoint(const PointReferenceT& p) const;

        /// \brief (see ReferenceGeometry.h)
        PointReference getCenter() const;

        /// \brief (see ReferenceGeometry.h)
        const PointReference& getNode(const IndexT& i) const;

        /// \brief (see ReferenceGeometry.h)
        String getName() const
        {
            return "ReferencePoint";
        }
        
        /// \brief Given a face index, and an index of the node position relative to the face,
        /// return the local index of the node. Dummy function, doesn't make sense for point
        std::size_t getLocalNodeIndex(std::size_t face, std::size_t node) const
        {
            std::cout << "WARNING: ReferencePoint::getLocalNodeIndex might have unexpected behaviour." << std::endl;
            return -1;
        }
        
        // ================================== Codimension 0 ========================================
        
        /// \brief Returns 0.
        std::size_t getCodim0MappingIndex(const ListOfIndexesT&, const ListOfIndexesT&) const;

        /// \brief Returns 0.
        const MappingReferenceToReference* getCodim0MappingPtr(const IndexT a) const;

        using MappingCodimensions::getCodim0MappingPtr;

        // =============================== Refinement mappings =====================================
        
        /// \brief Transform a reference point using refinement mapping
        void refinementTransform(int refineType, std::size_t subElementIdx, const PointReferenceT& p, PointReferenceT& pMap) const
        {
        }
        
        /// \brief Transformation matrix of this refinement when located on the LEFT side
        void getRefinementMappingMatrixL(int refineType, std::size_t subElementIdx, LinearAlgebra::Matrix& Q) const
        {
        }
        
        /// \brief Transformation matrix of this refinement when located on the RIGHT side
        void getRefinementMappingMatrixR(int refineType, std::size_t subElementIdx, LinearAlgebra::Matrix& Q) const
        {
        }
        
        /// \brief Refinement mapping on codim1 for a given refinement on codim0
        /// Note: this should also applied on other dimensions
        void getCodim1RefinementMappingMatrixL(int refineType, std::size_t subElementIdx, std::size_t faLocalIndex, LinearAlgebra::Matrix& Q) const
        {
        }
        
        /// \brief Refinement mapping on codim1 for a given refinement on codim0
        /// \brief Note: this should also applied on other dimensions
        void getCodim1RefinementMappingMatrixR(int refineType, std::size_t subElementIdx, std::size_t faLocalIndex, LinearAlgebra::Matrix& Q) const
        {
        }
        
        /// \brief List of valid quadrature rules for this reference geometry
        std::vector<QuadratureRules::GaussQuadratureRule*> lstGaussQuadratureRules_;

        //! Codimension 1 mappings, from a line to a line. TODO: Where is this used? clarify here.
        const MappingReferenceToReference* mappingsPointToPoint_; //!< codim0
    };
}
;
#endif
