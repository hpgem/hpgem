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
#include "Logger.h"

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
        using ReferenceGeometryT::VectorOfReferencePointsT;
        using ListOfIndexesT = std::vector<std::size_t>;
        
        static ReferencePoint& Instance()
        {
            static ReferencePoint theInstance;
            return theInstance;
        }

        ReferencePoint(const ReferencePoint&) = delete;

        /// \brief Return true.
        bool isInternalPoint(const PointReference& p) const override final;

        /// \brief (see ReferenceGeometry.h)
        PointReference getCenter() const override final;
        
        // ================================== Codimension 0 ========================================
        
        /// \brief Returns 0.
        std::size_t getCodim0MappingIndex(const ListOfIndexesT&, const ListOfIndexesT&) const override final;

        /// \brief Returns 0.
        const MappingReferenceToReference* getCodim0MappingPtr(const std::size_t a) const override final;

        using MappingCodimensions::getCodim0MappingPtr;

        // =============================== Refinement mappings =====================================
        
        /// \brief Transform a reference point using refinement mapping
        void refinementTransform(int refineType, std::size_t subElementIdx, const PointReference& p, PointReference& pMap) const override final
        {
        }
        
        /// \brief Transformation matrix of this refinement when located on the LEFT side
        void getRefinementMappingMatrixL(int refineType, std::size_t subElementIdx, LinearAlgebra::Matrix& Q) const override final
        {
        }
        
        /// \brief Transformation matrix of this refinement when located on the RIGHT side
        void getRefinementMappingMatrixR(int refineType, std::size_t subElementIdx, LinearAlgebra::Matrix& Q) const override final
        {
        }
        
        /// \brief Refinement mapping on codim1 for a given refinement on codim0
        /// Note: this should also applied on other dimensions
        void getCodim1RefinementMappingMatrixL(int refineType, std::size_t subElementIdx, std::size_t faLocalIndex, LinearAlgebra::Matrix& Q) const override final
        {
        }
        
        /// \brief Refinement mapping on codim1 for a given refinement on codim0
        /// \brief Note: this should also applied on other dimensions
        void getCodim1RefinementMappingMatrixR(int refineType, std::size_t subElementIdx, std::size_t faLocalIndex, LinearAlgebra::Matrix& Q) const override final
        {
        }
        
    private:
        /// \brief List of valid quadrature rules for this reference geometry
        std::vector<QuadratureRules::GaussQuadratureRule*> lstGaussQuadratureRules_;

        //! Codimension 1 mappings, from a line to a line. TODO: Where is this used? clarify here.
        const MappingReferenceToReference* mappingsPointToPoint_; //!< codim0
        
        ReferencePoint();
    };
}
#endif
