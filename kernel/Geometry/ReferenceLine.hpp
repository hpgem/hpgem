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

#ifndef ____ReferenceLine__
#define ____ReferenceLine__

#include "ReferenceGeometry.hpp"

#include <iostream>
#include <vector>

namespace Geometry
{
    /* Behold the reference line:
     *
     * (-1) 0---1---1 (+1)
     *
     */
    class ReferenceLine: public ReferenceGeometry
    {

    public:
        typedef ReferenceGeometry ReferenceGeometryT;
        
        using ReferenceGeometryT::PointReferenceT;
        using ReferenceGeometryT::VectorOfReferencePointsT;
        using ReferenceGeometryT::IndexT;
        using ReferenceGeometryT::String;
        using ReferenceGeometryT::const_iterator;
        typedef std::vector<unsigned int>                 ListOfIndexesT;
        typedef MappingReferenceToReference    Ref1ToRef1MappingT; // Numbers indicate dim.
        typedef MappingReferenceToReference    Ref0ToRef1MappingT;

    public:
        static ReferenceLine& Instance()
        {
            static ReferenceLine theInstance;
            return theInstance;
        }

    private:

        ReferenceLine();

        ReferenceLine(const ReferenceLine&);
        
    public:

        //! (see ReferenceGeometry.hpp)
        bool            isInternalPoint(const PointReferenceT&) const;

        //! (see ReferenceGeometry.hpp)
        void            getCenter(PointReferenceT&) const;

        //! (see ReferenceGeometry.hpp)
        void            getNode(const IndexT& i, PointReferenceT& point) const;

        //! (see ReferenceGeometry.hpp)
        String          getName() const {return "ReferenceLine";}

        //! Given a face index, and an index of the node position relative to the face,
        //! return the local index of the node.
        int             getLocalNodeIndex(int face, int node) const {return localNodeIndexes_[face][node];}///\bug this name makes me expect a map from a global node index to a local one

        //! Output routine.
        friend std::ostream& operator<<(std::ostream& os, const ReferenceLine& point);

        // ================================== Codimension 0 ========================================

        //! (see MappingCodimensions.hpp)
        int                                         getCodim0MappingIndex(const ListOfIndexesT&, const ListOfIndexesT&) const;

        //! (see MappingCodimensions.hpp)
        const MappingReferenceToReference*    getCodim0MappingPtr(const IndexT) const;

        using MappingCodimensions::getCodim0MappingPtr;

        // ================================== Codimension 1 ========================================

        //! (see MappingCodimensions.hpp)
        unsigned int                                getNrOfCodim1Entities() const {return 2;}

        //! (see MappingCodimensions.hpp)
        void                                        getCodim1EntityLocalIndices(const IndexT, ListOfIndexesT& faceNodesLocal) const;

        //! (see MappingCodimensions.hpp)
        const MappingReferenceToReference*    getCodim1MappingPtr(const IndexT) const;

        //! (see MappingCodimensions.hpp)
        const ReferenceGeometry*                 getCodim1ReferenceGeometry(const IndexT) const;

        // =============================== Refinement mappings =====================================
        
        //! Transform a reference point using refinement mapping
        void refinementTransform(int refineType, int subElementIdx, 
                      const PointReferenceT& p, PointReferenceT& pMap) const {}

        //! Transformation matrix of this refinement when located on the LEFT side
        void getRefinementMappingMatrixL(int refineType, int subElementIdx, 
                    LinearAlgebra::Matrix& Q) const {}

        //! Transformation matrix of this refinement when located on the RIGHT side
        void getRefinementMappingMatrixR(int refineType, int subElementIdx, 
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
        //! See top comment for the corresponding numbering. (In a line, the 'faces' are nodes,
        //! and nodes are just a coordinate).
        static int                                          localNodeIndexes_[2][1];

        //! Codimension 1 mappings, from a line to a line. TODO: Where is this used? clarify here.
        const Ref1ToRef1MappingT*                           mappingsLineToLine_[2]; //!< codim0

        //! Codimension 0 mappings, from a point to a line. TODO: Where is this used? clarify here. This is the 1D face->element map
        const Ref0ToRef1MappingT*                           mappingsPointToLine_[2];

        //! Pointer to the Codimension 1 reference geometry, in this case, to ReferencePoint.
        ReferenceGeometry* const                         referenceGeometryCodim1Ptr_;
        
        //! List of valid quadrature rules for this reference geometry
        std::vector<QuadratureRules::GaussQuadratureRule*> lstGaussQuadratureRules_;

    };
};
#endif
