//
//  ReferenceLine.h
//  
//
//  Created by Shavarsh Nurijanyan on 2/4/13.
//
//

#ifndef ____ReferenceLine__
#define ____ReferenceLine__

#include <iostream>
using std::ostream;

#include "ReferencePoint.hpp"
#include "ReferenceGeometry.hpp"
#include "GlobalNamespaceGeometry.hpp"
#include "Mappings/MappingReferenceToReference.hpp"
#include "Integration/QuadratureRules/GaussQuadratureRule.hpp"

namespace Geometry
{
    /* Behold the reference line:
     *
     * (-1) 0---1---1 (+1)
     *
     */
    class ReferenceLine: public ReferenceGeometry<OneD>
    {

    public:
        typedef ReferenceGeometry<OneD> ReferenceGeometryT;
        
        using ReferenceGeometryT::PointReferenceT;
        using ReferenceGeometryT::VectorOfReferencePointsT;
        using ReferenceGeometryT::IndexT;
        using ReferenceGeometryT::String;
        using ReferenceGeometryT::const_iterator;
        typedef std::vector<IndexT> ListOfIndexesT;
        typedef MappingReferenceToReference<1,1> Ref1ToRef1MappingT; // Numbers indicate dim.
        typedef MappingReferenceToReference<0,1> Ref0ToRef1MappingT;

    public:
        static ReferenceLine& Instance()
        {
            static ReferenceLine theInstance;
            return theInstance;
        }

        ReferenceLine();

        ReferenceLine(const ReferenceLine&);
        
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
        int             getLocalNodeIndex(int face, int node) const {return localNodeIndexes_[face][node];}

        //! Output routine.
        friend ostream& operator<<(ostream& os, const ReferenceLine& point);

        // ================================== Codimension 0 ========================================

        //! (see MappingCodimensions.hpp)
        int                                         getCodim0MappingIndex(const ListOfIndexesT&, const ListOfIndexesT&) const;

        //! (see MappingCodimensions.hpp)
        const MappingReferenceToReference<1, 1>*    getCodim0MappingPtr(const IndexT) const;

        // ================================== Codimension 1 ========================================

        //! (see MappingCodimensions.hpp)
        unsigned int                                getNrOfCodim1Entities() const {return 2;}

        //! (see MappingCodimensions.hpp)
        void                                        getCodim1EntityLocalIndices(const IndexT, ListOfIndexesT& faceNodesLocal) const;

        //! (see MappingCodimensions.hpp)
        const MappingReferenceToReference<0, 1>*    getCodim1MappingPtr(const IndexT) const;

        //! (see MappingCodimensions.hpp)
        const ReferenceGeometry<0>*                 getCodim1ReferenceGeometry(const IndexT) const;
        
        // ================================== Quadrature rules =====================================

        /// Add a quadrature rule into the list of valid quadrature rules for this geometry.
        virtual void addGaussQuadratureRule(QuadratureRules::GaussQuadratureRule<1>* const qr);

        /// Get a valid quadrature for this geometry.
        virtual QuadratureRules::GaussQuadratureRule<1>* const getGaussQuadratureRule(int order) const;

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

        //! Codimension 0 mappings, from a point to a line. TODO: Where is this used? clarify here.
        const Ref0ToRef1MappingT*                           mappingsPointToLine_[2];

        //! Pointer to the Codimension 1 reference geometry, in this case, to ReferencePoint.
        ReferenceGeometry<0>* const                         referenceGeometryCodim1Ptr_;
        
        //! List of valid quadrature rules for this reference geometry
        std::list<QuadratureRules::GaussQuadratureRule<1>*> lstGaussQuadratureRules_;

    };
};
#endif
