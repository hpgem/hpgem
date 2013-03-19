//
//  ReferenceTriangle.h
//  
//
//  Created by Shavarsh Nurijanyan on 2/5/13.
//
//

#ifndef ____ReferenceTriangle__
#define ____ReferenceTriangle__

#include <iostream>
using std::ostream;

#include "GlobalNamespaceGeometry.hpp"
#include "ReferenceLine.hpp"
#include "Mappings/MappingToRefLineToTriangle.hpp"
#include "Mappings/MappingToRefTriangleToTriangle.hpp"
#include "Integration/QuadratureRules/GaussQuadratureRule.hpp"
#include "Integration/QuadratureRules/GaussQuadratureRulesForTriangle.hpp"//do not take this out, quadrature rules are


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
    class ReferenceTriangle: public ReferenceGeometry<TwoD>
    {

    public:
        typedef ReferenceGeometry<TwoD> ReferenceGeometryT;
        
        using ReferenceGeometryT::IndexT;
        using ReferenceGeometryT::PointReferenceT;
        using ReferenceGeometryT::String;
        using ReferenceGeometryT::ListOfIndexesT;

    public:
        static ReferenceTriangle& Instance()
        {
            static ReferenceTriangle theInstance;
            return theInstance;
        }

        ReferenceTriangle();

        ReferenceTriangle(const ReferenceTriangle& copy);
        
        /// /see (see ReferenceGeometry.hpp)
        bool            isInternalPoint(const PointReferenceT& point) const;
        
        //! (see ReferenceGeometry.hpp)
        void            getCenter(PointReferenceT& point) const;
        
        //! (see ReferenceGeometry.hpp)
        void            getNode(const IndexT& i, PointReferenceT& point) const;
        
        //! (see ReferenceGeometry.hpp)
        String          getName() const {return "ReferenceTriangle";}
        
        //! Given a face index, and an index of the node position relative to the face,
        //! return the local index of the node.
        int             getLocalNodeIndex(int face, int node) const {return localNodeIndexes_[face][node];}

        /// Output routine.
        friend std::ostream& operator<<(std::ostream& os, const ReferenceTriangle& point);

        // ================================== Codimension 0 ========================================

        //! (see MappingCodimensions.hpp)
        int                                      getCodim0MappingIndex(const ListOfIndexesT&, const ListOfIndexesT&) const;

        //! (see MappingCodimensions.hpp)
        const MappingReferenceToReference<2, 2>* getCodim0MappingPtr(const IndexT) const;

        // ================================== Codimension 1 ========================================

        //! (see MappingCodimensions.hpp)
        unsigned int getNrOfCodim1Entities() const {return 3;}

        //! (see MappingCodimensions.hpp)
        void                                     getCodim1EntityLocalIndices(const IndexT, ListOfIndexesT& faceNodesLocal) const;

        //! (see MappingCodimensions.hpp)
        const MappingReferenceToReference<1, 2>* getCodim1MappingPtr(const IndexT) const;

        //! (see MappingCodimensions.hpp)
        const ReferenceGeometry<1>*              getCodim1ReferenceGeometry(const IndexT) const;

        // ================================== Codimension 2 ========================================

        //! (see MappingCodimensions.hpp)
        int                                      getNrOfCodim2Entities() const {return 3;};

        void                                     getCodim2EntityLocalIndices(const int, std::vector<unsigned int>&) const {return;}
  
        
        // ================================== Quadrature rules =====================================

        /// Add a quadrature rule into the list of valid quadrature rules for this geometry.
        virtual void addGaussQuadratureRule(QuadratureRules::GaussQuadratureRule<2>* const qr);

        /// Get a valid quadrature for this geometry.
        virtual QuadratureRules::GaussQuadratureRule<2>* const getGaussQuadratureRule(int order) const;
            
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
        //! See top comment for the corresponding numbering.
        static int                                          localNodeIndexes_[3][2];

        //! Codimension 1 mappings, from a line to a triangle. TODO: Where is this used? clarify here.
        const MappingReferenceToReference<1, 2>*            mappingsLineToTriangle_[3];

        //! Codimension 0 mappings, from a triangle to a triangle. TODO: Where is this used? clarifiy here.
        const MappingReferenceToReference<2, 2>*            mappingsTriangleToTriangle_[6];

        //! Pointer to the Codimension 1 reference geometry, in this case, to ReferenceLine.
        ReferenceGeometry<1>* const                         referenceGeometryCodim1Ptr_;
        
        //! List of valid quadrature rules for this reference geometry
        std::list<QuadratureRules::GaussQuadratureRule<2>*> lstGaussQuadratureRules_;
    };
};
#endif /* defined(____ReferenceTriangle__) */
