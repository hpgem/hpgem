//
//  ReferenceRectangle.hpp
//  
//
//  Created by Shavarsh Nurijanyan on 2/5/13.
//
//

#ifndef ____ReferenceSquare__
#define ____ReferenceSquare__

#include <list>
#include <iostream>
using std::ostream;

#include "GlobalNamespaceGeometry.hpp"
#include "ReferenceLine.hpp"
#include "Mappings/MappingToRefLineToSquare.hpp"
#include "Mappings/MappingToRefSquareToSquare.hpp"
#include "Geometry/Mappings/RefinementMapping.hpp"
#include "Integration/QuadratureRules/GaussQuadratureRule.hpp"

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
    class ReferenceSquare: public ReferenceGeometry
    {

    public:
        typedef ReferenceGeometry                     ReferenceGeometryT;
//         typedef QuadratureRules::GaussQuadratureRule<TwoD>  GaussQuadratureRuleT;
//         typedef std::list<QuadratureRules::GaussQuadratureRule<TwoD>*>            ListOfGaussQuadratureRulePtrT;
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
        
        ReferenceSquare();

        ReferenceSquare(const ReferenceSquare& copy);

        //! (see ReferenceGeometry.hpp)
        bool            isInternalPoint(const PointReferenceT& point) const;

        //! (see ReferenceGeometry.hpp)
        void            getCenter(PointReferenceT& point) const;

        //! (see ReferenceGeometry.hpp)
        void            getNode(const IndexT& i, PointReferenceT& point) const;

        //! (see ReferenceGeometry.hpp)
        String          getName() const {return "ReferenceSquare";}
        
        //! Given a face index, and an index of the node position relative to the face,
        //! return the local index of the node.
        int             getLocalNodeIndex(int face, int node) const{return localNodeIndexes_[face][node];}

            //! (see ReferenceGeometry.hpp) //! (see ReferenceGeometry.hpp) duplicating the referenceGeometry.getNumberofNodes(), thus commented out
            //IndexT          getId() const {return 4;}

        //! Output routine.
        friend std::ostream& operator<<(std::ostream& os, const ReferenceSquare& point);

        // ================================== Codimension 0 ========================================

        //! (see MappingCodimensions.hpp)
        int                                      getCodim0MappingIndex(const ListOfIndexesT&, const ListOfIndexesT&) const;

        //! (see MappingCodimensions.hpp)
        const MappingReferenceToReference* getCodim0MappingPtr(const IndexT) const;

        // ================================== Codimension 1 ========================================

        //! (see MappingCodimensions.hpp)
        unsigned int                             getNrOfCodim1Entities() const {return 4;}

        //! (see MappingCodimensions.hpp)
        void                                     getCodim1EntityLocalIndices(const IndexT, ListOfIndexesT& faceNodesLocal) const;

        //! (see MappingCodimensions.hpp)
        const MappingReferenceToReference* getCodim1MappingPtr(const IndexT) const;
        
        //! (see MappingCodimensions.hpp)
        const ReferenceGeometry*              getCodim1ReferenceGeometry(const IndexT) const;

        // ================================== Codimension 2 ========================================

        //! (see MappingCodimensions.hpp)
        unsigned int                                      getNrOfCodim2Entities() const {return 4;};

        void                                     getCodim2EntityLocalIndices(const int, std::vector<unsigned int>&) const {return;}

        // ========================= Quadrature rules for this geometry  ===========================
        
        /// Add a quadrature rule into the list of valid quadrature rules for this geometry.
        virtual void addGaussQuadratureRule(QuadratureRules::GaussQuadratureRule* const qr);

        /// Get a valid quadrature for this geometry.
        virtual QuadratureRules::GaussQuadratureRule* const getGaussQuadratureRule(int order) const;
        
        // =============================== Refinement mappings =====================================
        
        //! Transform a reference point using refinement mapping
        void refinementTransform(int refineType, int subElementIdx, 
                      const PointReferenceT& p, PointReferenceT& pMap) const;

        //! Transformation matrix of this refinement when located on the LEFT side
        void getRefinementMappingMatrixL(int refineType, int subElementIdx, 
                    LinearAlgebra::Matrix& Q) const;

        //! Transformation matrix of this refinement when located on the RIGHT side
        void getRefinementMappingMatrixR(int refineType, int subElementIdx, 
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
        static int                                  localNodeIndexes_[4][2];

        //! Codimension 1 mappings, from a line to a square. TODO: Where is this used? clarify here.
        const MappingReferenceToReference*    mappingsLineToSquare_[4];

        //! Codimension 0 mappings, from a square to a square. TODO: Where is this used? clarifiy here.
        const MappingReferenceToReference*    mappingsSquareToSquare_[8];

        //! Pointer to the Codimension 1 reference geometry, in this case, to ReferenceLine.
        ReferenceGeometry* const                 referenceGeometryCodim1Ptr_;
        
        //! List of valid quadrature rules for this reference geometry
        std::list<QuadratureRules::GaussQuadratureRule*> lstGaussQuadratureRules_;
    };
};
#endif
