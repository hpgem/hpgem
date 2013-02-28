//
//  ReferenceTetrahedron.h
//  
//
//  Created by Shavarsh Nurijanyan on 2/5/13.
//
//

#ifndef ____ReferenceTetrahedron__
#define ____ReferenceTetrahedron__

#include <iostream>

#include "ReferenceGeometry.hpp"
#include "ReferenceLine.hpp"
#include "ReferenceTriangle.hpp"
#include "GlobalNamespaceGeometry.hpp"
#include "Mappings/MappingToRefTriangleToTetrahedron.hpp"
#include "Integration/QuadratureRules/GaussQuadratureRule.hpp"

namespace Geometry
{
    /* The ordering of the vertex and faces in a tetrahedron:
     *
     * 3 o
     *   |\
     *   |  \2
     *   |  / \
     *   |/     \
     * 0 o--------o 1
     *
     */
    class ReferenceTetrahedron: public ReferenceGeometry<ThreeD>
    {
    public:
        typedef ReferenceGeometry<ThreeD> ReferenceTetrahedronT;

        using ReferenceTetrahedronT::PointReferenceT;
        using ReferenceTetrahedronT::IndexT;
        using ReferenceTetrahedronT::String;
        using ReferenceTetrahedronT::const_iterator;

    public:
        static ReferenceTetrahedron& Instance()
        {
            static ReferenceTetrahedron theInstance;
            return theInstance;
        }

        ReferenceTetrahedron();

        ReferenceTetrahedron(const ReferenceTetrahedron& copy);
        
        //! (see ReferenceGeometry.hpp)
        bool            isInternalPoint(const PointReferenceT& point) const;
        
        //! (see ReferenceGeometry.hpp)
        void            getCenter(PointReferenceT& point) const;
        
        //! (see ReferenceGeometry.hpp)
        void            getNode(const IndexT& i, PointReferenceT& point) const;
        
        //! (see ReferenceGeometry.hpp)
        String          getName() const {return "ReferenceTetrahedron";}
        
        //! Given a face index, and an index of the node position relative to the face,
        //! return the local index of the node.
        int             getLocalNodeIndex(int face, int node) const {return localNodeIndexes_[face][node];}

        /// Output routine.
        friend std::ostream& operator<<(std::ostream& os, const ReferenceTetrahedron& point);
        
        // ================================== Codimension 0 ========================================

        //! (see MappingCodimensions.hpp)
        int                                      getCodim0MappingIndex(const ListOfIndexesT&, const ListOfIndexesT&) const;

        //! (see MappingCodimensions.hpp)
        const MappingReferenceToReference<3, 3>* getCodim0MappingPtr(const IndexT) const;

        // ================================== Codimension 1 ========================================

        //! (see MappingCodimensions.hpp)
        unsigned int                             getNrOfCodim1Entities() const {return 4;}

        //! (see MappingCodimensions.hpp)
        void                                     getCodim1EntityLocalIndices(const IndexT, ListOfIndexesT& faceNodesLocal) const;

        //! (see MappingCodimensions.hpp)
        const MappingReferenceToReference<2, 3>* getCodim1MappingPtr(const IndexT) const;

        //! (see MappingCodimensions.hpp)
        const ReferenceGeometry<2>*              getCodim1ReferenceGeometry(const IndexT) const;

        // ================================== Codimension 2 ========================================

        //! (see MappingCodimensions.hpp)
        unsigned int                             getNrOfCodim2Entities() const {return 6;};

        //! (see MappingCodimensions.hpp)
        void                                     getCodim2EntityLocalIndices(const IndexT, ListOfIndexesT& faceNodesLocal) const;

        //! (see MappingCodimensions.hpp)
        const MappingReferenceToReference<1, 3>* getCodim2MappingPtr(const IndexT) const;

        //! (see MappingCodimensions.hpp)
        const ReferenceGeometry<1>*              getCodim2ReferenceGeometry(const IndexT) const;

        // ================================== Codimension 3 ========================================

        //! (see MappingCodimensions.hpp)
        unsigned int                             getNrOfCodim3Entities() const {return 4;};

        //! (see MappingCodimensions.hpp)
        void                                     getCodim3EntityLocalIndices(const IndexT, ListOfIndexesT& faceNodesLocal) const;

        // ================================== Quadrature rules =====================================

        /// Add a quadrature rule into the list of valid quadrature rules for this geometry.
        virtual void addGaussQuadratureRule(QuadratureRules::GaussQuadratureRule<3>* const qr);

        /// Get a valid quadrature for this geometry.
        virtual QuadratureRules::GaussQuadratureRule<3>* const getGaussQuadratureRule(int order) const;
            
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
        static int                                          localNodeIndexes_[4][3];
        static int                                          localNodesOnEdge_[6][2];
        
        //! Codimension 1 mappings, from a square to a tetrahedron face. TODO: Where is this used? clarify here.
        const MappingReferenceToReference<2, 3>*            mappingsTriangleToTetrahedron_[4];
        const MappingReferenceToReference<3, 3>*            mappingsTetrahedronToTetrahedron_[1];

        //! Pointer to the Codimension 1 reference geometry.
        ReferenceGeometry<2>* const                         referenceGeometryCodim1Ptr_;
        ReferenceGeometry<1>* const                         referenceGeometryCodim2Ptr_;
        
        //! List of valid quadrature rules for this reference geometry
        std::list<QuadratureRules::GaussQuadratureRule<3>*> lstGaussQuadratureRules_;
    };
};

#endif /* defined(____ReferenceTetrahedron__) */
