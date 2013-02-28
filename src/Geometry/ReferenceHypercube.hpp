//
//  ReferenceHypercube.h
//  
//
//  Created by Shavarsh Nurijanyan on 2/7/13.
//
//

#ifndef REFERENCEHYPERCUBE_HH
#define REFERENCEHYPERCUBE_HH

#include "ReferenceGeometry.hpp"
#include "ReferenceLine.hpp"
#include "ReferenceSquare.hpp"
#include "ReferenceCube.hpp"
#include "Mappings/MappingReferenceToReference.hpp"
#include "Mappings/MappingToRefCubeToHypercube.hpp"
#include "GlobalNamespaceGeometry.hpp"
#include "Integration/QuadratureRules/GaussQuadratureRule.hpp"

namespace Geometry
{
    // Reference Hypercube.
    class ReferenceHypercube: public ReferenceGeometry<FourD>
    {
    public:
        typedef ReferenceGeometry<FourD> ReferenceGeometryT;
        
        using typename ReferenceGeometryT::IndexT;
        using typename ReferenceGeometryT::PointReferenceT;
        using typename ReferenceGeometryT::VectorOfReferencePointsT;
        using ReferenceGeometryT::String;
        
    public:
        static ReferenceHypercube& Instance()
        {
            static ReferenceHypercube theInstance;
            return theInstance;
        }
        
        ReferenceHypercube();

        ReferenceHypercube(const ReferenceHypercube& copy);
        
        //! (see ReferenceGeometry.hpp)
        bool            isInternalPoint(const PointReferenceT& point) const;
        
        //! (see ReferenceGeometry.hpp)
        void            getCenter(PointReferenceT& point) const;
        
        //! (see ReferenceGeometry.hpp)
        void            getNode(const IndexT& i, PointReferenceT& point) const;
        
        //! (see ReferenceGeometry.hpp)
        String          getName() const {return "ReferenceHypercube";}
        
        //! Given a face index, and an index of the node position relative to the face,
        //! return the local index of the node.
        int             getLocalNodeIndex(int face, int node) const{return localNodeIndexes_[face][node];}
        
        /// Output routine.
        friend std::ostream& operator<<(std::ostream& os, const ReferenceHypercube& point);
        
        // ================================== Codimension 0 ========================================

        //! (see MappingCodimensions.hpp)
        int                                      getCodim0MappingIndex(const ListOfIndexesT&, const ListOfIndexesT&) const;

        //! (see MappingCodimensions.hpp)
        const MappingReferenceToReference<4, 4>* getCodim0MappingPtr(const IndexT) const;

        // ================================== Codimension 1 ========================================

        //! (see MappingCodimensions.hpp)
        unsigned int                             getNrOfCodim1Entities() const {return 8;} // 'faces' (cubes)

        //! (see MappingCodimensions.hpp)
        void                                     getCodim1EntityLocalIndices(const IndexT, ListOfIndexesT& faceNodesLocal) const;

        //! (see MappingCodimensions.hpp)
        const MappingReferenceToReference<3, 4>* getCodim1MappingPtr(const IndexT) const;

        //! (see MappingCodimensions.hpp)
        const ReferenceGeometry<3>*              getCodim1ReferenceGeometry(const IndexT) const;

        // ================================== Codimension 2 ========================================

        //! (see MappingCodimensions.hpp)
        unsigned int                             getNrOfCodim2Entities() const {return 24;} // 'edges' (faces)

        //! (see MappingCodimensions.hpp)
        void                                     getCodim2EntityLocalIndices(const IndexT, ListOfIndexesT& faceNodesLocal) const;

        //! (see MappingCodimensions.hpp)
        const MappingReferenceToReference<2, 4>* getCodim2MappingPtr(const IndexT) const;

        //! (see MappingCodimensions.hpp)
        const ReferenceGeometry<2>*              getCodim2ReferenceGeometry(const IndexT) const;

        // ================================== Codimension 3 ========================================

        //! (see MappingCodimensions.hpp)
        unsigned int                             getNrOfCodim3Entities() const {return 16;} // 'vertices' (edges)

        //! (see MappingCodimensions.hpp)
        void                                     getCodim3EntityLocalIndices(const IndexT, ListOfIndexesT&) const;

        // ================================== Quadrature rules =====================================

        /// Add a quadrature rule into the list of valid quadrature rules for this geometry.
        virtual void addGaussQuadratureRule(QuadratureRules::GaussQuadratureRule<4>* const qr);

        /// Get a valid quadrature for this geometry.
        virtual QuadratureRules::GaussQuadratureRule<4>* const getGaussQuadratureRule(int order) const;

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
        static int                                          localNodeIndexes_[8][8]; // 8 'faces' (cubes) with 8 vertex.

        //! The nodes on edge contains the local index of the two nodes in every edge.
        static int                                          localNodesOnEdge_[24][2]; //!< 24 edges with 2 nodes

        //! Codimension 1 mappings, from a cube to a hypercube. TODO: Where is this used? clarify here.
        const MappingReferenceToReference<3, 4>*            mappingsCubeToHypercube_[8];

        //! Only requiered for 5D elements
        //const MappingReferenceToReference<4, 4>* mappingsHypercubeToHypercube_[8];

        //! Pointer to the Codimension 1 reference geometry, in this case, to ReferenceSquare.
        ReferenceGeometry<3>* const                         referenceGeometryCodim1Ptr_;

        //! Pointer to the Codimension 2 reference geometry, in this case, to ReferenceSquare.
        ReferenceGeometry<2>* const                         referenceGeometryCodim2Ptr_;

        //! Pointer to the Codimension 3 reference geometry, in this case, to ReferenceLine.
        ReferenceGeometry<1>* const                         referenceGeometryCodim3Ptr_;
        
        //! List of valid quadrature rules for this reference geometry
        std::list<QuadratureRules::GaussQuadratureRule<4>*> lstGaussQuadratureRules_;
    };
};
#endif
