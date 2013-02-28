//
//  ReferenceCube.h
//
//
//  Created by nicorivas on 2/16/13.
//
//

#ifndef REFERENCECUBE_HH
#define REFERENCECUBE_HH

#include "GlobalNamespaceGeometry.hpp"
#include "ReferenceGeometry.hpp"
#include "Mappings/MappingReferenceToReference.hpp"
#include "Integration/QuadratureRules/GaussQuadratureRule.hpp"

namespace Geometry
{
    /*! The ordering of the vertex and faces in a cube:
     *
     *     6o---------o7
     *     /|        /|
     *    / |       / |
     *  4o---------o5 |
     *   | 2o------|--o3
     *   | /       | /
     *   |/        |/
     *  0o---------o1
     *
     */
    class ReferenceCube: public ReferenceGeometry<ThreeD>
    {
    public:
        typedef ReferenceGeometry<ThreeD> ReferenceGeometryT;

        using typename ReferenceGeometryT::IndexT;
        using typename ReferenceGeometryT::PointReferenceT;
        using typename ReferenceGeometryT::VectorOfReferencePointsT;
        using ReferenceGeometryT::String;

    public:

        static ReferenceCube& Instance()
        {
            static ReferenceCube theInstance;
            return theInstance;
        }

        ReferenceCube();

        ReferenceCube(const ReferenceCube& other);

        //! (see ReferenceGeometry.hpp)
        bool        isInternalPoint(const PointReferenceT& p) const;

        //! (see ReferenceGeometry.hpp)
        void        getCenter(PointReferenceT& p) const;

        //! (see ReferenceGeometry.hpp)
        void        getNode(const IndexT& i, PointReferenceT& point) const;

        //! (see ReferenceGeometry.hpp)
        String      getName() const {return "ReferenceCube";}

        //! Given a face index, and an index of the node position relative to the face,
        //! return the local index of the node.
        int         getLocalNodeIndex(int face, int node) const
                    {return localNodeIndexes_[face][node];}

        //! Output routine
        friend std::ostream& operator<<(std::ostream& os, const ReferenceCube& cube);

        // ================================== Codimension 0 ========================================

        //! (see MappingCodimensions.hpp)
        int getCodim0MappingIndex(const ListOfIndexesT&, const ListOfIndexesT&) const;

        //! (see MappingCodimensions.hpp)
        const MappingReferenceToReference<3, 3>* getCodim0MappingPtr(const IndexT) const;

        // ================================== Codimension 1 ========================================

        //! (see MappingCodimensions.hpp)
        unsigned int                                getNrOfCodim1Entities() const {return 6;}

        //! (see MappingCodimensions.hpp)
        void                                        getCodim1EntityLocalIndices(const IndexT, ListOfIndexesT& faceNodesLocal) const;

        //! (see MappingCodimensions.hpp)
        const MappingReferenceToReference<2, 3>*    getCodim1MappingPtr(const IndexT) const;

        //! (see MappingCodimensions.hpp)
        const ReferenceGeometry<2>*                 getCodim1ReferenceGeometry(const IndexT) const;

        // ================================== Codimension 2 ========================================

        //! (see MappingCodimensions.hpp)
        unsigned int getNrOfCodim2Entities() const {return 12;}

        //! (see MappingCodimensions.hpp)
        void                                        getCodim2EntityLocalIndices(const IndexT, ListOfIndexesT& faceNodesLocal) const;

        //! (see MappingCodimensions.hpp)
        const MappingReferenceToReference<1, 3>*    getCodim2MappingPtr(const IndexT) const;

        //! (see MappingCodimensions.hpp)
        const ReferenceGeometry<1>*                 getCodim2ReferenceGeometry(const IndexT) const;

        // ================================== Codimension 3 ========================================

        //! (see MappingCodimensions.hpp)
        unsigned int    getNrOfCodim3Entities() const {return 8;};

        //! (see MappingCodimensions.hpp)
        void            getCodim3EntityLocalIndices(const unsigned int, std::vector<unsigned int>&) const {return;}

        // ================================== Quadrature rules =====================================

        /// Add a quadrature rule into the list of valid quadrature rules for this geometry.
        virtual void addGaussQuadratureRule(QuadratureRules::GaussQuadratureRule<3>* const qr); 

        /// Get a valid quadrature for this geometry.
        virtual QuadratureRules::GaussQuadratureRule<3>* const getGaussQuadratureRule(int order) const;

        // =============================== Refinement mappings =====================================
        
        //! Transform a reference point using refinement mapping
        virtual void refinementTransform(int refineType, int subElementIdx, 
                      const PointReferenceT& p, PointReferenceT& pMap) const;

        //! Transformation matrix of this refinement when located on the LEFT side
        virtual void getRefinementMappingMatrixL(int refineType, int subElementIdx, 
                    LinearAlgebra::Matrix& Q) const;

        //! Transformation matrix of this refinement when located on the RIGHT side
        virtual void getRefinementMappingMatrixR(int refineType, int subElementIdx, 
                    LinearAlgebra::Matrix& Q) const;

        //! Refinement mapping on codim1 for a given refinement on codim0
        virtual void getCodim1RefinementMappingMatrixL(int refineType, DimT subElementIdx, 
                                DimT faLocalIndex, LinearAlgebra::Matrix& Q) const;

        //! Refinement mapping on codim1 for a given refinement on codim0
        virtual void getCodim1RefinementMappingMatrixR(int refineType, DimT subElementIdx, 
                                DimT faLocalIndex, LinearAlgebra::Matrix& Q) const;

      private:
        //! Local node indexes contains the numbering of the vertex of the shape, ordered by faces.
        //! See top comment for the corresponding numbering.
        static int                                          localNodeIndexes_[6][4];

        //! The nodes on edge contains the local index of the two nodes in every edge.
        static int                                          localNodesOnEdge_[12][2]; //!< 12 edges with 2 nodes

        //! Codimension 1 mappings, from a line to a square. TODO: Where is this used? clarify here.
        const MappingReferenceToReference<2, 3>*            mappingsSquareToCube_[6];

        //! Codimension 0 mappings, from a square to a square. TODO: Where is this used? clarifiy here.
        const MappingReferenceToReference<3, 3>*            mappingsCubeToCube_[8];

        //! Pointer to the Codimension 1 reference geometry, in this case, to ReferenceSquare.
        ReferenceGeometry<2>* const                         referenceGeometryCodim1Ptr_;

        //! Pointer to the Codimension 2 reference geometry, in this case, to ReferenceLine.
        ReferenceGeometry<1>* const                         referenceGeometryCodim2Ptr_;
        
        //! List of valid quadrature rules for this reference geometry
        std::list<QuadratureRules::GaussQuadratureRule<3>*> lstGaussQuadratureRules_;
    };
};
#endif
