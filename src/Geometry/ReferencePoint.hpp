//
//  ReferencePoint.h
//  
//
//  Created by Shavarsh Nurijanyan on 2/4/13.
//
//

#ifndef ____ReferencePoint__
#define ____ReferencePoint__

#include <iostream>
using std::ostream;

#include "PointReference.hpp"
#include "ReferenceGeometry.hpp"
#include "GlobalNamespaceGeometry.hpp"
#include "Mappings/MappingReferenceToReference.hpp"
#include "Integration/QuadratureRules/GaussQuadratureRule.hpp"

namespace Geometry
{
    /// \class ReferencePoint
    /// \brief Reference geometry of dimensions 0.
    /// \details
    /// Most of the codimension functions are not well defined for dimension 0, so be careful.
    class ReferencePoint : public ReferenceGeometry
    {
    public:
        typedef ReferenceGeometry ReferenceGeometryT;
        
        using ReferenceGeometryT::PointReferenceT;
        using ReferenceGeometryT::VectorOfReferencePointsT;
        using ReferenceGeometryT::IndexT;
        typedef std::vector<IndexT> ListOfIndexesT;

    public:
        static ReferencePoint& Instance()
        {
            static ReferencePoint theInstance;
            return theInstance;
        }

        ReferencePoint();

        ReferencePoint(const ReferencePoint&);
        
        /// \brief Return true. TODO: Why?
        bool                                    isInternalPoint(const PointReferenceT& p) const;

        /// \brief (see ReferenceGeometry.hpp)
        void                                    getCenter(PointReferenceT&) const;

        /// \brief (see ReferenceGeometry.hpp)
        void                                    getNode(const IndexT& i, PointReferenceT& point) const;

        /// \brief (see ReferenceGeometry.hpp)
        String                                  getName() const {return "ReferencePoint";}

        /// \brief Given a face index, and an index of the node position relative to the face,
        /// return the local index of the node. Dummy function, doesn't make sense for point
        int                                     getLocalNodeIndex(int face, int node) const{return -1;}
        
        // ================================== Codimension 0 ========================================

        /// \brief Returns 0.
        int                                         getCodim0MappingIndex(const ListOfIndexesT&, const ListOfIndexesT&) const;

        /// \brief Returns 0.
        const MappingReferenceToReference*     getCodim0MappingPtr(const IndexT a) const;
        
        // ================================== Quadrature rules =====================================

        /// \brief Add a quadrature rule into the list of valid quadrature rules for this geometry.
        virtual void addGaussQuadratureRule(QuadratureRules::GaussQuadratureRule* const qr);

        /// \brief Get a valid quadrature for this geometry.
        virtual QuadratureRules::GaussQuadratureRule* const getGaussQuadratureRule(int order) const;

        // =============================== Refinement mappings =====================================
        
        /// \brief Transform a reference point using refinement mapping
        void refinementTransform(int refineType, int subElementIdx, 
                      const PointReferenceT& p, PointReferenceT& pMap) const {}

        /// \brief Transformation matrix of this refinement when located on the LEFT side
        void getRefinementMappingMatrixL(int refineType, int subElementIdx, 
                    LinearAlgebra::Matrix& Q) const {}

        /// \brief Transformation matrix of this refinement when located on the RIGHT side
        void getRefinementMappingMatrixR(int refineType, int subElementIdx, 
                    LinearAlgebra::Matrix& Q) const {}

        /// \brief Refinement mapping on codim1 for a given refinement on codim0
        /// Note: this should also applied on other dimensions
        void getCodim1RefinementMappingMatrixL(int refineType, DimT subElementIdx, 
                                DimT faLocalIndex, LinearAlgebra::Matrix& Q) const {}

        /// \brief Refinement mapping on codim1 for a given refinement on codim0
        /// \brief Note: this should also applied on other dimensions
        void getCodim1RefinementMappingMatrixR(int refineType, DimT subElementIdx, 
                                DimT faLocalIndex, LinearAlgebra::Matrix& Q) const {}

        /// \brief List of valid quadrature rules for this reference geometry
        std::list<QuadratureRules::GaussQuadratureRule*> lstGaussQuadratureRules_;

        //! Codimension 1 mappings, from a line to a line. TODO: Where is this used? clarify here.
        const MappingReferenceToReference*                           mappingsPointToPoint_; //!< codim0
    };
};
#endif
