//
//  ElementGeometry_Impl.hpp
//  
//
//  Created by Shavarsh Nurijanyan on 2/24/13.
//
//

#ifndef _ElementGeometry_Impl_hpp
#define _ElementGeometry_Impl_hpp

#include "ReferenceTetrahedron.hpp"
#include "ReferenceLine.hpp"
#include "ReferenceSquare.hpp"
#include "ReferenceTriangle.hpp"
#include "ReferencePyramid.hpp"
#include "ReferenceTriangularPrism.hpp"
#include "ReferenceCube.hpp"
#include "ReferenceHypercube.hpp"

#include "PhysicalTetrahedron.hpp"
#include "PhysicalLine.hpp"
#include "PhysicalQuadrilateral.hpp"
#include "PhysicalTriangle.hpp"
#include "PhysicalPyramid.hpp"
#include "PhysicalTriangularPrism.hpp"
#include "PhysicalHexahedron.hpp"
    //#include "PhysicalOcot...hpp"

#include "Mappings/MappingReferenceToPhysical.hpp"
#include "Mappings/MappingToPhysHypercubeLinear.hpp"
#include "Mappings/MappingToPhysSimplexLinear.hpp"
#include "Mappings/MappingToPhysPyramid.hpp"
#include "Mappings/MappingToPhysTriangularPrism.hpp"

#include "RefinementLine.hpp"
#include "RefinementTriangle.hpp"
#include "RefinementQuadrilateral.hpp"
#include "RefinementTetrahedron.hpp"
#include "RefinementPyramid.hpp"
#include "RefinementTriangularPrism.hpp"
#include "RefinementHexahedron.hpp"
#include "RefinementHypercube.hpp"




namespace Geometry
{
    template<unsigned int DIM>
    class ElementGeometry;

    template<unsigned int DIM>
    ElementGeometry<DIM>::ElementGeometry(const VectorOfPointIndexesT&          globalNodeIndexes,
                                          const VectorOfPhysicalPointsT&        nodes):
        referenceGeometry_(ElementGeometry<DIM>::createReferenceGeometry(globalNodeIndexes.size())),
        physicalGeometry_(ElementGeometry<DIM>::createPhysicalGeometry(globalNodeIndexes, nodes, referenceGeometry_)),
        referenceToPhysicalMapping_(ElementGeometry<DIM>::createMappings(globalNodeIndexes.size(), physicalGeometry_)),
        refinementGeometry_(NULL)//refinement is turned off by default, to  enable it one needs to call enableRefinement
    {
    }
    
        /// Copy constructor
    template<unsigned int DIM>
    ElementGeometry<DIM>::ElementGeometry(const ElementGeometry& other):
        referenceGeometry_(other.referenceGeometry_),
        physicalGeometry_(ElementGeometry<DIM>::createPhysicalGeometry(other.physicalGeometry_->getNodeIndexes(), other.physicalGeometry_->getNodes(), referenceGeometry_)),
        referenceToPhysicalMapping_(ElementGeometry<DIM>::createMappings(other.physicalGeometry_->getNodeIndexes().size(), physicalGeometry_)),
        refinementGeometry_(other.refinementGeometry_)//refinement is turned off by default, to  enable it one needs to call enableRefinement
    {

    }
    
    template<unsigned int DIM>
    ElementGeometry<DIM>::~ElementGeometry()
    {
        delete physicalGeometry_;
        delete referenceToPhysicalMapping_;
    }

        /// Returns a pointer to the referenceToPhysicalMapping
    template<unsigned int DIM>
    const MappingReferenceToPhysical<DIM,DIM>* const
    ElementGeometry<DIM>::getReferenceToPhysicalMap() const
    {
        return referenceToPhysicalMapping_;
    }
    
        /// Returns a pointer to the physicalGeometry object.
    template<unsigned int DIM>
    const PhysicalGeometry<DIM>* const
    ElementGeometry<DIM>::getPhysicalGeometry() const
    {
        return physicalGeometry_;
    }
    
        /// Returns a pointer to the referenceGeometry object.
    template<unsigned int DIM>
    const ReferenceGeometry<DIM>* const
    ElementGeometry<DIM>::getReferenceGeometry() const
    {
        return referenceGeometry_;
    }
    
        /// Returns a pointer to the refinementGeometry object.
    template<unsigned int DIM>
    const RefinementGeometry<DIM>* 
    ElementGeometry<DIM>::getRefinementGeometry() const
    {
        return refinementGeometry_;
    }
    
        /// This method gets a PointReference, which specifies a coordinate in the ReferenceGeometry,
        /// and returns a PointPhysical which is the corresponding point in the PhysicalGeometry,
        /// given the mapping.
    template<unsigned int DIM>
    void
    ElementGeometry<DIM>::referenceToPhysical(const PointReferenceT& pointReference,
                             PointPhysicalT& pointPhysical)
    {
        referenceToPhysicalMapping_->transform(pointReference, pointPhysical);
    }
    
    template<unsigned int DIM>
    void
    ElementGeometry<DIM>::referenceToPhysical(const PointReferenceT& pointReference,
                                              PointPhysicalT& pointPhysical)const
    {
        referenceToPhysicalMapping_->transform(pointReference, pointPhysical);
    }
    
        /// This method gets a PointReference and returns the corresponding jacobian of the
        /// referenceToPhysicalMapping.

    template<unsigned int DIM>
    void
    ElementGeometry<DIM>::calcJacobian(const PointReferenceT& pointReference, JacobianT& jacobian) const
    {
        referenceToPhysicalMapping_->calcJacobian(pointReference,jacobian);
    }

    
}

#endif
