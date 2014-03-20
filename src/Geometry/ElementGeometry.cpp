//
//
//  ElementGeometry.cpp
//
//
//  Created by Shavarsh Nurijanyan on 7/2/13.
//
//

#include <ElementGeometry.hpp>


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
    class ElementGeometry;

    const ReferenceGeometry* const
    ElementGeometry::createReferenceGeometry(unsigned int size,unsigned int DIM)
    {///\todo check for consistency of pairs (size,DIM)
    	switch(size){//select a proper type based on the number of nodes a reference geometry should have
    	case 2:
//        std::cout <<"I am a line" << std::endl;
    		return &ReferenceLine::Instance();
    	case 3:
//            std::cout <<"I am a triangle" << std::endl;
            return &ReferenceTriangle::Instance();
    	case 4:
    		if(DIM==2){
	//            std::cout <<"I am a Ref square" << std::endl;
				return &ReferenceSquare::Instance();
    		}else if(DIM==3){
	//            std::cout <<"I am a tetrahedron" << std::endl;
				return &ReferenceTetrahedron::Instance();
    		}else{
    			throw "This DIMension does not contain entities with 4 nodes";
    		}
    	case 5:
//            std::cout <<"I am a pyramid" << std::endl;
            return &ReferencePyramid::Instance();
    	case 6:
//            std::cout <<"I am a triangularPrism" << std::endl;
            return &ReferenceTriangularPrism::Instance();
    	case 8:
//            std::cout <<"I am a cube" << std::endl;
            return &ReferenceCube::Instance();
    	case 16:
//            std::cout <<"I am a hypercube" << std::endl;
            return &ReferenceHypercube::Instance();
    	default:
    		throw "No known entities contain this many nodes";
    	}
    }

    const PhysicalGeometry* const
    ElementGeometry::createPhysicalGeometry(const VectorOfPointIndexesT&      globalNodeIndexes,
                                               const VectorOfPhysicalPointsT&    nodes,
                                               const ReferenceGeometryT* const   geo)
    {
    	/*switch(typeid(*geo)){
    	case typeid(ReferenceLine::Instance()):
    			...
    	case ...
    	}*/
    	switch(globalNodeIndexes.size()){
    	case 2:
	//        std::cout <<"I am a line" << std::endl;
			return new Geometry::PhysicalLine(globalNodeIndexes, nodes, static_cast<const ReferenceLine* const>(geo));
    	case 3:
//            std::cout <<"I am a triangle" << std::endl;
            return new Geometry::PhysicalTriangle(globalNodeIndexes, nodes, static_cast<const ReferenceTriangle* const>(geo));
    	case 4:
    		if(nodes[0].size()==2){
	//            std::cout <<"I am a physcial square" << std::endl;
				return new Geometry::PhysicalQuadrilateral(globalNodeIndexes, nodes, static_cast<const ReferenceSquare* const>(geo));
    		}else if(nodes[0].size()==3){
	//            std::cout <<"I am a tetrahedron" << std::endl;
				return new Geometry::PhysicalTetrahedron(globalNodeIndexes, nodes, static_cast<const ReferenceTetrahedron* const>(geo));
    		}else{
    			throw "This DIMension does not contain entities with 4 nodes";
    		}
    	case 5:
//            std::cout <<"I am a pyramid" << std::endl;
            return new Geometry::PhysicalPyramid(globalNodeIndexes, nodes, static_cast<const ReferencePyramid* const>(geo));
    	case 6:
//            std::cout <<"I am a triangularPrism" << std::endl;
            return new Geometry::PhysicalTriangularPrism(globalNodeIndexes, nodes, static_cast<const ReferenceTriangularPrism* const>(geo));
    	case 8:
//            std::cout <<"I am a cube" << std::endl;
            return new Geometry::PhysicalHexahedron(globalNodeIndexes, nodes, static_cast<const ReferenceCube* const>(geo));
    	case 16:
//            std::cout <<"I am a hypercube" << std::endl;
    		///\todo can not find this one
            //return new Geometry::PhysicalHypercube(globalNodeIndexes, nodes, static_cast<const ReferenceHypercube* const>(geo));
    	default:
    		throw "No known entities contain this many nodes";

    	}
    }

    const MappingReferenceToPhysical* const
    ElementGeometry::createMappings(unsigned int size,unsigned int DIM, const PhysicalGeometryT* const pGeo)
    {
    	switch(size){
    	case 2:
//			std::cout <<"I am a line" << std::endl;
			return new Geometry::MappingToPhysHypercubeLinear<1>(pGeo);
    	case 3:
//            std::cout <<"I am a triangle" << std::endl;
            return new  Geometry::MappingToPhysSimplexLinear<2>(pGeo);
    	case 4:
    		if(DIM==2){
//				std::cout <<"I am a square" << std::endl;
				return new  Geometry::MappingToPhysHypercubeLinear<2>(pGeo);
    		}else if(DIM==3){
	//            std::cout <<"I am a tetrahedron" << std::
				return new  Geometry::MappingToPhysSimplexLinear<3>(pGeo);
    		}else{
    			throw "This DIMension does not contain entities with 4 nodes";
    		}
    	case 5:
//            std::cout <<"I am a pyramid" << std::endl;
            return new  Geometry::MappingToPhysPyramid(pGeo);
    	case 6:
//            std::cout <<"I am a triangularPrism" << std::endl;
            return new  Geometry::MappingToPhysTriangularPrism(pGeo);
    	case 8:
//            std::cout <<"I am a cube" << std::endl;
            return new  Geometry::MappingToPhysHypercubeLinear<3>(pGeo);
    	case 16:
//            std::cout <<"I am a hypercube" << std::endl;
            return new  Geometry::MappingToPhysHypercubeLinear<4>(pGeo);
    	default:
    		throw "No known entities contain this many nodes";
    	}
    }

    ElementGeometry::ElementGeometry(const VectorOfPointIndexesT&          globalNodeIndexes,
                                          const VectorOfPhysicalPointsT&        nodes):
        referenceGeometry_(ElementGeometry::createReferenceGeometry(globalNodeIndexes.size(),nodes[0].size())),
        physicalGeometry_(ElementGeometry::createPhysicalGeometry(globalNodeIndexes, nodes, referenceGeometry_)),
        referenceToPhysicalMapping_(ElementGeometry::createMappings(globalNodeIndexes.size(),nodes[0].size(), physicalGeometry_)),
        refinementGeometry_(NULL)//refinement is turned off by default, to  enable it one needs to call enableRefinement
    {
    }
    
        /// Copy constructor
    ElementGeometry::ElementGeometry(const ElementGeometry& other):
        referenceGeometry_(other.referenceGeometry_),
        physicalGeometry_(ElementGeometry::createPhysicalGeometry(other.physicalGeometry_->getNodeIndexes(), other.physicalGeometry_->getNodes(), referenceGeometry_)),
        referenceToPhysicalMapping_(ElementGeometry::createMappings(other.physicalGeometry_->getNodeIndexes().size(),other.physicalGeometry_->getNodePtr(0)->size(), physicalGeometry_)),
        refinementGeometry_(other.refinementGeometry_)//refinement is turned off by default, to  enable it one needs to call enableRefinement
    {

    }
    
    ElementGeometry::~ElementGeometry()
    {
        delete physicalGeometry_;
        delete referenceToPhysicalMapping_;
    }

        /// Returns a pointer to the referenceToPhysicalMapping
    const MappingReferenceToPhysical* const
    ElementGeometry::getReferenceToPhysicalMap() const
    {
        return referenceToPhysicalMapping_;
    }
    
        /// Returns a pointer to the physicalGeometry object.
    const PhysicalGeometry* const
    ElementGeometry::getPhysicalGeometry() const
    {
        return physicalGeometry_;
    }
    
        /// Returns a pointer to the referenceGeometry object.
    const ReferenceGeometry* const
    ElementGeometry::getReferenceGeometry() const
    {
        return referenceGeometry_;
    }
    
        /// Returns a pointer to the refinementGeometry object.
    const RefinementGeometry*
    ElementGeometry::getRefinementGeometry() const
    {
        return refinementGeometry_;
    }
    
        /// This method gets a PointReference, which specifies a coordinate in the ReferenceGeometry,
        /// and returns a PointPhysical which is the corresponding point in the PhysicalGeometry,
        /// given the mapping.
    void
    ElementGeometry::referenceToPhysical(const PointReferenceT& pointReference,
                             PointPhysicalT& pointPhysical)
    {
        referenceToPhysicalMapping_->transform(pointReference, pointPhysical);
    }
    
    void
    ElementGeometry::referenceToPhysical(const PointReferenceT& pointReference,
                                              PointPhysicalT& pointPhysical)const
    {
        referenceToPhysicalMapping_->transform(pointReference, pointPhysical);
    }
    
        /// This method gets a PointReference and returns the corresponding jacobian of the
        /// referenceToPhysicalMapping.

    void
    ElementGeometry::calcJacobian(const PointReferenceT& pointReference, JacobianT& jacobian) const
    {
        referenceToPhysicalMapping_->calcJacobian(pointReference,jacobian);
    }

    unsigned int
    ElementGeometry::getNrOfNodes() const
    {
    	return physicalGeometry_->getNumberOfNodes();
    }
    
}

#endif
