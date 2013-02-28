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
    
    template<>
    const ReferenceGeometry<1>* const
    ElementGeometry<1>::CreateReferenceGeometry(unsigned int size)
    {
         std::cout <<"I am a line" << std::endl;
        return &ReferenceLine::Instance();
    }

    template<>
    const ReferenceGeometry<2>* const
    ElementGeometry<2>::CreateReferenceGeometry(unsigned int size)
    {
        if (size==4)
        {
            std::cout <<"I am a Ref square" << std::endl;
            return &ReferenceSquare::Instance();
        }
        else if (size==3)
        {
            std::cout <<"I am a triangle" << std::endl;
            return &ReferenceTriangle::Instance();
        }
        else
        {
            std::cout <<"WOW! Shit!!!!s" << std::endl;
            return NULL;
        }
    }
    
    template<>
    const ReferenceGeometry<3>* const
    ElementGeometry<3>::CreateReferenceGeometry(unsigned int size)
    {
         if (size==4)
        {
            std::cout <<"I am a tetrahedron" << std::endl;
            return &ReferenceTetrahedron::Instance();
        }
        else if (size==5)
        {
            std::cout <<"I am a pyramid" << std::endl;
            return &ReferencePyramid::Instance();
        }
        else if (size==6)
        {
            std::cout <<"I am a triangularPrism" << std::endl;
            return &ReferenceTriangularPrism::Instance();
        }
        else if (size==8)
        {
            std::cout <<"I am a cube" << std::endl;
            return &ReferenceCube::Instance();
        }
        else
        {
            std::cout <<"WOW! Shit!!!!s" << std::endl;
            return NULL;
        }
    }

    
    template<>
    const ReferenceGeometry<4>* const
    ElementGeometry<4>::CreateReferenceGeometry(unsigned int size)
    {
          if (size==16)
        {
            std::cout <<"I am a hypercube" << std::endl;
            return &ReferenceHypercube::Instance();
        }
        else
        {
            std::cout <<"WOW! Shit!!!!s" << std::endl;
            return NULL;
        }
    }
    template<>
    const PhysicalGeometry<1>* const
    ElementGeometry<1>::CreatePhysicalGeometry(const VectorOfPointIndexesT&      globalNodeIndexes,
                                                 const VectorOfPhysicalPointsT&    nodes,
                                                 const ReferenceGeometryT* const   geo)
    {
        std::cout <<"I am a line" << std::endl;
        return new Geometry::PhysicalLine(globalNodeIndexes, nodes, static_cast<const ReferenceLine* const>(geo));
    }
    template<>
    const PhysicalGeometry<2>* const
    ElementGeometry<2>::CreatePhysicalGeometry(const VectorOfPointIndexesT&      globalNodeIndexes,
                                                 const VectorOfPhysicalPointsT&    nodes,
                                                 const ReferenceGeometryT* const   geo)
    {
        cout << "globalNodeIndexes()="<<globalNodeIndexes.size()<<endl;
        if (globalNodeIndexes.size()==4)
        {
            std::cout <<"I am a physcial square" << std::endl;
            return new Geometry::PhysicalQuadrilateral(globalNodeIndexes, nodes, static_cast<const ReferenceSquare* const>(geo));
        }
        else if (globalNodeIndexes.size()==3)
        {
            std::cout <<"I am a triangle" << std::endl;
            return new Geometry::PhysicalTriangle(globalNodeIndexes, nodes, static_cast<const ReferenceTriangle* const>(geo));
        }
        else
        {
            std::cout <<"WOW! Shit!!!!s" << std::endl;
            return NULL;
        }
     }
    template<>
    const PhysicalGeometry<3>* const
    ElementGeometry<3>::CreatePhysicalGeometry(const VectorOfPointIndexesT&      globalNodeIndexes,
                                                 const VectorOfPhysicalPointsT&    nodes,
                                                 const ReferenceGeometryT* const   geo)
    {
      if (globalNodeIndexes.size()==4)
        {
            std::cout <<"I am a tetrahedron" << std::endl;
            return new Geometry::PhysicalTetrahedron(globalNodeIndexes, nodes, static_cast<const ReferenceCube* const>(geo));
        }
        else if (globalNodeIndexes.size()==5)
        {
            std::cout <<"I am a pyramid" << std::endl;
                //return new Geometry::PhysicalPyramid(globalNodeIndexes, nodes, static_cast<const ReferencePyramid* const>(geo));
            
        }
        else if (globalNodeIndexes.size()==6)
        {
            std::cout <<"I am a triangularPrism" << std::endl;
                //return new Geometry::PhysicalTriangularPrism(globalNodeIndexes, nodes, static_cast<const ReferenceTriangularPrism* const>(geo));
        }
        else if (globalNodeIndexes.size()==8)
        {
            std::cout <<"I am a cube" << std::endl;
            return new Geometry::PhysicalHexahedron(globalNodeIndexes, nodes, static_cast<const ReferenceCube* const>(geo));
        }
        else
        {
            std::cout <<"WOW! Shit!!!!s" << std::endl;
            return NULL;
        }
    }

    template<>
    const PhysicalGeometry<4>* const
    ElementGeometry<4>::CreatePhysicalGeometry(const VectorOfPointIndexesT&      globalNodeIndexes,
                                                 const VectorOfPhysicalPointsT&    nodes,
                                                 const ReferenceGeometryT* const   geo)
    {
        if (globalNodeIndexes.size()==16)
        {
            std::cout <<"I am a hypercube" << std::endl;
                // return new Geometry::PhysicalHypercube(globalNodeIndexes, nodes, static_cast<const ReferenceHypercube* const>(geo));
        }
        else
        {
            std::cout <<"WOW! Shit!!!!s" << std::endl;
            return NULL;
        }
    }
    
    template<>
    const MappingReferenceToPhysical<1, 1>* const
    ElementGeometry<1>::CreateMappings(unsigned int size, const PhysicalGeometryT* const pGeo)
    {
        std::cout <<"I am a line" << std::endl;
        return new Geometry::MappingToPhysHypercubeLinear<1>(pGeo);
    }
    
    template<>
    const MappingReferenceToPhysical<2, 2>* const
    ElementGeometry<2>::CreateMappings(unsigned int size, const PhysicalGeometryT* const pGeo)
    {
        if (size==4)
        {
            std::cout <<"I am a square" << std::endl;
            return new  Geometry::MappingToPhysHypercubeLinear<2>(pGeo);
        }
        else if (size==3)
        {
            std::cout <<"I am a triangle" << std::endl;
            return new  Geometry::MappingToPhysSimplexLinear<2>(pGeo);
        }
        else
        {
            std::cout <<"WOW! Shit!!!!s" << std::endl;
            return NULL;
        }
    }
    
    template<>
    const MappingReferenceToPhysical<3, 3>* const
    ElementGeometry<3>::CreateMappings(unsigned int size, const PhysicalGeometryT* const pGeo)
    {
        if (size==4)
        {
            std::cout <<"I am a tetrahedron" << std::endl;
                //return &ReferenceTetrahedron::Instance(); can not find this one
        }
        else if (size==5)
        {
            std::cout <<"I am a pyramid" << std::endl;
            return new  Geometry::MappingToPhysPyramid(pGeo);
        }
        else if (size==6)
        {
            std::cout <<"I am a triangularPrism" << std::endl;
            return new  Geometry::MappingToPhysTriangularPrism(pGeo);
        }
        else if (size==8)
        {
            std::cout <<"I am a cube" << std::endl;
            return new  Geometry::MappingToPhysHypercubeLinear<3>(pGeo);
        }
        else
        {
            std::cout <<"WOW! Shit!!!!s" << std::endl;
            return NULL;
        }
    }
    
    template<>
    const MappingReferenceToPhysical<4, 4>* const
    ElementGeometry<4>::CreateMappings(unsigned int size, const PhysicalGeometryT* const pGeo)
    {
        if (size==16)
        {
            std::cout <<"I am a hypercube" << std::endl;
            return new  Geometry::MappingToPhysHypercubeLinear<4>(pGeo);
        }
        else
        {
            std::cout <<"WOW! Shit!!!!s" << std::endl;
            return NULL;
        }
    }

    template<unsigned int DIM>
    ElementGeometry<DIM>::ElementGeometry(const VectorOfPointIndexesT&          globalNodeIndexes,
                                          const VectorOfPhysicalPointsT&        nodes):
        referenceGeometry_(ElementGeometry<DIM>::CreateReferenceGeometry(globalNodeIndexes.size())),
        physicalGeometry_(ElementGeometry<DIM>::CreatePhysicalGeometry(globalNodeIndexes, nodes, referenceGeometry_)),
        referenceToPhysicalMapping_(ElementGeometry<DIM>::CreateMappings(globalNodeIndexes.size(), physicalGeometry_)),
        refinementGeometry_(NULL)//refinement is turned off by default, to  enable it one needs to call enableRefinement
    {
    }
    
        /// Copy constructor
    template<unsigned int DIM>
    ElementGeometry<DIM>::ElementGeometry(const ElementGeometry& other):
        referenceGeometry_(other.referenceGeometry_),
        physicalGeometry_(ElementGeometry<DIM>::CreatePhysicalGeometry(other.physicalGeometry_->getNodeIndexes(), other.physicalGeometry_->getNodes(), referenceGeometry_)),
        referenceToPhysicalMapping_(ElementGeometry<DIM>::CreateMappings(other.physicalGeometry_->getNodeIndexes().size(), physicalGeometry_)),
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
