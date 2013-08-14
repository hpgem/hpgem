//
//  ElementGeometry.cpp
//  
//
//  Created by Shavarsh Nurijanyan on 7/2/13.
//
//

#include <ElementGeometry.hpp>

namespace Geometry
{
    template<unsigned int DIM>
    class ElementGeometry;
    
    template<>
    const ReferenceGeometry<1>* const
    ElementGeometry<1>::createReferenceGeometry(unsigned int size)
    {
//        std::cout <<"I am a line" << std::endl;
        return &ReferenceLine::Instance();
    }
    
    template<>
    const ReferenceGeometry<2>* const
    ElementGeometry<2>::createReferenceGeometry(unsigned int size)
    {
        if (size==4)
        {
//            std::cout <<"I am a Ref square" << std::endl;
            return &ReferenceSquare::Instance();
        }
        else if (size==3)
        {
//            std::cout <<"I am a triangle" << std::endl;
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
    ElementGeometry<3>::createReferenceGeometry(unsigned int size)
    {
        if (size==4)
        {
//            std::cout <<"I am a tetrahedron" << std::endl;
            return &ReferenceTetrahedron::Instance();
        }
        else if (size==5)
        {
//            std::cout <<"I am a pyramid" << std::endl;
            return &ReferencePyramid::Instance();
        }
        else if (size==6)
        {
//            std::cout <<"I am a triangularPrism" << std::endl;
            return &ReferenceTriangularPrism::Instance();
        }
        else if (size==8)
        {
//            std::cout <<"I am a cube" << std::endl;
            return &ReferenceCube::Instance();
        }
        else
        {
//            std::cout <<"WOW! Shit!!!!s" << std::endl;
            return NULL;
        }
    }
    
    
    template<>
    const ReferenceGeometry<4>* const
    ElementGeometry<4>::createReferenceGeometry(unsigned int size)
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
    ElementGeometry<1>::createPhysicalGeometry(const VectorOfPointIndexesT&      globalNodeIndexes,
                                               const VectorOfPhysicalPointsT&    nodes,
                                               const ReferenceGeometryT* const   geo)
    {
        std::cout <<"I am a line" << std::endl;
        return new Geometry::PhysicalLine(globalNodeIndexes, nodes, static_cast<const ReferenceLine* const>(geo));
    }
    template<>
    const PhysicalGeometry<2>* const
    ElementGeometry<2>::createPhysicalGeometry(const VectorOfPointIndexesT&      globalNodeIndexes,
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
    ElementGeometry<3>::createPhysicalGeometry(const VectorOfPointIndexesT&      globalNodeIndexes,
                                               const VectorOfPhysicalPointsT&    nodes,
                                               const ReferenceGeometryT* const   geo)
    {
        if (globalNodeIndexes.size()==4)
        {
//            std::cout <<"I am a tetrahedron" << std::endl;
            return new Geometry::PhysicalTetrahedron(globalNodeIndexes, nodes, static_cast<const ReferenceCube* const>(geo));
        }
        else if (globalNodeIndexes.size()==5)
        {
//            std::cout <<"I am a pyramid" << std::endl;
                //return new Geometry::PhysicalPyramid(globalNodeIndexes, nodes, static_cast<const ReferencePyramid* const>(geo));
            
        }
        else if (globalNodeIndexes.size()==6)
        {
//            std::cout <<"I am a triangularPrism" << std::endl;
                //return new Geometry::PhysicalTriangularPrism(globalNodeIndexes, nodes, static_cast<const ReferenceTriangularPrism* const>(geo));
        }
        else if (globalNodeIndexes.size()==8)
        {
//            std::cout <<"I am a cube" << std::endl;
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
    ElementGeometry<4>::createPhysicalGeometry(const VectorOfPointIndexesT&      globalNodeIndexes,
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
    ElementGeometry<1>::createMappings(unsigned int size, const PhysicalGeometryT* const pGeo)
    {
        std::cout <<"I am a line" << std::endl;
        return new Geometry::MappingToPhysHypercubeLinear<1>(pGeo);
    }
    
    template<>
    const MappingReferenceToPhysical<2, 2>* const
    ElementGeometry<2>::createMappings(unsigned int size, const PhysicalGeometryT* const pGeo)
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
    ElementGeometry<3>::createMappings(unsigned int size, const PhysicalGeometryT* const pGeo)
    {
        if (size==4)
        {
//            std::cout <<"I am a tetrahedron" << std::endl;
                //return &ReferenceTetrahedron::Instance(); can not find this one
        }
        else if (size==5)
        {
//            std::cout <<"I am a pyramid" << std::endl;
            return new  Geometry::MappingToPhysPyramid(pGeo);
        }
        else if (size==6)
        {
//            std::cout <<"I am a triangularPrism" << std::endl;
            return new  Geometry::MappingToPhysTriangularPrism(pGeo);
        }
        else if (size==8)
        {
//            std::cout <<"I am a cube" << std::endl;
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
    ElementGeometry<4>::createMappings(unsigned int size, const PhysicalGeometryT* const pGeo)
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
}