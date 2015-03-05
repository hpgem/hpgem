/*
 This file forms part of hpGEM. This package has been developed over a number of years by various people at the University of Twente and a full list of contributors can be found at
 http://hpgem.org/about-the-code/team
 
 This code is distributed using BSD 3-Clause License. A copy of which can found below.
 
 
 Copyright (c) 2014, University of Twente
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 
 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef _ElementGeometry_Impl_hpp
#define _ElementGeometry_Impl_hpp

#include <ElementGeometry.hpp>

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
#include "PhysicalOctachoron.hpp"

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

#include "PointReference.hpp"



namespace Geometry
{
    class ElementGeometry;

    const ReferenceGeometry * const
    ElementGeometry::createReferenceGeometry(std::size_t size, std::size_t DIM)
    {///\todo check for consistency of pairs (size,DIM)
        switch (size)
        {//select a proper type based on the number of nodes a reference geometry should have
            case 2:
                //        std::cout <<"I am a line" << std::endl;
                return &ReferenceLine::Instance();
            case 3:
                //            std::cout <<"I am a triangle" << std::endl;
                return &ReferenceTriangle::Instance();
            case 4:
                if (DIM == 2)
                {
                    //            std::cout <<"I am a Ref square" << std::endl;
                    return &ReferenceSquare::Instance();
                }
                else if (DIM == 3)
                {
                    //            std::cout <<"I am a tetrahedron" << std::endl;
                    return &ReferenceTetrahedron::Instance();
                }
                else
                {
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

    const PhysicalGeometry * const
    ElementGeometry::createPhysicalGeometry(const VectorOfPointIndexesT& globalNodeIndexes,
                                            const VectorOfPhysicalPointsT& nodes,
                                            const ReferenceGeometryT * const geo)
    {
        /*switch(typeid(*geo)){
        case typeid(ReferenceLine::Instance()):
                ...
        case ...
        }*/
        switch (globalNodeIndexes.size())
        {
            case 2:
                //        std::cout <<"I am a line" << std::endl;
                return new Geometry::PhysicalLine(globalNodeIndexes, nodes);
            case 3:
                //            std::cout <<"I am a triangle" << std::endl;
                return new Geometry::PhysicalTriangle(globalNodeIndexes, nodes);
            case 4:
                if (nodes[0].size() == 2)
                {
                    //            std::cout <<"I am a physcial square" << std::endl;
                    return new Geometry::PhysicalQuadrilateral(globalNodeIndexes, nodes);
                }
                else if (nodes[0].size() == 3)
                {
                    //            std::cout <<"I am a tetrahedron" << std::endl;
                    return new Geometry::PhysicalTetrahedron(globalNodeIndexes, nodes);
                }
                else
                {
                    throw "This DIMension does not contain entities with 4 nodes";
                }
            case 5:
                //            std::cout <<"I am a pyramid" << std::endl;
                return new Geometry::PhysicalPyramid(globalNodeIndexes, nodes);
            case 6:
                //            std::cout <<"I am a triangularPrism" << std::endl;
                return new Geometry::PhysicalTriangularPrism(globalNodeIndexes, nodes);
            case 8:
                //            std::cout <<"I am a cube" << std::endl;
                return new Geometry::PhysicalHexahedron(globalNodeIndexes, nodes);
            case 16:
                //            std::cout <<"I am a hypercube" << std::endl;
                return new Geometry::PhysicalOctachoron(globalNodeIndexes, nodes);
            default:
                throw "No known entities contain this many nodes";

        }
    }

    const MappingReferenceToPhysical * const
    ElementGeometry::createMappings(std::size_t size, std::size_t DIM, const PhysicalGeometryT * const pGeo)
    {
        switch (size)
        {
            case 2:
                //			std::cout <<"I am a line" << std::endl;
                return new Geometry::MappingToPhysHypercubeLinear<1>(pGeo);
            case 3:
                //            std::cout <<"I am a triangle" << std::endl;
                return new Geometry::MappingToPhysSimplexLinear<2>(pGeo);
            case 4:
                if (DIM == 2)
                {
                    //				std::cout <<"I am a square" << std::endl;
                    return new Geometry::MappingToPhysHypercubeLinear<2>(pGeo);
                }
                else if (DIM == 3)
                {
                    //            std::cout <<"I am a tetrahedron" << std::
                    return new Geometry::MappingToPhysSimplexLinear<3>(pGeo);
                }
                else
                {
                    throw "This DIMension does not contain entities with 4 nodes";
                }
            case 5:
                //            std::cout <<"I am a pyramid" << std::endl;
                return new Geometry::MappingToPhysPyramid(pGeo);
            case 6:
                //            std::cout <<"I am a triangularPrism" << std::endl;
                return new Geometry::MappingToPhysTriangularPrism(pGeo);
            case 8:
                //            std::cout <<"I am a cube" << std::endl;
                return new Geometry::MappingToPhysHypercubeLinear<3>(pGeo);
            case 16:
                //            std::cout <<"I am a hypercube" << std::endl;
                return new Geometry::MappingToPhysHypercubeLinear<4>(pGeo);
            default:
                throw "No known entities contain this many nodes";
        }
    }

    std::ostream& operator<<(std::ostream& os, const ElementGeometry& elementGeometry)
    {
        os << "PhysicalGeometry={";
        for (std::size_t i = 0; i < elementGeometry.physicalGeometry_->getNumberOfNodes(); i++)
        {
            os << (elementGeometry.physicalGeometry_)->getNodeIndex(i) << " ";
        }
        os << '}' << std::endl;
        return os;
    }

    ElementGeometry::ElementGeometry(const VectorOfPointIndexesT& globalNodeIndexes,
                                     const VectorOfPhysicalPointsT& nodes) :
    referenceGeometry_(ElementGeometry::createReferenceGeometry(globalNodeIndexes.size(), nodes[0].size())),
    physicalGeometry_(ElementGeometry::createPhysicalGeometry(globalNodeIndexes, nodes, referenceGeometry_)),
    referenceToPhysicalMapping_(ElementGeometry::createMappings(globalNodeIndexes.size(), nodes[0].size(), physicalGeometry_)),
    refinementGeometry_(nullptr)//refinement is turned off by default, to  enable it one needs to call enableRefinement
    { }

    /// Copy constructor

    ElementGeometry::ElementGeometry(const ElementGeometry& other) :
    referenceGeometry_(other.referenceGeometry_),
    physicalGeometry_(ElementGeometry::createPhysicalGeometry(other.physicalGeometry_->getNodeIndexes(), other.physicalGeometry_->getNodes(), referenceGeometry_)),
    referenceToPhysicalMapping_(ElementGeometry::createMappings(other.physicalGeometry_->getNodeIndexes().size(), other.physicalGeometry_->getNodePtr(0)->size(), physicalGeometry_)),
    refinementGeometry_(other.refinementGeometry_)//refinement is turned off by default, to  enable it one needs to call enableRefinement
    { }

    ElementGeometry::~ElementGeometry()
    {
        delete physicalGeometry_;
        delete referenceToPhysicalMapping_;
    }

    /// Returns a pointer to the referenceToPhysicalMapping

    const MappingReferenceToPhysical * const
    ElementGeometry::getReferenceToPhysicalMap() const
    {
        return referenceToPhysicalMapping_;
    }

    /// Returns a pointer to the physicalGeometry object.

    const PhysicalGeometry * const
    ElementGeometry::getPhysicalGeometry() const
    {
        return physicalGeometry_;
    }

    /// Returns a pointer to the referenceGeometry object.

    const ReferenceGeometry * const
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

    PointPhysical
    ElementGeometry::referenceToPhysical(const PointReferenceT& pointReference) const
    {
        return referenceToPhysicalMapping_->transform(pointReference);
    }

    /// This method gets a PointReference and returns the corresponding jacobian of the
    /// referenceToPhysicalMapping.

    Jacobian
    ElementGeometry::calcJacobian(const PointReferenceT& pointReference) const
    {
        return referenceToPhysicalMapping_->calcJacobian(pointReference);
    }

    std::size_t
    ElementGeometry::getNrOfNodes() const
    {
        return physicalGeometry_->getNumberOfNodes();
    }

}

#endif
