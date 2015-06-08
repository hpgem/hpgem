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

#ifndef ELEMENTGEOMETRY_H_
#define ELEMENTGEOMETRY_H_

#include <vector>
#include <iostream>
#include "Point.h"
#include "PointPhysical.h"
#include "Jacobian.h"
#include "Mappings/MappingReferenceToPhysical.h"

#include "ReferenceTetrahedron.h"
#include "ReferenceLine.h"
#include "ReferenceSquare.h"
#include "ReferenceTriangle.h"
#include "ReferencePyramid.h"
#include "ReferenceTriangularPrism.h"
#include "ReferenceCube.h"
#include "ReferenceHypercube.h"

#include "PhysicalTetrahedron.h"
#include "PhysicalLine.h"
#include "PhysicalQuadrilateral.h"
#include "PhysicalTriangle.h"
#include "PhysicalPyramid.h"
#include "PhysicalTriangularPrism.h"
#include "PhysicalHexahedron.h"
#include "PhysicalOctachoron.h"

#include "Mappings/MappingReferenceToPhysical.h"
#include "Mappings/MappingToPhysHypercubeLinear.h"
#include "Mappings/MappingToPhysSimplexLinear.h"
#include "Mappings/MappingToPhysPyramid.h"
#include "Mappings/MappingToPhysTriangularPrism.h"

#include "RefinementLine.h"
#include "RefinementTriangle.h"
#include "RefinementQuadrilateral.h"
#include "RefinementTetrahedron.h"
#include "RefinementPyramid.h"
#include "RefinementTriangularPrism.h"
#include "RefinementHexahedron.h"
#include "RefinementHypercube.h"

#include "PointReference.h"

namespace Geometry
{
    template<std::size_t DIM>
    class PointReference;
    template<std::size_t DIM>
    class PointPhysical;
    class MappingReferenceToPhysical;
    class PhysicalGeometryBase;
    class ReferenceGeometry;
    class RefinementGeometry;
    template<std::size_t dimFrom, std::size_t dimTo>
    class Jacobian;
    
    class ElementGeometry
    {
    public:
        using PointIndexT = std::size_t;
        using VectorOfPointIndexesT = std::vector<PointIndexT>;
    public:
        
        /// New style constructor with one less pass
        template<std::size_t DIM>
        ElementGeometry(const VectorOfPointIndexesT& globalNodeIndexes, std::vector<PointPhysical<DIM> >& nodes);

        /// Copy constructor
        ElementGeometry(const ElementGeometry& other);

        virtual ~ElementGeometry();

        /// Returns a pointer to the referenceToPhysicalMapping
        const MappingReferenceToPhysical* getReferenceToPhysicalMap() const;
        MappingReferenceToPhysical* getReferenceToPhysicalMap();

        /// Returns a pointer to the physicalGeometry object.
        const PhysicalGeometryBase* getPhysicalGeometry() const;
        /// Returns a pointer to the physicalGeometry object.
        PhysicalGeometryBase* getPhysicalGeometry();
        /// Returns a pointer to the physicalGeometry object.
        std::size_t getNrOfNodes() const;
        /// Returns a pointer to the referenceGeometry object.
        const ReferenceGeometry* getReferenceGeometry() const;
        ReferenceGeometry* getReferenceGeometry();
        /// Returns a pointer to the refinementGeometry object.
        const RefinementGeometry* getRefinementGeometry() const;
        /// This method gets a PointReference, which specifies a coordinate in the ReferenceGeometry,
        /// and returns a PointPhysical which is the corresponding point in the PhysicalGeometry,
        /// given the mapping.
        template<std::size_t DIM>
        PointPhysical<DIM> referenceToPhysical(const PointReference<DIM>& pointReference) const;

        /// This method gets a PointReference and returns the corresponding jacobian of the
        /// referenceToPhysicalMapping.
        template<std::size_t DIM>
        Jacobian<DIM, DIM> calcJacobian(const PointReference<DIM>& pointReference) const;

        void enableRefinement();

    public:
        /// Output operator.
        friend std::ostream& operator <<(std::ostream& os, const ElementGeometry& elementGeometry);
        
    private:

        template<std::size_t DIM>
        static ReferenceGeometry* createReferenceGeometry(std::size_t size);

        template<std::size_t DIM>
        static PhysicalGeometry<DIM>* createPhysicalGeometry(const VectorOfPointIndexesT& globalNodeIndexes, std::vector<PointPhysical<DIM> >& nodes, const ReferenceGeometry* const geo);

        template<std::size_t DIM>
        static MappingReferenceToPhysical* createMappings(std::size_t size, const PhysicalGeometry<DIM>* const pGeo);

    protected:
        /// The corresponding referenceGeometry object, for integration.
        ReferenceGeometry* const referenceGeometry_;

        /// The physicalGeometry object contains pointers to the actual physical points, and
        /// a container of global node indexes.
        PhysicalGeometryBase* physicalGeometry_;

        /// The referenceToPhysicalMapping relates the coordinates of the reference object to the
        /// physical object; basically a matrix transformation.
        MappingReferenceToPhysical* referenceToPhysicalMapping_;

        /// The corresponding refinementGeometry object
        RefinementGeometry* refinementGeometry_;
    };



    /// This method gets a PointReference, which specifies a coordinate in the ReferenceGeometry,
    /// and returns a PointPhysical which is the corresponding point in the PhysicalGeometry,
    /// given the mapping.

    template<std::size_t DIM>
    PointPhysical<DIM> ElementGeometry::referenceToPhysical(const PointReference<DIM>& pointReference) const
    {
        return referenceToPhysicalMapping_->transform(pointReference);
    }

    /// This method gets a PointReference and returns the corresponding jacobian of the
    /// referenceToPhysicalMapping.

    template<std::size_t DIM>
    Jacobian<DIM, DIM> ElementGeometry::calcJacobian(const PointReference<DIM>& pointReference) const
    {
        return referenceToPhysicalMapping_->calcJacobian(pointReference);
    }

    template<std::size_t DIM>
    ReferenceGeometry *
    ElementGeometry::createReferenceGeometry(std::size_t size)
    {
        switch (size)
        {
            //select a proper type based on the number of nodes a reference geometry should have
            case 2:
                logger.assert(DIM==1, "This Dimension does not contain entities with 2 nodes");
                logger(VERBOSE, "ElementGeometry created a reference line.");
                return &ReferenceLine::Instance();
            case 3:
                logger.assert(DIM==2, "This Dimension does not contain entities with 3 nodes");
                logger(VERBOSE, "ElementGeometry created a reference triangle.");
                return &ReferenceTriangle::Instance();
            case 4:
                if (DIM == 2)
                {
                    logger(VERBOSE, "ElementGeometry created a reference square.");
                    return &ReferenceSquare::Instance();
                }
                else if (DIM == 3)
                {
                    logger(VERBOSE, "ElementGeometry created a reference tetrahedron.");
                    return &ReferenceTetrahedron::Instance();
                }
                else
                {
                    logger(ERROR, "This dimension does not contain entities with 4 nodes. \n");
                }
                break;
            case 5:
                logger.assert(DIM==3, "This Dimension does not contain entities with 5 nodes");
                logger(VERBOSE, "ElementGeometry created a reference pyramid.");
                return &ReferencePyramid::Instance();
            case 6:
                logger.assert(DIM==3, "This Dimension does not contain entities with 6 nodes");
                logger(VERBOSE, "ElementGeometry created a reference triangular prism.");
                return &ReferenceTriangularPrism::Instance();
            case 8:
                logger.assert(DIM==3, "This Dimension does not contain entities with 8 nodes");
                logger(VERBOSE, "ElementGeometry created a reference cube.");
                return &ReferenceCube::Instance();
            case 16:
                logger.assert(DIM==4, "This Dimension does not contain entities with 16 nodes");
                logger(VERBOSE, "ElementGeometry created a reference hypercube.");
                return &ReferenceHypercube::Instance();
            default:
                logger(FATAL, "No know entities contain this many nodes. \n");
        }
        return 0;
    }

    template<std::size_t DIM>
    PhysicalGeometry<DIM> *
    ElementGeometry::createPhysicalGeometry(const VectorOfPointIndexesT& globalNodeIndexes, std::vector<PointPhysical<DIM> >& nodes, const ReferenceGeometry * const geo)
    {
        logger.assert(geo!=nullptr, "Invalid reference geometry passed");
        return new PhysicalGeometry<DIM>(globalNodeIndexes, nodes, geo);
    }

    template<std::size_t DIM>
    MappingReferenceToPhysical *
    ElementGeometry::createMappings(std::size_t size, const PhysicalGeometry<DIM> * const pGeo)
    {
        logger(ERROR, "DIM may range from 1 to 4, but it seems to be %", DIM);
    }

    template<>
    inline MappingReferenceToPhysical *
    ElementGeometry::createMappings(std::size_t size, const PhysicalGeometry<1> * const pGeo)
    {
        logger.assert(pGeo!=nullptr, "Invalid physical geometry passed");
        logger.assert(size == 2, "1D can only map to a line");
        logger(VERBOSE, "ElementGeometry created a mapping for a line.");
        return new Geometry::MappingToPhysHypercubeLinear<1>(pGeo);
    }

    template<>
    inline MappingReferenceToPhysical *
    ElementGeometry::createMappings(std::size_t size, const PhysicalGeometry<2> * const pGeo)
    {
        logger.assert(pGeo!=nullptr, "Invalid physical geometry passed");
        switch (size)
        {
            case 3:
                logger(VERBOSE, "ElementGeometry created a mapping for a triangle.");
                return new Geometry::MappingToPhysSimplexLinear<2>(pGeo);
            case 4:
                logger(VERBOSE, "ElementGeometry created a mapping for a square.");
                return new Geometry::MappingToPhysHypercubeLinear<2>(pGeo);
        }
        logger(FATAL, "No know entities contain this many nodes. \n");
        return nullptr;
    }

    template<>
    inline MappingReferenceToPhysical *
    ElementGeometry::createMappings(std::size_t size, const PhysicalGeometry<3> * const pGeo)
    {
        logger.assert(pGeo!=nullptr, "Invalid physical geometry passed");
        switch (size)
        {
            case 4:
                logger(VERBOSE, "ElementGeometry created a mapping for a tetrahedron.");
                return new Geometry::MappingToPhysSimplexLinear<3>(pGeo);
            case 5:
                logger(VERBOSE, "ElementGeometry created a mapping for a pyramid.");
                return new Geometry::MappingToPhysPyramid(pGeo);
            case 6:
                logger(VERBOSE, "ElementGeometry created a mapping for a triangular prism.");
                return new Geometry::MappingToPhysTriangularPrism(pGeo);
            case 8:
                logger(VERBOSE, "ElementGeometry created a mapping for a cube.");
                return new Geometry::MappingToPhysHypercubeLinear<3>(pGeo);
        }
        logger(FATAL, "No know entities contain this many nodes. \n");
        return nullptr;
    }

    template<>
    inline MappingReferenceToPhysical *
    ElementGeometry::createMappings(std::size_t size, const PhysicalGeometry<4> * const pGeo)
    {
        logger.assert(pGeo!=nullptr, "Invalid physical geometry passed");
        logger.assert(size == 16, "4D can only map to a hypercube");
        logger(VERBOSE, "ElementGeometry created a mapping for a hypercube.");
        return new Geometry::MappingToPhysHypercubeLinear<4>(pGeo);
    }

    template<std::size_t DIM>
    ElementGeometry::ElementGeometry(const VectorOfPointIndexesT& globalNodeIndexes, std::vector<PointPhysical<DIM> >& nodes)
            : referenceGeometry_(ElementGeometry::createReferenceGeometry<DIM>(globalNodeIndexes.size())),
        physicalGeometry_(ElementGeometry::createPhysicalGeometry(globalNodeIndexes, nodes, referenceGeometry_)),
        referenceToPhysicalMapping_(ElementGeometry::createMappings<DIM>(globalNodeIndexes.size(), static_cast<PhysicalGeometry<DIM>*>(physicalGeometry_))),
        refinementGeometry_(nullptr) //refinement is turned off by default, to  enable it one needs to call enableRefinement
    {
    }
}
#endif
