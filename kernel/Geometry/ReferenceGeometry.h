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

#ifndef REFERENCEGEOMETRY_HH
#define REFERENCEGEOMETRY_HH

#include "Geometry/Mappings/MappingCodimensions.h"
#include "Geometry/Mappings/RefinementMapping.h"
#include "PointReference.h"

#include <map>
#include <unordered_map>
#include <iostream>
#include <vector>

namespace LinearAlgebra
{
    class MiddleSizeVector;
}

namespace QuadratureRules
{
    // forward declaration
    class GaussQuadratureRule;
}

namespace Base
{
    class BaseBasisFunction;
}

namespace Geometry
{
    
    enum class ReferenceGeometryType
    {
        POINT, LINE, TRIANGLE, SQUARE, TETRAHEDRON, PYRAMID, CUBE, TRIANGULARPRISM, HYPERCUBE
    };

    /*! \class ReferenceGeometry
     * \brief ReferenceGeometry stores a the information of a unitary geometry where the integration is made.
     * \details
     * ReferenceGeometry stores the information of the corresponding reference object, which
     * are used for integration routines. It is a pure virtual class; an interface for every
     * particular reference geometry. The information is a container of the actual reference
     * points.
     *
     * Constructors are protected, to avoid the creation of two identical physical geometries.
     * Only one of every needed type is necessary.
     */
    class ReferenceGeometry : public RefinementMapping, public MappingCodimensions
    {
    public:
        using String = std::string;
        using ListOfIndexesT = std::vector<std::size_t>;

    public:
        
        virtual ~ReferenceGeometry() = default;
        
        ReferenceGeometry(const ReferenceGeometry& other) = delete;

        /// \brief Check whether a given point is within the ReferenceGeometry.
        virtual bool isInternalPoint(const PointReference<0>& point) const
        {
            logger(ERROR, "You passed a point of the wrong dimension");
            return false;
        }

        /// \brief Check whether a given point is within the ReferenceGeometry.
        virtual bool isInternalPoint(const PointReference<1>& point) const
        {
            logger(ERROR, "You passed a point of the wrong dimension");
            return false;
        }

        /// \brief Check whether a given point is within the ReferenceGeometry.
        virtual bool isInternalPoint(const PointReference<2>& point) const
        {
            logger(ERROR, "You passed a point of the wrong dimension");
            return false;
        }

        /// \brief Check whether a given point is within the ReferenceGeometry.
        virtual bool isInternalPoint(const PointReference<3>& point) const
        {
            logger(ERROR, "You passed a point of the wrong dimension");
            return false;
        }

        /// \brief Check whether a given point is within the ReferenceGeometry.
        virtual bool isInternalPoint(const PointReference<4>& point) const
        {
            logger(ERROR, "You passed a point of the wrong dimension");
            return false;
        }

        /// \brief Each reference geometry knows its center of mass.
        virtual const PointReferenceBase& getCenter() const = 0;

        /// \brief Return number of nodes of the reference shape.
        virtual std::size_t getNumberOfNodes() const = 0;

        ReferenceGeometryType getGeometryType() const
        {
            return geometryType_;
        }
        
        /// \brief Given a local index, return (assign to point) the corresponding node.
        virtual const PointReferenceBase& getReferenceNodeCoordinate(const std::size_t& localIndex) const = 0;

        std::size_t getLocalNodeIndexFromFaceAndIndexOnFace(std::size_t face, std::size_t node) const
        {
            logger.assert(face < getNumberOfCodim1Entities(), "Asked for face %, but there are only % faces", face, getNumberOfCodim1Entities());
            logger.assert(node < getCodim1ReferenceGeometry(face)->getNumberOfNodes(), "Asked for node %, but there are only % nodes", node, getCodim1ReferenceGeometry(face)->getNumberOfNodes());
            return getCodim1EntityLocalIndices(face)[node];
        }

        /// \brief For debugging and checkpointing: a human-readable name.
        std::string getName() const
        {
            return name;
        }

        // ================================== Quadrature rules =====================================
        
        /// \brief Get a valid quadrature for this geometry.
        const QuadratureRules::GaussQuadratureRule* getGaussQuadratureRule(std::size_t order) const;
        
    protected:
        ReferenceGeometry(std::size_t numberOfNodes, std::size_t DIM, const ReferenceGeometryType& geoT, std::initializer_list<double> center);

        /// An identifier of the type of referenceGeometry, that some say shouldn't be used.
        const ReferenceGeometryType geometryType_;

        std::string name;
        
    };

}

#endif
