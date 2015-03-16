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
#include "ReferenceGeometry.h"
#include "Integration/QuadratureRules/AllGaussQuadratureRules.h"
#include "Geometry/PointReference.h"
#include "Base/BaseBasisFunction.h"

#include <stdexcept>

namespace Geometry
{
    
    ReferenceGeometry::ReferenceGeometry(const ReferenceGeometryType& geoT)
            : geometryType_(geoT)
    {
        
    }
    
    ReferenceGeometry::ReferenceGeometry(std::size_t numberOfNodes, std::size_t DIM, const ReferenceGeometryType& geoT)
            : points_(numberOfNodes, DIM), geometryType_(geoT)
    {
        
    }
    
    const QuadratureRules::GaussQuadratureRule* const ReferenceGeometry::getGaussQuadratureRule(std::size_t order) const
    {
        return QuadratureRules::AllGaussQuadratureRules::instance().getRule(this, order);
    }
    
    ReferenceGeometry::ReferenceGeometry(const ReferenceGeometry& other)
            : points_(other.points_), geometryType_(other.geometryType_)
    {
        
    }
    
    double ReferenceGeometry::getBasisFunctionValue(const Base::BaseBasisFunction* function, const PointReference& p)
    {
        logger.assert(function!=nullptr, "Invalid basis function passed");
        try
        {
            return basisfunctionValues_[function].at(p);
        }
        catch (std::out_of_range&)
        {
            basisfunctionValues_[function][p] = function->eval(p);
            return basisfunctionValues_[function].at(p);
        }
        return function->eval(p);
    }
    
    LinearAlgebra::NumericalVector&
    ReferenceGeometry::getBasisFunctionDerivative(const Base::BaseBasisFunction* function, const PointReference& p)
    {
        logger.assert(function!=nullptr, "Invalid basis function passed");
        try
        {
            return basisfunctionDerivatives_[function].at(p);
        }
        catch (std::out_of_range&)
        {
            basisfunctionDerivatives_[function][p] = function->evalDeriv(p);
            return basisfunctionDerivatives_[function].at(p);
        }
    }
    
    const PointReference& ReferenceGeometry::getNode(const std::size_t& localIndex) const
    {
        logger.assert(localIndex<getNumberOfNodes(), "Asked for node %, but there are only % nodes", localIndex, getNumberOfNodes());
        return points_[localIndex];
    }

}

std::size_t Geometry::PointHasher::operator()(const Geometry::PointReference& point) const
{
    std::hash<double> hasher;
    std::size_t ret = 0;
    for (std::size_t i = 0; i < point.size(); ++i)
    {
        ret ^= hasher(point[i]) + 0x9e3779b9 + (ret << 6) + (ret >> 2);
    }
    return ret;
}
