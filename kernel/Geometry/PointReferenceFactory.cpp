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

#include <PointReferenceFactory.h>
#include <limits>
#include <algorithm>

namespace Geometry
{
    
    PointReferenceFactory::PointReferenceFactory()
    {
        points_.push_back(new PointReference{std::numeric_limits<double>::quiet_NaN()});
    }

    //developer note: during normal computation this routine is called a few times when the mappings and the quadrature rules are set up
    //during the unit tests this routine is called many times, with points different from the quadrature rules, because most unit tests
    //should function independent of the quadrature rules and they should try to make sure everything works fine even for reference points
    //that are not listed in the quardature rules.This means the quadrature rules create many reference points and use them only once.
    //this is likely to be slow with the current set up
    const PointReference* PointReferenceFactory::makePoint(const Point& p)
    {
        //cannot find nan back using bounded difference
        if(p.size() == 1 && std::isnan(p[0]))
        {
            return points_[0];
        }
        auto existingPoint = std::find_if(points_.begin(), points_.end(), [&p](const PointReference* other){
            //slow version: return p.size() == other->size() && std::sqrt((p - (*other)) * (p - (*other))) < 1e-12 (can be used for stack based points for less size checking and less moving up and down the stack)
            if(p.size() != other->size())
            {
                return false;
            }
            double squaredSize = 0.;
            for(std::size_t i = 0; i < p.size(); ++i)
            {
                squaredSize += (p[i] - (*other)[i]) * (p[i] - (*other)[i]);
            }
            return std::sqrt(squaredSize) < 1e-12;
        });
        if(existingPoint!=points_.end())
        {
            return *existingPoint;
        }
        else
        {
            points_.push_back(new PointReference(p));
            return points_.back();
        }
    }
    
    PointReferenceFactory::~PointReferenceFactory()
    {
        for(const PointReference* point : points_)
        {
            delete point;
        }
    }

    void PointReferenceFactory::removeBasisFunctionData(const Base::BaseBasisFunction* function)
    {
        for(PointReference* point : points_)
        {
            point->removeBasisFunctionData(function);
        }
    }

} /* namespace Geometry */
