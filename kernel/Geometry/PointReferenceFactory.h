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

#ifndef POINTREFERENCEFACTORY_H_
#define POINTREFERENCEFACTORY_H_

#include "Geometry/PointReference.h"
#include "Base/BaseBasisFunction.h"

namespace Geometry
{
    class PointReference;
    
    ///\brief make unique reference coordinates
    ///\details we want to evaluate basis functions only once for each reference point
    /// but it is quite difficult in places where they are needed to check whether or not a reference point has already been created
    /// It is assumed reference points will not be created very often in a simulation, so this factory is allowed to be relatively slow
    /// returns pointer to const PointReference, because the PointReference might be used in a ton of other places (and I cant imagine why anyone would want to move them about)
    /// There is a tolerance for when points are considered the same in makePoint(Point)
    class PointReferenceFactory
    {
    public:
        static PointReferenceFactory* instance()
        {
            static PointReferenceFactory theInstance;
            return &theInstance;
        }

        void removeBasisFunctionData(const Base::BaseBasisFunction* function);

        const PointReference* makePoint(const Point& p);

        const PointReference* makePoint(std::size_t DIM)
        {
            return makePoint(Point(DIM));
        }


        const PointReference* makePoint(std::initializer_list<double> data)
        {
            return makePoint(Point(data));
        }

        const PointReference* makePoint(double coords[], std::size_t DIM)
        {
            return makePoint(Point(coords, DIM));
        }

        const PointReference* makePoint(const LinearAlgebra::NumericalVector& coord)
        {
            return makePoint(Point(coord));
        }

    private:
        PointReferenceFactory();
        ~PointReferenceFactory();
        PointReferenceFactory(const PointReferenceFactory&) = delete;
        PointReferenceFactory(PointReferenceFactory&&) = delete;

        std::vector<PointReference*> points_;
    };

} /* namespace Geometry */

#endif /* POINTREFERENCEFACTORY_H_ */
