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

#ifndef BaseBasisFunction_h
#define BaseBasisFunction_h

#include "LinearAlgebra/SmallVector.h"

namespace Geometry
{
    template<std::size_t DIM>
    class PointReference;
}

namespace Base
{
    
    class BaseBasisFunction
    {
    public:
        BaseBasisFunction() = default;
        BaseBasisFunction(const BaseBasisFunction &other) = default;
        BaseBasisFunction& operator=(const BaseBasisFunction &other) = default;

        virtual ~ BaseBasisFunction()
        {
        }

        //we have to manually specify reasonable choices for template parameters in virtual functions
        virtual double eval(const Geometry::PointReference<1>& p) const
        {
            logger(ERROR, "The reference point you passed has the wrong dimension");
            return 0;
        }

        virtual double eval(const Geometry::PointReference<2>& p) const
        {
            logger(ERROR, "The reference point you passed has the wrong dimension");
            return 0;
        }

        virtual double eval(const Geometry::PointReference<3>& p) const
        {
            logger(ERROR, "The reference point you passed has the wrong dimension");
            return 0;
        }

        virtual double eval(const Geometry::PointReference<4>& p) const
        {
            logger(ERROR, "The reference point you passed has the wrong dimension");
            return 0;
        }

        virtual void eval(const Geometry::PointReference<2>& p, LinearAlgebra::SmallVector<2>& ret) const
        {
            logger(ERROR, "The reference point you passed has the wrong dimension");
        }

        virtual void eval(const Geometry::PointReference<3>& p, LinearAlgebra::SmallVector<3>& ret) const
        {
            logger(ERROR, "The reference point you passed has the wrong dimension");
        }

        virtual void eval(const Geometry::PointReference<4>& p, LinearAlgebra::SmallVector<4>& ret) const
        {
            logger(ERROR, "The reference point you passed has the wrong dimension");
        }

        virtual double evalDeriv0(const Geometry::PointReference<1>& p) const
        {
            logger(ERROR, "The reference point you passed has the wrong dimension");
            return 0;
        }

        virtual double evalDeriv0(const Geometry::PointReference<2>& p) const
        {
            logger(ERROR, "The reference point you passed has the wrong dimension");
            return 0;
        }

        virtual double evalDeriv0(const Geometry::PointReference<3>& p) const
        {
            logger(ERROR, "The reference point you passed has the wrong dimension");
            return 0;
        }

        virtual double evalDeriv0(const Geometry::PointReference<4>& p) const
        {
            logger(ERROR, "The reference point you passed has the wrong dimension");
            return 0;
        }

        virtual double evalDeriv1(const Geometry::PointReference<2>& p) const
        {
            logger(ERROR, "The reference point you passed has the wrong dimension");
            return 0;
        }

        virtual double evalDeriv1(const Geometry::PointReference<3>& p) const
        {
            logger(ERROR, "The reference point you passed has the wrong dimension");
            return 0;
        }

        virtual double evalDeriv1(const Geometry::PointReference<4>& p) const
        {
            logger(ERROR, "The reference point you passed has the wrong dimension");
            return 0;
        }

        virtual double evalDeriv2(const Geometry::PointReference<3>& p) const
        {
            logger(ERROR, "The reference point you passed has the wrong dimension");
            return 0;
        }

        virtual double evalDeriv2(const Geometry::PointReference<4>& p) const
        {
            logger(ERROR, "The reference point you passed has the wrong dimension");
            return 0;
        }

        virtual double evalDeriv3(const Geometry::PointReference<4>& p) const
        {
            logger(ERROR, "The reference point you passed has the wrong dimension");
            return 0;
        }

        virtual LinearAlgebra::SmallVector<2> evalCurl(const Geometry::PointReference<2>& p) const
        {
            logger(ERROR, "The curl of a scalar valued basis function is not implemented. Perhaps you meant evalDeriv?");
            return LinearAlgebra::SmallVector<2>();
        }
        
        virtual LinearAlgebra::SmallVector<3> evalCurl(const Geometry::PointReference<3>& p) const
        {
            logger(ERROR, "The curl of a scalar valued basis function is not implemented. Perhaps you meant evalDeriv?");
            return LinearAlgebra::SmallVector<3>();
        }

        virtual double evalDiv(const Geometry::PointReference<2>& p) const
        {
            logger(ERROR, "The divergence of a scalar valued basis function is not implemented. Perhaps you meant evalDeriv?");
            return 0;
        }

        virtual double evalDiv(const Geometry::PointReference<3>& p) const
        {
            logger(ERROR, "The divergence of a scalar valued basis function is not implemented. Perhaps you meant evalDeriv?");
            return 0;
        }

        virtual double evalDiv(const Geometry::PointReference<4>& p) const
        {
            logger(ERROR, "The divergence of a scalar valued basis function is not implemented. Perhaps you meant evalDeriv?");
            return 0;
        }

        virtual LinearAlgebra::SmallVector<1> evalDeriv(const Geometry::PointReference<1>& p) const;
        virtual LinearAlgebra::SmallVector<2> evalDeriv(const Geometry::PointReference<2>& p) const;
        virtual LinearAlgebra::SmallVector<3> evalDeriv(const Geometry::PointReference<3>& p) const;
        virtual LinearAlgebra::SmallVector<4> evalDeriv(const Geometry::PointReference<4>& p) const;
    };

}

#endif // BaseBasisFunction_h
