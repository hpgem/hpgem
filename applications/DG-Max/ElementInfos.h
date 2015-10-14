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

//this file has a container data structure for everything you want to know on a per element basis
#ifndef Elementinfos_h
#define Elementinfos_h

#include "Base/UserData.h"
#include "Base/Element.h"
#include "Base/Face.h"
#include "BasisFunctionCollection_Curl.h"
#include "Base/GlobalData.h"
#include "Geometry/Jacobian.h"


//one there is a default way in hpGEM to configure code without haveing to recompile this and some other things should be grouped in another file
/**
 * Stores some parameters that should be available everywhere.
 * Attributes need not be constant since this struct will be stored in a constant field.
 */
struct MaxwellData : public Base::GlobalData
{
    double Sigma_; //allows for lossy media - untested except maybe by Domokos
    static double StabCoeff_; //a_F for IP or 1+eta_F/4 for BR
    int NumberOfIntervals_;
    int PolynomialOrder_;
    double StartTime_;
    double EndTime_;

    MaxwellData(int numberOfIntervals, int polynomialOrder);
};



/**
 * store some usefull information that needs to be computed everytime at the beginning of an integrand
 * specialized for tetrahedra
 */


class ElementInfos : public UserElementData
{
public:
    double determinant_;
    Geometry::Jacobian<DIM, DIM> Jacobian_, inverse_;
    double epsilon_;

    ElementInfos(const Base::Element& element);
};

#endif  