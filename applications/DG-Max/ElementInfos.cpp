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

#include "ElementInfos.h"
#include "Geometry/PointReference.h"
#include "Geometry/ReferenceGeometry.h"
#include "Geometry/PointPhysical.h"
#include "Base/ElementCacheData.h"
#include "Base/FaceCacheData.h"


double MaxwellData::StabCoeff_;

ElementInfos::ElementInfos(const Base::Element& element)
{
    const Geometry::PointReference<DIM>& p = element.getReferenceGeometry()->getCenter();
    Geometry::PointPhysical<DIM> pPhys = element.referenceToPhysical(p);
    
    /*
    if((pPhys[1]-0.5)*(pPhys[1]-0.5)+(pPhys[2]-0.5)*(pPhys[2]-0.5)<.25*.25){
    if((pPhys[0]<0.3)||pPhys[0]>0.7||pPhys[1]<0.3||pPhys[1]>0.7){
    epsilon_ = 1;
    }else{
    epsilon_=1;
    }
     */
    
    epsilon_=1;
    Jacobian_ = element.calcJacobian(p);
    determinant_ = Jacobian_.determinant();
}



MaxwellData::MaxwellData(int numberOfIntervals, int polynomialOrder)
        : Sigma_(0), StartTime_(0), EndTime_(0.1), NumberOfIntervals_(numberOfIntervals), PolynomialOrder_(polynomialOrder)
{
    StabCoeff_ = numberOfIntervals * (polynomialOrder + 1) * (polynomialOrder + 3);
}
