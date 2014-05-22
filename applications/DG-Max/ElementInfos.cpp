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

#include "ElementInfos.hpp"
#include "Geometry/PointReference.hpp"
#include "Geometry/ReferenceGeometry.hpp"
#include "Geometry/PointPhysical.hpp"
#include "Base/ElementCacheData.hpp"
#include "Base/FaceCacheData.hpp"

void InvertAndTranspose(Geometry::Jacobian& orig, Geometry::Jacobian& inverse){
    //direct computation using the definitions of the inverse and the transpose
    inverse(0,0)=orig(1,1)*orig(2,2)-orig(1,2)*orig(2,1);
    inverse(1,0)=orig(2,1)*orig(0,2)-orig(2,2)*orig(0,1);
    inverse(2,0)=orig(0,1)*orig(1,2)-orig(0,2)*orig(1,1);
    inverse(0,1)=orig(1,2)*orig(2,0)-orig(1,0)*orig(2,2);
    inverse(1,1)=orig(2,2)*orig(0,0)-orig(2,0)*orig(0,2);
    inverse(2,1)=orig(0,2)*orig(1,0)-orig(0,0)*orig(1,2);
    inverse(0,2)=orig(1,0)*orig(2,1)-orig(1,1)*orig(2,0);
    inverse(1,2)=orig(2,0)*orig(0,1)-orig(2,1)*orig(0,0);
    inverse(2,2)=orig(0,0)*orig(1,1)-orig(0,1)*orig(1,0);
    inverse/=orig.determinant();
}

void FunctionCache::getFunctionValuesVector(const Base::Element* element, const PointElementReferenceT& point, std::vector< LinearAlgebra::NumericalVector >& values){
    values=valueCache_[point];
    if(values.empty()){
	LinearAlgebra::NumericalVector value(3);
	for(int j=0;j<element->getNrOfBasisFunctions();++j){
	    element->basisFunction(j,point,value);
	    values.push_back(value);
	} 
	valueCache_[point]=values;
    }
}

void FunctionCache::getFunctionCurlsVector(const Base::Element* element, const PointElementReferenceT& point, std::vector< LinearAlgebra::NumericalVector >& curls){
    curls=curlCache_[point];
    if(curls.empty()){
	LinearAlgebra::NumericalVector curl(3);
	for(int j=0;j<element->getNrOfBasisFunctions();++j){
	    element->basisFunctionCurl(j,point,curl);
	    curls.push_back(curl);
	}
	curlCache_[point]=curls;
    }
}

myMap FunctionCache::valueCache_;
myMap FunctionCache::curlCache_;

ElementInfos::ElementInfos(const Base::Element& element):inverse_(3,3),Jacobian_(3,3){
    PointElementReferenceT p(3);
    Geometry::PointPhysical pPhys(3);
    element.getReferenceGeometry()->getCenter(p);
    element.referenceToPhysical(p,pPhys);
    //not quite sure about the best way to implement this; this works for now
    if((pPhys[1]-0.5)*(pPhys[1]-0.5)+(pPhys[2]-0.5)*(pPhys[2]-0.5)<.25*.25){
    ///\bug Does not check that there are element boundaries at any of the discontinuities
    //if((pPhys[0]<0.3)||pPhys[0]>0.7||pPhys[1]<0.3||pPhys[1]>0.7){
	epsilon_=1;
    }else{
	epsilon_=1;
    }
    //the jacobian of a tetrahedron is constant.
    element.calcJacobian(p,Jacobian_);
    determinant_=Jacobian_.determinant();
    InvertAndTranspose(Jacobian_,inverse_);
    //cout<<"element "<<element.getID()<<" has jacobean determinant "<<determinant_<<endl;
}

void ElementInfos::makeFunctionValuesVector(const Base::Element* element, const PointElementReferenceT& point, std::vector< LinearAlgebra::NumericalVector >& values){
    FunctionCache::getFunctionValuesVector(element,point,values);
    for(int j=0;j<element->getNrOfBasisFunctions();++j){
	//3D coordinate transformations -- something for a GPU?
	values[j]=inverse_*values[j];
    }        
}

void ElementInfos::makeFunctionCurlsVector(const Base::Element* element, const PointElementReferenceT& point, std::vector< LinearAlgebra::NumericalVector >& curls){
    FunctionCache::getFunctionCurlsVector(element,point,curls);
    for(int j=0;j<element->getNrOfBasisFunctions();++j){
	//3D coordinate transformations -- something for a GPU?
	curls[j]=Jacobian_*curls[j]/determinant_;
    } 
}

MaxwellData::MaxwellData(int numberOfIntervals, int polynomialOrder):Sigma_(0),
    StabCoeff_(3*numberOfIntervals*polynomialOrder*(polynomialOrder+2)+1),
    StartTime_(0),
    EndTime_(0.1),
    NumberOfIntervals_(numberOfIntervals),
    PolynomialOrder_(polynomialOrder){}
