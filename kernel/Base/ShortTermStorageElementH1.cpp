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


#include "Base/ShortTermStorageElementH1.hpp"
#include "ElementCacheData.hpp"

void Base::ShortTermStorageElementH1::computeData() {
    ShortTermStorageElementBase::computeData();
    basisFunctionValues_.resize(element_->getNrOfBasisFunctions());
    basisFunctionDerivatives_.resize(element_->getNrOfBasisFunctions());
    basisFunctionIndividualDerivatives_.resize(element_->getNrOfBasisFunctions());
    for(int i=0;i<element_->getNrOfBasisFunctions();++i){
        basisFunctionValues_[i].resize(1);
        basisFunctionValues_[i][0]=element_->basisFunction(i,currentPoint_);
        basisFunctionDerivatives_[i].resize(currentPoint_.size());
        element_->basisFunctionDeriv(i,currentPoint_,basisFunctionDerivatives_[i],this);
        basisFunctionIndividualDerivatives_[i].resize(currentPoint_.size());
        for(int j=0;j<currentPoint_.size();++j){
            basisFunctionIndividualDerivatives_[i][j]=element_->basisFunctionDeriv(i,j,currentPoint_);
        }
    }
}


double Base::ShortTermStorageElementH1::basisFunction(unsigned int i, const PointReferenceT& p) {
    if(!(p==currentPoint_)){
        currentPoint_=p;
        computeData();
    }
    return basisFunctionValues_[i][0];
}


double Base::ShortTermStorageElementH1::basisFunction(unsigned int i, const PointReferenceT& p) const {
	if(!(p==currentPoint_)){
		std::cout<<"WARNING: you are using a slow operator";
		return element_->basisFunction(i,p);
	}
	return basisFunctionValues_[i][0];
}

void Base::ShortTermStorageElementH1::basisFunction(unsigned int i, const PointReferenceT& p, LinearAlgebra::NumericalVector& ret) {
	if(!(p==currentPoint_)){
		currentPoint_=p;
		computeData();
	}
	ret=basisFunctionValues_[i];
}

void Base::ShortTermStorageElementH1::basisFunction(unsigned int i, const PointReferenceT& p, LinearAlgebra::NumericalVector& ret) const {
	ret=basisFunctionValues_[i];
	if(!(p==currentPoint_)){
		std::cout<<"WARNING: you are using a slow operator";
		element_->basisFunction(i,p,ret);
	}
}


void Base::ShortTermStorageElementH1::basisFunctionDeriv(unsigned int i, const PointReferenceT& p, LinearAlgebra::NumericalVector& ret,const Element* ) {
	if(!(p==currentPoint_)){
		currentPoint_=p;
		computeData();
	}
	ret=basisFunctionDerivatives_[i];
}

void Base::ShortTermStorageElementH1::basisFunctionDeriv(unsigned int i, const PointReferenceT& p, LinearAlgebra::NumericalVector& ret,const Element*) const {
	ret=basisFunctionDerivatives_[i];
	if(!(p==currentPoint_)){
		std::cout<<"WARNING: you are using a slow operator";
		element_->basisFunctionDeriv(i,p,ret,this);
	}
}

double Base::ShortTermStorageElementH1::basisFunctionDeriv(unsigned int i, unsigned int jDir, const PointReferenceT& p) {
        if(!(p==currentPoint_)){
                currentPoint_=p;
                computeData();
        }
        return basisFunctionIndividualDerivatives_[i][jDir];
}

double Base::ShortTermStorageElementH1::basisFunctionDeriv(unsigned int i, unsigned int jDir, const PointReferenceT& p) const {
        if(!(p==currentPoint_)){
		std::cout<<"WARNING: you are using a slow operator";
		return element_->basisFunctionDeriv(i,jDir,p);
	}
        return basisFunctionIndividualDerivatives_[i][jDir];
}


