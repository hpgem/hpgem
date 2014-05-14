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

#include"Integration/ReturnTrait1.hpp"
#include"Base/Element.hpp"
#include "Base/ShortTermStorageElementH1.hpp"
#include "Base/ShortTermStorageFaceH1.hpp"

namespace Base
{
    class HpgemUISimplified;
    
    HpgemUISimplified::HpgemUISimplified(unsigned int DIM,int polynomialOrder):HpgemUI(new GlobalData,new ConfigurationData(DIM,1,polynomialOrder,1)){}
    
    bool
    HpgemUISimplified::solve()
    {
        initialise();
        checkInitialisation();
        for (int i=0; i < HpgemUI::meshes_.size();i++)
        {
            HpgemUI::meshes_[i]->move(); // just for testing
        }
        doAllElementIntegration();
        doAllFaceIntegration();
        
        return true;
    }
    template<class T>
    struct ABC
    {
        void aas(T& a){}
    };
    
    void
    HpgemUISimplified::doAllFaceIntegration(unsigned int meshID)
    {
        bool useCache   = false;
        
        LinearAlgebra::Matrix fMatrixData;
        LinearAlgebra::NumericalVector fVectorData;
        FaceIntegralT   faceIntegral(useCache);
        faceIntegral.setStorageWrapper(new Base::ShortTermStorageFaceH1(meshes_[meshID]->dimension()));
        
        for (MeshManipulator::FaceIterator citFe = Base::HpgemUI::faceColBegin(); citFe != Base::HpgemUI::faceColEnd(); ++citFe)
        {
        	int n=(*citFe)->getPtrElementLeft()->getNrOfUnknows()*(*citFe)->getPtrElementLeft()->getNrOfBasisFunctions();
        	if((*citFe)->isInternal())
        		n+=(*citFe)->getPtrElementRight()->getNrOfUnknows()*(*citFe)->getPtrElementRight()->getNrOfBasisFunctions();
        	fMatrixData.resize(n,n);
        	fVectorData.resize(n);
            faceIntegral.integrate<LinearAlgebra::Matrix>((*citFe), this, fMatrixData);
            (*citFe)->setFaceMatrix(fMatrixData);
            faceIntegral.integrate<LinearAlgebra::NumericalVector>((*citFe), this, fVectorData);
            (*citFe)->setFaceVector(fVectorData);
        }
    }
    
    void
    HpgemUISimplified::doAllElementIntegration(unsigned int meshID)
    {
        unsigned int ndof = HpgemUI::configData_->numberOfBasisFunctions_*HpgemUI::configData_->numberOfUnknowns_;
        LinearAlgebra::Matrix  	eMatrixData(ndof, ndof);
        LinearAlgebra::NumericalVector eVectorData(ndof);
        
        bool isUseCache(false);
        Integration::ElementIntegral 	elIntegral(isUseCache);
        elIntegral.setStorageWrapper(new ShortTermStorageElementH1(meshes_[meshID]->dimension()));
        
        for (ElementIterator it=HpgemUI::meshes_[meshID]->elementColBegin(); it!= HpgemUI::meshes_[meshID]->elementColEnd(); ++it)
        {
            elIntegral.integrate<LinearAlgebra::Matrix>((*it), this, eMatrixData);
            (*it)->setElementMatrix(eMatrixData);
            elIntegral.integrate<LinearAlgebra::NumericalVector>((*it), this, eVectorData);
            (*it)->setElementVector(eVectorData);

            //cout << result;
            //cout<< "#####################################END of ELEMENT######"<<endl;
        }
    }
    
    bool
    HpgemUISimplified::checkInitialisation()
    {
        if (HpgemUI::meshes_.size()==0)
        {

            std::cerr << "Error no mesh created : You need to create at least one mesh to solve a problem" << std::endl;

        }
    }
}
