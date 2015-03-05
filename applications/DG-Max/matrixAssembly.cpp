
/*Alpha
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
#define hpGEM_INCLUDE_PETSC_SUPPORT

#include "matrixAssembly.hpp"
#include "BaseExtended.hpp"
#include "Base/HpgemUI.hpp"
#include "Base/HpgemUISimplified.hpp"
#include "Utilities/GlobalMatrix.hpp"
#include "Utilities/GlobalVector.hpp"
#include "petscksp.h"
#include "BaseExtended.hpp"
#include "Geometry/PointPhysical.hpp"
#include "Base/ConfigurationData.hpp"
#include "Integration/FaceIntegral.hpp"
#include "Integration/ElementIntegral.hpp"
#include "Base/FaceCacheData.hpp"
#include "Base/ElementCacheData.hpp"
#include "Base/GlobalData.hpp"
#include "Base/ConfigurationData.hpp"
#include "Base/Element.hpp"
#include "Base/Face.hpp"
#include "Base/ShortTermStorageElementHcurl.hpp"
#include "Base/ShortTermStorageFaceHcurl.hpp"



//cannot understand about the elementMatrixId and faceMatrisID as well as elementVectorId and faceVectorID
// since in DGMax, we  have n number of elements from o to 1 and the number n has to be supplied by the User


//Matrix and vector have been defined, but their storage has still to be done...also for SourceTerm and Initial Conditions , the integrand has to be rechecked

/*
void MatrixAssemblyIP::fillMatrics(hpGemExtensions* matrixContainer)
{
    time_t oldTime, newTime;
    time(&oldTime);
    std::cout<<"using an IP-DG method with penalty parameter"<<matrixContainer->getData()->StabCoeff->std::endl;
    
    //Calculation of number of local processors and the total amount of processor
    int TotalAmoutOfProcessors, localProcessorNumber, numberOfElements;
    MPI_Comm_size(PETSC_COMM_WORLD, &TotalAmoutOfProcessors);
    MPI_Comm_rank(PETSC_COMM_WORLD, &localProcessorNumber);
    numberOfElements = matrixContainer->getNumberOfElements(0);
    
    
    Utilities::GlobalPetscMatrix M_(HpgemUI::meshes_[0], 1, 0), S_(HpgemUI::meshes_[0], 1, 1);
    //No of element Matrices in M_ = 1 and No of face Matrices in M_ = 0
    //No of element Matrices in S_ = 1 and No of face Matrices in S_ = 1 :- still doubtfull about number of faceMatrices since it includes the boundary face term as well
    
    
    Utilities::GlobalPetscVector x_(HpgemUI::MESHES_[0]), RHS_(HpgemUI::meshes_[0]), derivative_(HpgemUI::MESHES_[0]);
    
    //Mass Matrix initiated
    
    Base::ShortTermStorageElementBase* localElement_;
    localElement_ = new Base::ShortTermStorageElementHcurl(3); // creating object for H curl transformation
    Integration::ElementIntegral elIntegral(false);
    elIntegral.setStorageWrapper(localElement_);
    
    LinearAlgebra::Matrix matrix(1, 1);
    LinearAlgebra::NumericalVector vector(1);
    
    time(&newTime);
    std::cout<<"Matrices getting filled, preparation took "<<difftime(newTime, oldTime)<<" seconds"<<std::endl;
    oldTime = newTime;
    
    for(hpGemUIExtentions::ElementIterator it=matrixContainer->elementColBegin(); it!=matrixContainer->elementColEnd(); ++it)
    {
        //The (*it)->getID has to be tested or not ???
        
        if(((*it)->getID())/(numberOfElements/TotalAmountOfProcessors +  1) == localProcessorNumber)
        {
            matrix.resize((*it)->getNrOfBasisFunctions(), (*it)->getNrOfBasisFunctions());
            elIntegral.integrate<LinearAlgebra::Matrix>((*it), &DGMax::elementMassIntegrand, matrix);
            if(matrixContainer->MHasToBeInverted_)
            {
                matrix.inverse(matrix);
            }
            int places[] = {(int)(*it)->getID()};
            
            for(int i = =; i < (*it)->getNrOfBasisFunctions()*
      
        }
    }
    
    //Mass Matrix ended
    
    
    //Stiffness Matrix Initiated
    for(hpGemUIExtentions::ElementIterator it=matrixContainer->elementColBegin(); it!=matrixContainer->elementColEnd(); ++it)
    {
        matrix.resize((*it)->getNrOfBasisFunctions(), (*it)->getNrOfBasisFunctions());
        elIntegral.integrate<LinearAlgebra::Matrix>((*it), &DGMax::elementStiffnessIntegrand, matrix);
        
    }
    
    Base::ShortTermStorageFaceBase* localFace_;
    localFace_ = new Base::ShortTermStorageFaceHcurl(3); // creating object for H curl transformation
    
    Integration::FaceIntegral faIntegral(false);
    
    faIntegral.setStorageWrapper(localFace_);
    
    for(hpGemUIExtentions::FaceIterator it = matrixContainer->faceColBegin(); it!= matrixContainer->faceColEnd(); ++it)
    {
        if((*it)->isInternal())
        {
            matrix.resize((*it)->getPtrElementLeft()->getNrOfBasisFunctions()+(*it)->getPtrElementRight()->getNrOfBasisFunctions(),(*it)->getPtrElementLeft()->getNrOfBasisFunctions()+(*it)
            getPtrElementRight()->getNrOfBasisFunctions());
            
            int places[] = {(int)(*it)->getPtrElementLeft()->getID(), (int)(*it)->getPtrElementRight()->getID()};
            faIntegral.integrate<LinearAlgebra::Matrix>(*it, &DGMax::faceStiffnessIntegrand, matrix);
        }
        else
        {
         matrix.resize((*it)->getPtrElementLeft()->getNrOfBasisFunctions(),(*it)->getPtrElementLeft()->getNrOfBasisFunctions());
        int places[]= {(int)(*it)->getPtrElementLeft()->getID()};
        faIntegral.integrate<LinearAlgebra::Matrix>(*it, &DGMax::faceStiffnessIntegrand, matrix);
        }
    }
    
    //Stiffness Matrix ended
    
    //face matrix in Stiffness matrix for IP part starts
    for(hpGemUIExtentions::FaceIterator it=matrixContainer->faceColBegin(); it!=matrixContainer->faceColEnd(); ++it)
    {
    if((*it)->isInternal())
    {
    matrix.resize((*it)->getPtrElementLeft()->getNrOfBasisFunctions()+(*it)->getPtrElementRight()->getNrOfBasisFunctions(),(*it)->getPtrElementLeft()->getNrOfBasisFunctions()+(*it)->getPtrElementRight()
    getNrOfBasisFunctions());
        
    int places[]= {(int)(*it)->getPtrElementLeft()->getID(),(int)(*it)->getPtrElementRight()->getID()};
    faIntegral.integrate<LinearAlgebra::Matrix>(*it,&DGMax::faceStiffnessIntegrandIP, matrix);
    }
    else
    {
        matrix.resize((*it)->getPtrElementLeft()->getNrOfBasisFunctions(),(*it)->getPtrElementLeft()->getNrOfBasisFunctions());
        int places[]= {(int)(*it)->getPtrElementLeft()->getID()};
        faIntegral.integrate<LinearAlgebra::Matrix>(*it,&DGMax::faceStiffnessIntegrandIP, matrix);
    }
    }
    //face integrand in Stiffness matrix for IP part ends
    
    //face integrand in RHS vector starts
    for(hpGemUIExtentions::ElementIterator it=matrixContainer->elementColBegin(); it!=matrixContainer->elementColEnd(); ++it)
    {
        vector.resize((*it)->getNrOfBasisFunctions());
        elIntegral.integrate<LinearAlgebra::NumericalVector>((*it), &DGMax::elementSpaceIntegrand, vector);
        int places[] = {(int)(*it->getID()};
    }
    
    //face integrand in RHS vector ends
    //face integrand in initial conditions starts
    for(hpGemUIExtentions::ElementIterator it=matrixContainer->elementColBegin(); it!=matrixContainer->elementColEnd(); ++it)
    {
        vector.resize((*it)->getNrOfBasisFunctions());
        elIntegral.integrate<LinearAlgebra::NumericalVector>((*it), &DGMax::elementSpaceIntegrand, vector);
    }
    //face integrand in initial conditions ends
    
    
}
*/
/*
void MatrixAssemblyBR::fillMatrics(hpGemExtensions* matrixContainer)
{
    LinearAlgebra::Matrix matrix(1,1);
    LinearAlgebra::NumericalVector vector(1);
    
    Base::ShortTermStorageElementBase* localElement_;
    localElement_ = new Base::ShortTermStorageElementHcurl(3); // creating object for H curl transformation
    Integration::ElementIntegral elIntegral(false);
    elIntegral.setStorageWrapper(localElement_);
    
    Base::ShortTermStorageFaceBase* localFace_;
    localFace_ = new Base::ShortTermStorageFaceHcurl(3); // creating object for H curl transformation
    Integration::FaceIntegral faIntegral(false);
    faIntegral.setStorageWrapper(localFace_);
    
    //Element Mass matrix definition started
    
    for(hpGemUIExtentions::ElementIterator it=matrixContainer->elementColBegin(); it!=matrixContainer->elementColEnd(); ++it)
    {
        matrix.resize((*it)->getNrOfBasisFunctions(), (*it)->getNrOfBasisFunctions());
        elIntegral.integrate<LinearAlgebra::Matrix>((*it), &DGMax::elementMassIntegrand, matrix);
    }
    
    
    //Element Mass matrix definition ended
    
    //Element stiffness matrix definition started
    for(hpGemUIExtentions::ElementIterator it = matrixContainer->elementColBegin(); it! = matrixContainer->elementColEnd(); ++it)
    {
        matrix.resize((*it)->getNrOfBasisFunctions(), (*it)->getNrOfBasisFunctions());
        elIntegral.integrate<LinearAlgebra::Matrix>((*it), &DGMax::elementStiffnessIntegrand, matrix);
    }
    
    for(hpGemUIExtentions::FaceIterator it = matrixContainer->faceColBegin(); it! = matrixContainer->faceColEnd(); ++it))
    {
        if(*it)->isInternal())
        {
            matrix.resize((*it)->getPtrElementLeft()->getNrOfBasisFunctions()+(*it)->getPtrElementRight()->getNrOfBasisFunctions(), (*it)->getPtrElementLeft()->getNrOfBasisFunctions()+(*it)->getPtrElementRight()->getNrOfBasisFunctions());
            int places[] = {(int)(*it)->getPtrElementLeft()->getID(), (int)(*it)->getPtrElementRight()->getID()};
            faIntegral.integrate<LinearAlgebra::Matrix>(*it, &DGMax::faceStiffnessIntegrand, matrix);
        }
        else
        {
            matrix.resize((*it)->getPtrElementLeft()->getNrOfBasisFunctions(), (*it)->getPtrElementLeft()->getNrOfBasisFunctions());
            int places[] = {(int)(*it)->getPtrElementLeft()->getID()};
            faIntegral.integrate<LinearAlgebra::Matrix>(*it, &DGMax::faceStiffnessIntegrand, matrix);
  
        }
    }
    
    //Element stiffness matrix definition ended
    
    
}


            
void MatrixAssemblyIP::CompleteElementIntegrationIP(DGMax* matrixContainer)
        {
            Base::ShortTermStorageElementBase* localElement_;
            localElement_ = new Base::ShortTermStorageElementHcurl(3); // creating object for H curl transformation
            Integration::ElementIntegral elIntegral(false);
            elIntegral.setStorageWrapper(localElement_);
            
            LinearAlgebra::Matrix matrix1(1, 1), matrix2(1, 1);
            LinearAlgebra::NumericalVector vector1(1), vector2(1);
            
            for(hpGemUIExtentions::ElementIterator it=matrixContainer->elementColBegin(); it!=matrixContainer->elementColEnd(); ++it)
            {
                matrix1.resize((*it)->getNrOfBasisFunctions(), (*it)->getNrOfBasisFunctions());
                elIntegral.integrate<LinearAlgebra::Matrix>((*it), &(matrixContainer->elementMassIntegrand), matrix1);
                (*it)->setElementMatrix(matrix1, 0);
                
                matrix2.resize((*it)->getNrOfBasisFunctions(), (*it)->getNrOfBasisFunctions());
                elIntegral.integrate((*it), &(matrixContainer->elementStiffnessIntegrand), matrix2);
                (*it)->setElementMatrix(matrix2, 1);
                
                vector1.resize((*it)->getNrOfBasisFunctions());
                elIntegral.integrate((*it), &(matrixContainer->elementSpaceIntegrand), vector1);
                (*it)->setElementVector(vector1, 0);
                
                vector2.resize((*it)->getNrOfBasisFunctions());
                elIntegral.integrate((*it), &(matrixContainer->elementSpaceIntegrand), vector2);
                (*it)->setElementVector(vector2, 1);
            }
        }
            
            
void MatrixAssemblyIP::CompleteFaceIntegrationIP(DGMax* matrixContainer)
        {
            Base::ShortTermStorageFaceBase* localFace_;
            localFace_ = new Base::ShortTermStorageFaceHcurl(3); // creating object for H curl transformation
            Integration::FaceIntegral faIntegral(false);
            faIntegral.setStorageWrapper(localFace_);
            
            LinearAlgebra::Matrix matrix1(1, 1), matrix2(1, 1);
            LinearAlgebra::NumericalVector vector(1);
            
            for(hpGemUIExtentions::FaceIterator it=matrixContainer->faceColBegin(); it!=matrixContainer->faceColEnd(); ++it)
            {
                if((*it)->isInternal())
                   {
                       matrix1.resize((*it)->getPtrElementLeft()->getNrOfBasisFunctions() + (*it)->getPtrElementRight()->getNrOfBasisFunctions(), (*it)->getPtrElementLeft()->getNrOfBasisFunctions() + (*it)->getPtrElementRight()->getNrOfBasisFunctions());
                       
                       faIntegral.integrate((*it), &(matrixContainer->faceStiffnessIntegrand), matrix1);
                       
                       (*it)->setFaceMatrix(matrix1, 1);
                       
                       matrix2.resize((*it)->getPtrElementLeft()->getNrOfBasisFunctions() + (*it)->getPtrElementRight()->getNrOfBasisFunctions(), (*it)->getPtrElementLeft()->getNrOfBasisFunctions() + (*it)->getPtrElementRight()->getNrOfBasisFunctions());
                       
                       faIntegral.integrate((*it), &(matrixContainer->faceStiffnessIntegrandIP), matrix2);
                       
                       (*it)->setFaceMatrix(matrix2, 1);
                   }
                   
                else
                   {
                       matrix1.resize((*it)->getPtrElementLeft()->getNrOfBasisFunctions(), (*it)->getPtrElementLeft()->getNrOfBasisFunctions());
                       faIntegral.integrate((*it), &(matrixContainer->faceStiffnessIntegrand), matrix1);
                       
                       (*it)->setFaceMatrix(matrix1, 1);
                       
                       matrix2.resize((*it)->getPtrElementLeft()->getNrOfBasisFunctions(), (*it)->getPtrElementLeft()->getNrOfBasisFunctions());
                       
                       faIntegral.integrate((*it), &(matrixContainer->faceStiffnessIntegrandIP), matrix2);
                       
                       (*it)->setFaceMatrix(matrix2, 1);
                       
                       vector.resize((*it)->getPtrElementLeft()->getNrOfBasisFunctions());
                       
                       faIntegral.integrate((*it), &(matrixContainer->faceSpaceIntegrandIP), vector);
                       
                       (*it)->setFaceVector(vector, 0);

                   }

            }
        }

void MatrixAssemblyIP::fillMatrices(hpGemUIExtentions* matrixContainer1)
        {
            DGMax* matrixContainer;
            matrixContainer = reinterpret_cast<DGMax*>(&matrixContainer1);
            
            // to assign base class pointer to derived class pointer...static_cast didnt work
            
            CompleteElementIntegrationIP(matrixContainer);
            CompleteFaceIntegrationIP(matrixContainer);
    
        }


*/







