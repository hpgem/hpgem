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

#include "fillMatrices.h"
#include "BaseExtended.h"
#include "Geometry/PointPhysical.h"
#include "Base/ConfigurationData.h"
#include "Integration/ElementIntegral.h"
#include "Base/FaceCacheData.h"
#include "Base/ElementCacheData.h"
#include "Base/ShortTermStorageElementHCurl.h"
#include "Base/ShortTermStorageFaceHCurl.h"
/*
 //  this is where you choose the solution of your problem
 //  this will only have an effect on the accuracy of your error estimates
 //  as a temporary solution remember to also update the exact solution in DG-Max.cpp
 
 void matrixFiller::initialExactSolution(const PointPhysicalT& p, LinearAlgebra::NumericalVector &ret){
 ret[0]=sin(M_PI*2*p[1])*sin(M_PI*2*p[2]);
 ret[1]=sin(M_PI*2*p[2])*sin(M_PI*2*p[0]);
 ret[2]=sin(M_PI*2*p[0])*sin(M_PI*2*p[1]);
 //        ret[0]=sin(M_PI*p[1])*sin(M_PI*p[2]);
 //        ret[1]=sin(M_PI*p[2])*sin(M_PI*p[0]);
 //        ret[2]=sin(M_PI*p[0])*sin(M_PI*p[1]);
 //            ret[0]=p[0]*(1-p[0]);
 //            ret[1]=0;
 // 	  		  ret[2]=0;
 }
 
 
 // Computes element contributions to the stiffness matrix i.e. (nabla x phi_i) * (nabla x phi_j)
 // returns the contibutions at this gauss point to the entire element matrix in one go
 
 void matrixFiller::elementIntegrand(const ElementT* element, const Geometry::PointReference& p, LinearAlgebra::Matrix& ret){
 //cout<<"\nIn the element integrand for the stiffness matrix for element id: "<<element->getID();
 ret.resize(element->getNrOfBasisFunctions(),element->getNrOfBasisFunctions());
 //ElementInfos* info = static_cast<ElementInfos*> (const_cast<ElementT*>(element)->getUserData());
 LinearAlgebra::NumericalVector phi_i(3),phi_j(3);
 //std::vector<LinearAlgebra::NumericalVector> functionCurls;
 //info->makeFunctionCurlsVector(element,p,functionCurls);
 for(int i=0;i<element->getNrOfBasisFunctions();++i){
 //phi_i=functionCurls[i];
 element->basisFunctionCurl(i, p, phi_i);
 for(int j=i;j<element->getNrOfBasisFunctions();++j){
 //phi_j=functionCurls[j];
 element->basisFunctionCurl(j, p, phi_j);
 ret(i,j)=phi_i[0]*phi_j[0]+phi_i[1]*phi_j[1]+phi_i[2]*phi_j[2];
 ret(j,i)=ret(i,j);
 }
 }
 }
 
 void matrixFiller::elementIntegrand(const Base::Element* element, const Geometry::PointReference& p, LinearAlgebra::NumericalVector& ret){
 ret.resize(element->getNrOfBasisFunctions());
 //ElementInfos* info = static_cast<ElementInfos*>(const_cast<ElementT*>(element)->getUserData());
 PointPhysicalT pPhys(3);
 element->referenceToPhysical(p,pPhys);
 LinearAlgebra::NumericalVector val(3),phi(3);
 //std::vector<LinearAlgebra::NumericalVector> functionValues;
 //info->makeFunctionValuesVector(element,p,functionValues);
 sourceTerm(pPhys,val);
 for(int i=0; i<element->getNrOfBasisFunctions(); ++i) {
 //phi=functionValues[i];
 element->basisFunction(i, p, phi);
 ret(i)=phi[0]*val[0]+phi[1]*val[1]+phi[2]*val[2];
 }
 }
 
 void matrixFillerIP::faceIntegrand(const FaceT* face, const LinearAlgebra::NumericalVector& normal, const Geometry::PointReference& p, LinearAlgebra::Matrix& ret){
 //cout<<"\nIn the face integrand for the stiffness matrix (IP-only part) for element id: "<<face->getPtrElementLeft()->getID();
 //int M = face->getNrOfBasisFunctions();
 
 int M = face->getPtrElementLeft()->getNrOfBasisFunctions() + face->getPtrElementRight()->getNrOfBasisFunctions();
 ret.resize(M, M);
 
 LinearAlgebra::NumericalVector phi_i(3), phi_j(3);
 for(int i = 0; i < M; ++i)
 {
 face->basisFunctionNormal(i, normal, p, phi_i);
 
 for (int j = i; j < M; ++j)
 {
 face->basisFunctionNormal(j, normal, p, phi_j);
 ret(i, j) = stabCoeff_ * (phi_i[0] * phi_j[0] + phi_i[1] * phi_j[1] + phi_i[2] * phi_j[2]);
 ret(j, i) = ret(i, j);
 
 }
 }
 
 }
 
 
 
 
 
 // Computes the bits of the face contributions that are only used in the IP method
 // i.e. eta_F( (n x phi_i) * (n x phi_j) )
 // returns the contibutions at this gauss point to the entire face matrix in one go
 
 void matrixFillerIP::faceIntegrand(const FaceT* face, const LinearAlgebra::NumericalVector& normal, const Geometry::PointReference& p, LinearAlgebra::Matrix& ret){
 //cout<<"\nIn the face integrand for the stiffness matrix (IP-only part) for element id: "<<face->getPtrElementLeft()->getID();
 ElementT* right;
 ElementT* left=const_cast<ElementT*>(face->getPtrElementLeft());
 //ElementInfos* leftInfo = static_cast<ElementInfos*> (left->getUserData());
 //ElementInfos* rightInfo;
 PointElementReferenceT pLeft(3),pRight(3);
 face->mapRefFaceToRefElemL(p,pLeft);
 //std::vector<LinearAlgebra::NumericalVector> leftValues,rightValues;
 //leftInfo->makeFunctionValuesVector(left,pLeft,leftValues);
 LinearAlgebra::NumericalVector normedNormal(3);
 normedNormal[0] = (normal*(1/Base::L2Norm(normal)))[0];
 normedNormal[1] = (normal*(1/Base::L2Norm(normal)))[1];
 normedNormal[2] = (normal*(1/Base::L2Norm(normal)))[2];
 int leftSize=left->getNrOfBasisFunctions();
 int dimension=left->getNrOfBasisFunctions();
 LinearAlgebra::NumericalVector phi_i(3),phi_j(3),dummy(3);
 if(face->isInternal()){
 //cout<<" and element id: "<<face->getPtrElementRight()->getID();
 right=const_cast<ElementT*>(face->getPtrElementRight());
 face->mapRefFaceToRefElemR(p,pRight);
 //rightInfo = static_cast<ElementInfos*> (right->getUserData());
 dimension+=right->getNrOfBasisFunctions();
 //rightInfo->makeFunctionValuesVector(right,pRight,rightValues);
 }
 ret.resize(dimension,dimension);
 for(int i=0;i<dimension;++i){
 if(i<leftSize){
 //dummy=leftValues[i];
 left->basisFunction(i, pLeft, dummy);
 }else{
 //dummy=rightValues[i-leftSize];
 right->basisFunction(i-leftSize, pRight, dummy);
 dummy*=-1;
 }
 OuterProduct(normedNormal,dummy,phi_i);
 for(int j=i;j<dimension;++j){
 if(j<leftSize){
 //dummy=leftValues[j];
 left->basisFunction(j, pLeft, dummy);
 }else{
 //dummy=rightValues[j-leftSize];
 right->basisFunction(j-leftSize, pRight, dummy);
 dummy*=-1;
 }
 OuterProduct(normedNormal,dummy,phi_j);
 ret(i,j)=stabCoeff_*(phi_i[0]*phi_j[0]+phi_i[1]*phi_j[1]+phi_i[2]*phi_j[2]);
 ret(j,i)=ret(i,j);
 }
 }
 }
 
 void matrixFillerIP::faceIntegrand(const Base::Face* face, const LinearAlgebra::NumericalVector& normal, const Geometry::PointReference& p, LinearAlgebra::NumericalVector& ret){
 ElementT* left=const_cast<ElementT*>(face->getPtrElementLeft());
 ElementInfos* info = static_cast<ElementInfos*>(left->getUserData());
 ret.resize(left->getNrOfBasisFunctions());
 Geometry::PointReference PLeft(3);
 face->mapRefFaceToRefElemL(p,PLeft);
 //std::vector<LinearAlgebra::NumericalVector> functionValues,functionCurls;
 //info->makeFunctionValuesVector(left,PLeft,functionValues);
 //info->makeFunctionCurlsVector(left,PLeft,functionCurls);
 PointPhysicalT PPhys(3);
 left->referenceToPhysical(PLeft,PPhys);
 LinearAlgebra::NumericalVector normedNormal(3);
 normedNormal[0] = (normal*(1/Base::L2Norm(normal)))[0];
 normedNormal[1] = (normal*(1/Base::L2Norm(normal)))[1];
 normedNormal[2] = (normal*(1/Base::L2Norm(normal)))[2];
 LinearAlgebra::NumericalVector val(3),phi(3),phi_curl(3),dummy(3);
 boundaryConditions(PPhys,dummy);//assumes the initial conditions and the boundary conditions match
 OuterProduct(normedNormal,dummy,val);
 for(int i=0; i<left->getNrOfBasisFunctions(); ++i) {
 //dummy=functionValues[i];
 //phi_curl=functionCurls[i];
 left->basisFunction(i, PLeft, dummy);
 left->basisFunctionCurl(i, PLeft, phi_curl);
 OuterProduct(normedNormal,dummy,phi);
 ret(i)=-(phi_curl[0]*val[0]+phi_curl[1]*val[1]+phi_curl[2]*val[2])+
 stabCoeff_*(phi[0]*val[0]+phi[1]*val[1]+phi[2]*val[2]);
 }
 }
 
 
 void matrixFillerIP::faceIntegrand(const Base::Face* face, const LinearAlgebra::NumericalVector& normal, const Geometry::PointReference& p, LinearAlgebra::NumericalVector& ret)
 {
 
 ElementT* left=const_cast<ElementT*>(face->getPtrElementLeft());
 Geometry::PointReference PLeft(3);
 face->mapRefFaceToRefElemL(p,PLeft);
 
 
 PointPhysicalT PPhys(3);
 left->referenceToPhysical(PLeft,PPhys);
 LinearAlgebra::NumericalVector normedNormal(3);
 
 normedNormal[0] = (normal * (1/Base::L2Norm(normal)))[0];
 normedNormal[1] = (normal * (1/Base::L2Norm(normal)))[1];
 normedNormal[2] = (normal * (1/Base::L2Norm(normal)))[2];
 
 LinearAlgebra::NumericalVector val(3),phi(3),phi_curl(3),dummy(3);
 
 boundaryConditions(PPhys,dummy);//assumes the initial conditions and the boundary conditions match
 
 OuterProduct(normedNormal,dummy,val);
 int n = face->getPtrElementLeft()->getNrOfBasisFunctions();
 ret.resize(n);
 
 for(int i = 0; i < n; ++i)
 {
 face->basisFunctionNormal(i, normal, p, phi);
 face->basisFunctionCurl(i, p, phi_curl);
 
 ret(i) = -(phi_curl[0] * val[0] + phi_curl[1] * val[1] + phi_curl[2] * val[2]) + stabCoeff_ * (phi[0] * val[0] + phi[1] * val[1] + phi[2] * val[2]);
 
 }
 
 }
 
 
 
 void matrixFillerIP::fillMatrixes(hpGemUIExtentions* matrixContainer){
 time_t oldTime,newTime;
 time(&oldTime);
 std::cout<<"using an IP-DG method with penalty parameter "<<matrixContainer->getData()->StabCoeff_<<std::endl;
 //learn your own processor number and the total amount of processors so there is no double work done while filling the matrices
 int TotalAmountOfProcessors,localProcessorNumber,numberOfElements;
 MPI_Comm_size(PETSC_COMM_WORLD,&TotalAmountOfProcessors);
 MPI_Comm_rank(PETSC_COMM_WORLD,&localProcessorNumber);
 numberOfElements=matrixContainer->getNumberOfElements(0);
 
 // 	if(matrixContainer->getNumberOfElements(0)%TotalAmountOfProcessors!=0){
 // 	    cout<<"WARNING: the case where the number of elements is not a multiple of the number of processors has not been toroughly tested; gliches may occur!";'
 // 	}
 
 matrixContainer->ierr_=VecSetSizes(matrixContainer->x_,PETSC_DECIDE,matrixContainer->getData()->numberOfUnknowns_);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=VecSetSizes(matrixContainer->RHS_,PETSC_DECIDE,matrixContainer->getData()->numberOfUnknowns_);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=VecSetSizes(matrixContainer->derivative_,PETSC_DECIDE,matrixContainer->getData()->numberOfUnknowns_);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=MatSetSizes(matrixContainer->M_,PETSC_DECIDE,PETSC_DECIDE,matrixContainer->getData()->numberOfUnknowns_,matrixContainer->getData()->numberOfUnknowns_);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=MatSetSizes(matrixContainer->S_,PETSC_DECIDE,PETSC_DECIDE,matrixContainer->getData()->numberOfUnknowns_,matrixContainer->getData()->numberOfUnknowns_);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=VecSetBlockSize(matrixContainer->x_,matrixContainer->getConfigData()->numberOfBasisFunctions_);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=VecSetBlockSize(matrixContainer->RHS_,matrixContainer->getConfigData()->numberOfBasisFunctions_);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=VecSetBlockSize(matrixContainer->derivative_,matrixContainer->getConfigData()->numberOfBasisFunctions_);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=MatSetBlockSize(matrixContainer->M_,matrixContainer->getConfigData()->numberOfBasisFunctions_);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=MatSetBlockSize(matrixContainer->S_,matrixContainer->getConfigData()->numberOfBasisFunctions_);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=VecSetUp(matrixContainer->x_);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=VecSetUp(matrixContainer->RHS_);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=VecSetUp(matrixContainer->derivative_);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=MatSetUp(matrixContainer->M_);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=MatSetUp(matrixContainer->S_);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 
 //once it is finished this code should make use of the global assembly interface
 LinearAlgebra::Matrix matrix(1,1);
 LinearAlgebra::NumericalVector vector(1);
 PetscScalar *tempComplexArray= new PetscScalar[4*matrixContainer->getConfigData()->numberOfBasisFunctions_*matrixContainer->getConfigData()->numberOfBasisFunctions_];
 //hpGemUIExtentions::ElementFunction elF = &hpGemUIExtentions::elementMassIntegrand;
 
 Base::ShortTermStorageElementBase* localElement_;
 localElement_ = new Base::ShortTermStorageElementHCurl(3); // creating object for H curl transformation
 Integration::ElementIntegral elIntegral(false);
 elIntegral.setStorageWrapper(localElement_);
 
 time(&newTime);
 std::cout<<"filling the matrices; preparation took "<<difftime(newTime,oldTime)<<" seconds"<<std::endl;
 oldTime=newTime;
 for(hpGemUIExtentions::ElementIterator it=matrixContainer->elementColBegin(); it!=matrixContainer->elementColEnd(); ++it) {
 if(((*it)->getID())/(numberOfElements/TotalAmountOfProcessors+1)==localProcessorNumber) {
 matrix.resize((*it)->getNrOfBasisFunctions(),(*it)->getNrOfBasisFunctions());
 elIntegral.integrate<LinearAlgebra::Matrix>((*it),matrixContainer,matrix);
 //MATRIX INVERSE IMPLEMENTED TODAY
 if(matrixContainer->MHasToBeInverted_) {
 matrix.inverse(matrix);
 }

 int places[]= {(int)(*it)->getID()};
 for(int i=0; i<(*it)->getNrOfBasisFunctions()*(*it)->getNrOfBasisFunctions(); ++i) {
 tempComplexArray[i]=matrix[i];//the real<->complex conflict should probably be solved in some other way
 }
 matrixContainer->ierr_=MatSetValuesBlocked(matrixContainer->M_,1,places,1,places,&tempComplexArray[0],INSERT_VALUES);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 }//the & and the [0] are still here for easy conversion back to &matrix[0]
 }
 //PETSc can start assembling M if it wants to; we wont set any extra antries anymore
 matrixContainer->ierr_=MatAssemblyBegin(matrixContainer->M_,MAT_FINAL_ASSEMBLY);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 time(&newTime);
 std::cout<<"filling the mass matrix took "<<difftime(newTime,oldTime)<<" seconds"<<std::endl;
 oldTime=newTime;
 //elF=&hpGemUIExtentions::elementStiffnessIntegrand;
 for(hpGemUIExtentions::ElementIterator it=matrixContainer->elementColBegin(); it!=matrixContainer->elementColEnd(); ++it) {
 if(((*it)->getID())/(numberOfElements/TotalAmountOfProcessors+1)==localProcessorNumber) {
 matrix.resize((*it)->getNrOfBasisFunctions(),(*it)->getNrOfBasisFunctions());
 elIntegral.integrate<LinearAlgebra::Matrix>((*it),this,matrix);
 int places[]= {(int)(*it)->getID()};
 for(int i=0; i<(*it)->getNrOfBasisFunctions()*(*it)->getNrOfBasisFunctions(); ++i) {
 tempComplexArray[i]=matrix[i];
 }
 matrixContainer->ierr_=MatSetValuesBlocked(matrixContainer->S_,1,places,1,places,&tempComplexArray[0],ADD_VALUES);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 }
 }
 //hpGemUIExtentions::FaceFunction faF = &hpGemUIExtentions::faceIntegrand;
 
 Base::ShortTermStorageFaceBase* localFace_;
 localFace_ = new Base::ShortTermStorageFaceHCurl(3); // creating object for H curl transformation
 
 Integration::FaceIntegral faIntegral(false);
 
 faIntegral.setStorageWrapper(localFace_);
 
 for(hpGemUIExtentions::FaceIterator it=matrixContainer->faceColBegin(); it!=matrixContainer->faceColEnd(); ++it) {
 //faces dont have an ID; pick an arbitrary processor to do the work
 if(((*it)->getPtrElementLeft()->getID())/(numberOfElements/TotalAmountOfProcessors+1)==localProcessorNumber) {
 if((*it)->isInternal()) {
 matrix.resize((*it)->getPtrElementLeft()->getNrOfBasisFunctions()+(*it)->getPtrElementRight()->getNrOfBasisFunctions(),(*it)->getPtrElementLeft()->getNrOfBasisFunctions()+(*it)->getPtrElementRight()->getNrOfBasisFunctions());
 int places[]= {(int)(*it)->getPtrElementLeft()->getID(),(int)(*it)->getPtrElementRight()->getID()};
 faIntegral.integrate<LinearAlgebra::Matrix>(*it,matrixContainer,matrix);
 for(int i=0; i<4*(*it)->getPtrElementLeft()->getNrOfBasisFunctions()*(*it)->getPtrElementLeft()->getNrOfBasisFunctions(); ++i) {
 tempComplexArray[i]=matrix[i];
 }
 matrixContainer->ierr_=MatSetValuesBlocked(matrixContainer->S_,2,places,2,places,&tempComplexArray[0],ADD_VALUES);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 } else {
 matrix.resize((*it)->getPtrElementLeft()->getNrOfBasisFunctions(),(*it)->getPtrElementLeft()->getNrOfBasisFunctions());
 int places[]= {(int)(*it)->getPtrElementLeft()->getID()};
 faIntegral.integrate<LinearAlgebra::Matrix>(*it,matrixContainer,matrix);
 for(int i=0; i<(*it)->getPtrElementLeft()->getNrOfBasisFunctions()*(*it)->getPtrElementLeft()->getNrOfBasisFunctions(); ++i) {
 tempComplexArray[i]=matrix[i];
 }
 matrixContainer->ierr_=MatSetValuesBlocked(matrixContainer->S_,1,places,1,places,&tempComplexArray[0],ADD_VALUES);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 }
 }
 }
 //faF = &hpGemUIExtentions::faceIntegrandIPPart;
 for(hpGemUIExtentions::FaceIterator it=matrixContainer->faceColBegin(); it!=matrixContainer->faceColEnd(); ++it) {
 if(((*it)->getPtrElementLeft()->getID())/(numberOfElements/TotalAmountOfProcessors+1)==localProcessorNumber) {
 if((*it)->isInternal()) {
 matrix.resize((*it)->getPtrElementLeft()->getNrOfBasisFunctions()+(*it)->getPtrElementRight()->getNrOfBasisFunctions(),(*it)->getPtrElementLeft()->getNrOfBasisFunctions()+(*it)->getPtrElementRight()->getNrOfBasisFunctions());
 int places[]= {(int)(*it)->getPtrElementLeft()->getID(),(int)(*it)->getPtrElementRight()->getID()};
 faIntegral.integrate<LinearAlgebra::Matrix>(*it,this,matrix);
 for(int i=0; i<4*(*it)->getPtrElementLeft()->getNrOfBasisFunctions()*(*it)->getPtrElementLeft()->getNrOfBasisFunctions(); ++i) {
 tempComplexArray[i]=matrix[i];
 }
 matrixContainer->ierr_=MatSetValuesBlocked(matrixContainer->S_,2,places,2,places,&tempComplexArray[0],ADD_VALUES);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 } else {
 matrix.resize((*it)->getPtrElementLeft()->getNrOfBasisFunctions(),(*it)->getPtrElementLeft()->getNrOfBasisFunctions());
 int places[]= {(int)(*it)->getPtrElementLeft()->getID()};
 faIntegral.integrate<LinearAlgebra::Matrix>(*it,this,matrix);
 for(int i=0; i<(*it)->getPtrElementLeft()->getNrOfBasisFunctions()*(*it)->getPtrElementLeft()->getNrOfBasisFunctions(); ++i) {
 tempComplexArray[i]=matrix[i];
 }
 matrixContainer->ierr_=MatSetValuesBlocked(matrixContainer->S_,1,places,1,places,&tempComplexArray[0],ADD_VALUES);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 }
 }
 }
 time(&newTime);
 std::cout<<"filling the stiffness matrix took "<<difftime(newTime,oldTime)<<" seconds"<<std::endl;
 oldTime=newTime;
 matrixContainer->ierr_=MatAssemblyBegin(matrixContainer->S_,MAT_FINAL_ASSEMBLY);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 //elF=&hpGemUIExtentions::sourceTerm;
 for(hpGemUIExtentions::ElementIterator it=matrixContainer->elementColBegin(); it!=matrixContainer->elementColEnd(); ++it) {
 if(((*it)->getID())/(numberOfElements/TotalAmountOfProcessors+1)==localProcessorNumber) {
 vector.resize((*it)->getNrOfBasisFunctions());
 elIntegral.integrate<LinearAlgebra::NumericalVector>((*it),this,vector);
 int places[]= {(int)(*it)->getID()};
 for(int i=0; i<(*it)->getNrOfBasisFunctions(); ++i) {
 tempComplexArray[i]=vector[i];
 }
 matrixContainer->ierr_=VecSetValuesBlocked(matrixContainer->RHS_,1,places,&tempComplexArray[0],ADD_VALUES);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 }
 }
 //elF=&hpGemUIExtentions::initialConditions;
 for(hpGemUIExtentions::ElementIterator it=matrixContainer->elementColBegin(); it!=matrixContainer->elementColEnd(); ++it) {
 if(((*it)->getID())/(numberOfElements/TotalAmountOfProcessors+1)==localProcessorNumber) {
 vector.resize((*it)->getNrOfBasisFunctions());
 elIntegral.integrate<LinearAlgebra::NumericalVector>((*it),matrixContainer,vector);
 int places[]= {(int)(*it)->getID()};
 for(int i=0; i<(*it)->getNrOfBasisFunctions(); ++i) {
 tempComplexArray[i]=vector[i];
 }
 matrixContainer->ierr_=VecSetValuesBlocked(matrixContainer->x_,1,places,&tempComplexArray[0],ADD_VALUES);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 }
 }
 matrixContainer->ierr_=VecAssemblyBegin(matrixContainer->x_);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 //elF=&hpGemUIExtentions::initialConditionsDeriv;
 for(hpGemUIExtentions::ElementIterator it=matrixContainer->elementColBegin(); it!=matrixContainer->elementColEnd(); ++it) {
 if(((*it)->getID())/(numberOfElements/TotalAmountOfProcessors+1)==localProcessorNumber) {
 vector.resize((*it)->getNrOfBasisFunctions());
 elIntegral.integrate((*it),elF,vector,matrixContainer);
 int places[]= {(*it)->getID()};
 for(int i=0; i<(*it)->getNrOfBasisFunctions(); ++i) {
 tempComplexArray[i]=matrix[i];
 }
 matrixContainer->ierr_=VecSetValuesBlocked(matrixContainer->derivative_,1,places,&tempComplexArray[0],ADD_VALUES);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 }
 }
 matrixContainer->ierr_=VecAssemblyBegin(matrixContainer->derivative_);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 //faF = &hpGemUIExtentions::sourceTermBoundaryIP;
 for(hpGemUIExtentions::FaceIterator it=matrixContainer->faceColBegin(); it!=matrixContainer->faceColEnd(); ++it) {
 if(((*it)->getPtrElementLeft()->getID())/(numberOfElements/TotalAmountOfProcessors+1)==localProcessorNumber) {
 if((*it)->isInternal()) {
 //internal faces dont produce boundary contributions to the RHS
 } else {
 vector.resize((*it)->getPtrElementLeft()->getNrOfBasisFunctions());
 int places[]= {(int)(*it)->getPtrElementLeft()->getID()};
 faIntegral.integrate<LinearAlgebra::NumericalVector>(*it,this,vector);
 for(int i=0; i<(*it)->getPtrElementLeft()->getNrOfBasisFunctions(); ++i) {
 tempComplexArray[i]=vector[i];
 }
 matrixContainer->ierr_=VecSetValuesBlocked(matrixContainer->RHS_,1,places,&tempComplexArray[0],ADD_VALUES);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 }
 }
 }
 time(&newTime);
 std::cout<<"filling the vectors took "<<difftime(newTime,oldTime)<<" seconds"<<std::endl;
 oldTime=newTime;
 matrixContainer->ierr_=VecAssemblyBegin(matrixContainer->RHS_);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 
 //this functions guarantees assembled matrices so lock untill assembly is finished
 matrixContainer->ierr_=VecAssemblyEnd(matrixContainer->x_);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=VecAssemblyEnd(matrixContainer->derivative_);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=MatAssemblyEnd(matrixContainer->M_,MAT_FINAL_ASSEMBLY);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=MatAssemblyEnd(matrixContainer->S_,MAT_FINAL_ASSEMBLY);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=VecAssemblyEnd(matrixContainer->RHS_);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 time(&newTime);
 std::cout<<"communication took "<<difftime(newTime,oldTime)<<" extra seconds"<<std::endl;
 delete[] tempComplexArray;
 }
 
 void matrixFillerBR::faceIntegrand(const FaceT* face, const  LinearAlgebra::NumericalVector& normal, const Geometry::PointReference& p, LinearAlgebra::Matrix& ret){
 //cout<<"\nIn the face integrand for the stiffness matrix (BR-only part) for element id: "<<face->getPtrElementLeft()->getID();
 ElementT* right;
 ElementT* left=const_cast<ElementT*>(face->getPtrElementLeft());
 ElementInfos* leftInfo = static_cast<ElementInfos*> (left->getUserData());
 ElementInfos* rightInfo;
 
 //PointElementReferenceT pLeft(3),pRight(3);
 //face->mapRefFaceToRefElemL(p,pLeft);
 
 double localepsilon;
 
 if(face->isInternal())
 {
 //cout<<" and element id: "<<face->getPtrElementRight()->getID();
 //face->mapRefFaceToRefElemR(p,pRight);
 right=const_cast<ElementT*>(face->getPtrElementRight());
 rightInfo = static_cast<ElementInfos*> (right->getUserData());
 }
 
 int n = face->getPtrElementLeft()->getNrOfBasisFunctions();
 int M = face->getPtrElementLeft()->getNrOfBasisFunctions() + face->getPtrElementRight()->getNrOfBasisFunctions();
 //int M = face->getNrOfBasisFunctions();
 
 ret.resize(M, M);
 
 LinearAlgebra::NumericalVector phi_i(3),phi_j(3);
 
 for(int i = 0; i < M; ++i)
 {
 if(i < n)
 {
 localepsilon=leftInfo->epsilon_;
 }
 else
 {
 localepsilon=rightInfo->epsilon_;
 }
 
 localepsilon=leftInfo->epsilon_;
 face->basisFunction(i, p, phi_i);
 //std::cout<<phi_i<<std::endl;
 for (int j = 0; j < M; ++j)
 {
 face->basisFunctionNormal(j, normal, p, phi_j);
 
 ret(j,i) = (face->isInternal() ? 1 : 2) * (phi_i[0] * phi_j[0] + phi_i[1] * phi_j[1] + phi_i[2] * phi_j[2]) * sqrt(localepsilon);
 //std::cout<<ret<<std::endl;
 }
 }
 
 }
 
 
 
 
 //  Computes the bits of the face contributions that are only used in the BR formulation
 //  more accurately only returns phi_i * (n x phi_j) the matrix product should be done elsewhere
 //  returns the contibutions at this gauss point to the entire face matrix in one go
 
 
 void matrixFillerBR::faceIntegrand(const FaceT* face, const  LinearAlgebra::NumericalVector& normal, const Geometry::PointReference& p, LinearAlgebra::Matrix& ret){
 //cout<<"\nIn the face integrand for the stiffness matrix (BR-only part) for element id: "<<face->getPtrElementLeft()->getID();
 ElementT* right;
 ElementT* left=const_cast<ElementT*>(face->getPtrElementLeft());
 ElementInfos* leftInfo = static_cast<ElementInfos*> (left->getUserData());
 ElementInfos* rightInfo;
 PointElementReferenceT pLeft(3),pRight(3);
 double localepsilon;
 face->mapRefFaceToRefElemL(p,pLeft);
 //std::vector<LinearAlgebra::NumericalVector> leftValues,rightValues;
 //leftInfo->makeFunctionValuesVector(left,pLeft,leftValues);
 LinearAlgebra::NumericalVector normedNormal(3);
 normedNormal[0] = (normal*(1/Base::L2Norm(normal)))[0];
 normedNormal[1] = (normal*(1/Base::L2Norm(normal)))[1];
 normedNormal[2] = (normal*(1/Base::L2Norm(normal)))[2];
 int leftSize=left->getNrOfBasisFunctions();
 int dimension=left->getNrOfBasisFunctions();
 LinearAlgebra::NumericalVector phi_i(3),phi_j(3),dummy(3);
 if(face->isInternal()){
 //cout<<" and element id: "<<face->getPtrElementRight()->getID();
 right=const_cast<ElementT*>(face->getPtrElementRight());
 face->mapRefFaceToRefElemR(p,pRight);
 rightInfo = static_cast<ElementInfos*> (right->getUserData());
 dimension+=right->getNrOfBasisFunctions();
 //rightInfo->makeFunctionValuesVector(right,pRight,rightValues);
 }
 ret.resize(dimension,dimension);
 for(int i=0;i<dimension;++i){
 if(i<leftSize){
 //phi_i=leftValues[i];
 left->basisFunction(i, pLeft, phi_i);
 localepsilon=leftInfo->epsilon_;
 }else{
 //phi_i=rightValues[i-leftSize];
 right->basisFunction(i - leftSize, pRight, phi_i);
 localepsilon=rightInfo->epsilon_;
 }
 for(int j=0;j<dimension;++j){
 if(j<leftSize){
 //dummy=leftValues[j];
 left->basisFunction(j, pLeft, dummy);
 
 }else{
 //dummy=rightValues[j-leftSize];
 right->basisFunction(j-leftSize, pRight, dummy);
 dummy*=-1;
 }
 OuterProduct(normedNormal,dummy,phi_j);
 ret(j,i)=(face->isInternal()?1:2)*(phi_i[0]*phi_j[0]+phi_i[1]*phi_j[1]+phi_i[2]*phi_j[2])*sqrt(localepsilon);
 }
 }
 }
 
 void matrixFillerBR::faceIntegrand(const Base::Face* face, const LinearAlgebra::NumericalVector& normal, const Geometry::PointReference& p, LinearAlgebra::NumericalVector& ret){
 ElementT* left=const_cast<ElementT*>(face->getPtrElementLeft());
 ElementInfos* info = static_cast<ElementInfos*>(left->getUserData());
 ret.resize(left->getNrOfBasisFunctions());
 PointElementReferenceT PLeft(3);
 face->mapRefFaceToRefElemL(p,PLeft);
 //std::vector<LinearAlgebra::NumericalVector> functionValues;
 //info->makeFunctionValuesVector(left,PLeft,functionValues);
 PointPhysicalT PPhys(3);
 left->referenceToPhysical(PLeft,PPhys);
 LinearAlgebra::NumericalVector normedNormal(3);
 normedNormal[0] = (normal*(1/Base::L2Norm(normal)))[0];
 normedNormal[1] = (normal*(1/Base::L2Norm(normal)))[1];
 normedNormal[2] = (normal*(1/Base::L2Norm(normal)))[2];
 LinearAlgebra::NumericalVector val(3),dummy(3),phi(3);
 boundaryConditions(PPhys,dummy);
 OuterProduct(normedNormal,dummy,val);
 for(int i=0; i<left->getNrOfBasisFunctions(); ++i) {
 //phi=functionValues[i];
 left->basisFunction(i, PLeft, phi);
 ret(i)=2*(phi[0]*val[0]+phi[1]*val[1]+phi[2]*val[2]);
 }
 }
 
 
 void matrixFillerBR::faceIntegrand(const Base::Face* face, const LinearAlgebra::NumericalVector& normal, const Geometry::PointReference& p, LinearAlgebra::NumericalVector& ret)
 {
 
 int n = face->getPtrElementLeft()->getNrOfBasisFunctions();
 
 ElementT* left=const_cast<ElementT*>(face->getPtrElementLeft());
 PointElementReferenceT PLeft(3);
 face->mapRefFaceToRefElemL(p,PLeft);
 
 PointPhysicalT PPhys(3);
 left->referenceToPhysical(PLeft,PPhys);
 LinearAlgebra::NumericalVector normedNormal(3);
 
 normedNormal[0] = (normal*(1/Base::L2Norm(normal)))[0];
 normedNormal[1] = (normal*(1/Base::L2Norm(normal)))[1];
 normedNormal[2] = (normal*(1/Base::L2Norm(normal)))[2];
 
 LinearAlgebra::NumericalVector val(3),dummy(3),phi(3);
 boundaryConditions(PPhys,dummy);
 OuterProduct(normedNormal,dummy,val);
 //std::cout<<val<<std::endl;
 ret.resize(n);
 for(int i = 0; i < n; ++i)
 {
 face->basisFunction(i, p, phi);
 
 //std::cout<<phi<<std::endl;
 
 ret(i) = 2*(phi[0]*val[0]+phi[1]*val[1]+phi[2]*val[2]);
 }
 }
 
 
 void matrixFillerBR::fillMatrixes(hpGemUIExtentions* matrixContainer){
 time_t oldTime,newTime;
 time(&oldTime);
 std::cout<<"using a Brezzi-flux with penalty parameter "<<matrixContainer->getData()->StabCoeff_<<std::endl;
 int TotalAmountOfProcessors,localProcessorNumber,numberOfElements;
 MPI_Comm_size(PETSC_COMM_WORLD,&TotalAmountOfProcessors);
 MPI_Comm_rank(PETSC_COMM_WORLD,&localProcessorNumber);
 numberOfElements=matrixContainer->getNumberOfElements(0);
 
 //     if(matrixContainer->getNumberOfElements(0)%TotalAmountOfProcessors!=0){
 // 	std::cout<<"WARNING: the case where the number of elements is not a multiple of the number of processors has not been toroughly tested; gliches may occur!";
 //     }
 
 matrixContainer->ierr_=VecSetSizes(matrixContainer->x_,PETSC_DECIDE,matrixContainer->getData()->numberOfUnknowns_);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=VecSetSizes(matrixContainer->RHS_,PETSC_DECIDE,matrixContainer->getData()->numberOfUnknowns_);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=VecSetSizes(matrixContainer->derivative_,PETSC_DECIDE,matrixContainer->getData()->numberOfUnknowns_);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=MatSetSizes(matrixContainer->M_,PETSC_DECIDE,PETSC_DECIDE,matrixContainer->getData()->numberOfUnknowns_,matrixContainer->getData()->numberOfUnknowns_);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=MatSetSizes(matrixContainer->S_,PETSC_DECIDE,PETSC_DECIDE,matrixContainer->getData()->numberOfUnknowns_,matrixContainer->getData()->numberOfUnknowns_);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=VecSetBlockSize(matrixContainer->x_,matrixContainer->getConfigData()->numberOfBasisFunctions_);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=VecSetBlockSize(matrixContainer->RHS_,matrixContainer->getConfigData()->numberOfBasisFunctions_);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=VecSetBlockSize(matrixContainer->derivative_,matrixContainer->getConfigData()->numberOfBasisFunctions_);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=MatSetBlockSize(matrixContainer->M_,matrixContainer->getConfigData()->numberOfBasisFunctions_);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=MatSetBlockSize(matrixContainer->S_,matrixContainer->getConfigData()->numberOfBasisFunctions_);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=VecSetUp(matrixContainer->x_);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=VecSetUp(matrixContainer->RHS_);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=VecSetUp(matrixContainer->derivative_);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=MatSetUp(matrixContainer->M_);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=MatSetUp(matrixContainer->S_);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 LinearAlgebra::Matrix matrix(1,1);
 LinearAlgebra::NumericalVector vector(1);
 //hpGemUIExtentions::ElementFunction elF = &hpGemUIExtentions::elementMassIntegrand;
 PetscScalar *tempComplexArray=new PetscScalar[4*matrixContainer->getConfigData()->numberOfBasisFunctions_*matrixContainer->getConfigData()->numberOfBasisFunctions_];
 Base::ShortTermStorageElementBase* localElement_;
 localElement_ = new Base::ShortTermStorageElementHCurl(3); // creating object for H curl transformation
 Integration::ElementIntegral elIntegral(false);
 
 elIntegral.setStorageWrapper(localElement_);
 time(&newTime);
 std::cout<<"preparation took "<<difftime(newTime,oldTime)<<" seconds"<<std::endl;
 oldTime=newTime;
 for(hpGemUIExtentions::ElementIterator it=matrixContainer->elementColBegin(); it!=matrixContainer->elementColEnd(); ++it) {
 //	if(((*it)->getID())/(numberOfElements/TotalAmountOfProcessors+1)==localProcessorNumber) {
 matrix.resize((*it)->getNrOfBasisFunctions(),(*it)->getNrOfBasisFunctions());
 elIntegral.integrate<LinearAlgebra::Matrix>((*it),matrixContainer,matrix);
 if(matrixContainer->MHasToBeInverted_) {
 matrix.inverse(matrix);
 }
 int places[]= {(int)(*it)->getID()};
 for(int i=0; i<(*it)->getNrOfBasisFunctions()*(*it)->getNrOfBasisFunctions(); ++i) {
 tempComplexArray[i]=matrix[i];
 }
 matrixContainer->ierr_=MatSetValuesBlocked(matrixContainer->M_,1,places,1,places,&tempComplexArray[0],ADD_VALUES);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 //	}
 }
 time(&newTime);
 std::cout<<"filling the mass matrix took "<<difftime(newTime,oldTime)<<" seconds"<<std::endl;
 oldTime=newTime;
 //we are done filling M. PETSc may start assembling M in the background
 matrixContainer->ierr_=MatAssemblyBegin(matrixContainer->M_,MAT_FINAL_ASSEMBLY);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 //elF=&hpGemUIExtentions::elementStiffnessIntegrand;
 for(hpGemUIExtentions::ElementIterator it=matrixContainer->elementColBegin(); it!=matrixContainer->elementColEnd(); ++it) {
 //	if(((*it)->getID())/(numberOfElements/TotalAmountOfProcessors+1)==localProcessorNumber) {
 matrix.resize((*it)->getNrOfBasisFunctions(),(*it)->getNrOfBasisFunctions());
 elIntegral.integrate<LinearAlgebra::Matrix>((*it),this,matrix);
 int places[]= {(int)(*it)->getID()};
 for(int i=0; i<(*it)->getNrOfBasisFunctions()*(*it)->getNrOfBasisFunctions(); ++i) {
 tempComplexArray[i]=matrix[i];
 }
 matrixContainer->ierr_=MatSetValuesBlocked(matrixContainer->S_,1,places,1,places,&tempComplexArray[0],ADD_VALUES);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 //	}
 }
 //hpGemUIExtentions::FaceFunction faF = &hpGemUIExtentions::faceIntegrand;
 Base::ShortTermStorageFaceBase* localFace_;
 localFace_ = new Base::ShortTermStorageFaceHCurl(3); // creating object for H curl transformation
 
 Integration::FaceIntegral faIntegral(false);
 
 faIntegral.setStorageWrapper(localFace_);
 
 for(hpGemUIExtentions::FaceIterator it=matrixContainer->faceColBegin(); it!=matrixContainer->faceColEnd(); ++it) {
 //	if(((*it)->getPtrElementLeft()->getID())/(numberOfElements/TotalAmountOfProcessors+1)==localProcessorNumber) {
 if((*it)->isInternal()) {
 matrix.resize((*it)->getPtrElementLeft()->getNrOfBasisFunctions()+(*it)->getPtrElementRight()->getNrOfBasisFunctions(),(*it)->getPtrElementLeft()->getNrOfBasisFunctions()+(*it)->getPtrElementRight()->getNrOfBasisFunctions());
 int places[]= {(int)(*it)->getPtrElementLeft()->getID(),(int)(*it)->getPtrElementRight()->getID()};
 faIntegral.integrate<LinearAlgebra::Matrix>(*it,matrixContainer,matrix);
 std::cout<<matrix<<std::endl;
 for(int i=0; i<4*(*it)->getPtrElementLeft()->getNrOfBasisFunctions()*(*it)->getPtrElementLeft()->getNrOfBasisFunctions(); ++i) {
 tempComplexArray[i]=matrix[i];
 }
 matrixContainer->ierr_=MatSetValuesBlocked(matrixContainer->S_,2,places,2,places,&tempComplexArray[0],ADD_VALUES);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 } else {
 matrix.resize((*it)->getPtrElementLeft()->getNrOfBasisFunctions(),(*it)->getPtrElementLeft()->getNrOfBasisFunctions());
 int places[]= {(int)(*it)->getPtrElementLeft()->getID()};
 faIntegral.integrate<LinearAlgebra::Matrix>(*it,matrixContainer,matrix);
 std::cout<<matrix<<std::endl;
 for(int i=0; i<(*it)->getPtrElementLeft()->getNrOfBasisFunctions()*(*it)->getPtrElementLeft()->getNrOfBasisFunctions(); ++i) {
 tempComplexArray[i]=matrix[i];
 }
 matrixContainer->ierr_=MatSetValuesBlocked(matrixContainer->S_,1,places,1,places,&tempComplexArray[0],ADD_VALUES);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 }
 //	}
 }
 //we will be needing a few entries of M so assembly should be done by now
 matrixContainer->ierr_=MatAssemblyEnd(matrixContainer->M_,MAT_FINAL_ASSEMBLY);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 
 //make some assistant matrices
 Mat DInternal,DBoundary,MLocal,dummy,dummy2;
 matrixContainer->ierr_=MatCreateDense(PETSC_COMM_SELF,PETSC_DECIDE,PETSC_DECIDE,2*matrixContainer->getConfigData()->numberOfBasisFunctions_,2*matrixContainer->getConfigData()->numberOfBasisFunctions_,NULL,&DInternal);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=MatCreateDense(PETSC_COMM_SELF,PETSC_DECIDE,PETSC_DECIDE,matrixContainer->getConfigData()->numberOfBasisFunctions_,matrixContainer->getConfigData()->numberOfBasisFunctions_,NULL,&DBoundary);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=MatSetUp(DInternal);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=MatSetUp(DBoundary);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 int localnElements;
 int dimension,skip(0);
 //faF = &faceIntegrandBRPart;
 for(hpGemUIExtentions::FaceIterator it=matrixContainer->faceColBegin(); it!=matrixContainer->faceColEnd(); ++it) {
 int placesBlocked[]= {(int)(*it)->getPtrElementLeft()->getID(),-1};
 if((*it)->isInternal()) {
 localnElements=2;
 dimension=(*it)->getPtrElementLeft()->getNrOfBasisFunctions()+(*it)->getPtrElementRight()->getNrOfBasisFunctions();
 placesBlocked[1]=(*it)->getPtrElementRight()->getID();
 } else {
 localnElements=1;
 dimension=(*it)->getPtrElementLeft()->getNrOfBasisFunctions();
 }
 
 IS ISplaces;
 matrixContainer->ierr_=ISCreateBlock(PETSC_COMM_WORLD,(*it)->getPtrElementLeft()->getNrOfBasisFunctions(),localnElements,placesBlocked,PETSC_COPY_VALUES,&ISplaces);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=MatGetSubMatrix(matrixContainer->M_,ISplaces,ISplaces,MAT_INITIAL_MATRIX,&MLocal);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 //	if(((*it)->getPtrElementLeft()->getID())/(numberOfElements/TotalAmountOfProcessors+1)==localProcessorNumber) {
 matrix.resize(dimension,dimension);
 faIntegral.integrate<LinearAlgebra::Matrix>(*it,this,matrix);
 std::cout<<matrix<<std::endl;
 for(int i=0; i<dimension*dimension; ++i) {
 tempComplexArray[i]=matrix[i];
 }
 
 //make sure to throw away the bogus parts again
 matrixContainer->ierr_=MatGetDiagonalBlock(MLocal,&dummy2);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 
 //this is done locally so make sure to use sequential data structures
 matrixContainer->ierr_=MatConvert(dummy2,"seqdense",MAT_REUSE_MATRIX,&dummy2);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 
 //the Brezzi flux really need M^-1, so invert M anyway, even if the inverse is not explicitly needed later on
 if(!matrixContainer->MHasToBeInverted_) {
 LinearAlgebra::Matrix localM(dimension,dimension);
 PetscScalar* data;
 matrixContainer->ierr_=MatDenseGetArray(dummy2,&data);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 for(int i=0; i<dimension*dimension; ++i) {
 localM[i]=data[i].real();//insert the real part of the data in the underlying datastructions of localM
 }
 localM.inverse(localM);
 for(int i=0; i<dimension*dimension; ++i) {
 data[i]=localM[i];
 }
 matrixContainer->ierr_=MatDenseRestoreArray(dummy2,&data);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 }
 int places[dimension];
 for(int i=0; i<dimension; ++i) {
 places[i]=i;
 }
 matrixContainer->ierr_=MatConvert(dummy2,"seqaij",MAT_REUSE_MATRIX,&dummy2);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 if((*it)->isInternal()) {
 matrixContainer->ierr_=MatSetValues(DInternal,dimension,places,dimension,places,&tempComplexArray[0],INSERT_VALUES);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=MatAssemblyBegin(DInternal,MAT_FINAL_ASSEMBLY);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=MatAssemblyEnd(DInternal,MAT_FINAL_ASSEMBLY);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=MatConvert(DInternal,"seqaij",MAT_REUSE_MATRIX,&DInternal);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 
 //computation of D^t M D
 matrixContainer->ierr_=MatTransposeMatMult(DInternal,dummy2,MAT_INITIAL_MATRIX,1.0,&dummy);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=MatDestroy(&dummy2);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=MatMatMult(dummy,DInternal,MAT_INITIAL_MATRIX,1.0,&dummy2);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 } else {
 matrixContainer->ierr_=MatSetValues(DBoundary,dimension,places,dimension,places,&tempComplexArray[0],INSERT_VALUES);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=MatAssemblyBegin(DBoundary,MAT_FINAL_ASSEMBLY);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=MatAssemblyEnd(DBoundary,MAT_FINAL_ASSEMBLY);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=MatConvert(DBoundary,"seqaij",MAT_REUSE_MATRIX,&DBoundary);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=MatTransposeMatMult(DBoundary,dummy2,MAT_INITIAL_MATRIX,1.0,&dummy);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=MatDestroy(&dummy2);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=MatMatMult(dummy,DBoundary,MAT_INITIAL_MATRIX,1.0,&dummy2);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 
 //not that we have D^T*M^-1 use it to also compute boundary terms
 Vec DRHS,dummy3;
 matrixContainer->ierr_=VecCreate(PETSC_COMM_SELF,&DRHS);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=VecSetSizes(DRHS,PETSC_DECIDE,dimension);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=VecSetUp(DRHS);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=VecDuplicate(DRHS,&dummy3);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 //faF = &hpGemUIExtentions::sourceTermBoundaryBR;
 vector.resize(dimension);
 faIntegral.integrate<LinearAlgebra::NumericalVector>(*it,this,vector);
 std::cout<<vector<<std::endl;
 for(int i=0; i<dimension; ++i) {
 tempComplexArray[i]=vector[i];
 }
 matrixContainer->ierr_=VecSetValues(DRHS,dimension,places,&tempComplexArray[0],ADD_VALUES);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=VecAssemblyBegin(DRHS);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=VecAssemblyEnd(DRHS);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=MatMult(dummy,DRHS,dummy3);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=VecScale(dummy3,matrixContainer->getData()->StabCoeff_);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 const PetscScalar* arrayLoc;
 matrixContainer->ierr_=VecGetArrayRead(dummy3,&arrayLoc);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=VecSetValuesBlocked(matrixContainer->RHS_,1,placesBlocked,arrayLoc,ADD_VALUES);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=VecRestoreArrayRead(dummy3,&arrayLoc);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=VecDestroy(&DRHS);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=VecDestroy(&dummy3);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 //faF = &hpGemUIExtentions::faceIntegrandBRPart;
 }
 matrixContainer->ierr_=MatConvert(dummy2,"seqdense",MAT_REUSE_MATRIX,&dummy2);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 PetscScalar* arrayLoc;
 matrixContainer->ierr_=MatScale(dummy2,matrixContainer->getData()->StabCoeff_);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=MatDenseGetArray(dummy2,&arrayLoc);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=MatSetValuesBlocked(matrixContainer->S_,localnElements,placesBlocked,localnElements,placesBlocked,arrayLoc,ADD_VALUES);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=MatDenseRestoreArray(dummy2,&arrayLoc);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=MatDestroy(&dummy);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=MatDestroy(&dummy2);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=ISDestroy(&ISplaces);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 //	}
 }
 time(&newTime);
 std::cout<<"filling the stiffness matrix took "<<difftime(newTime,oldTime)<<" seconds"<<std::endl;
 oldTime=newTime;
 //elF=&hpGemUIExtentions::sourceTerm;
 for(hpGemUIExtentions::ElementIterator it=matrixContainer->elementColBegin(); it!=matrixContainer->elementColEnd(); ++it) {
 //	if(((*it)->getID())/(numberOfElements/TotalAmountOfProcessors+1)==localProcessorNumber) {
 vector.resize((*it)->getNrOfBasisFunctions());
 elIntegral.integrate<LinearAlgebra::NumericalVector>((*it),this,vector);
 int places[]= {(int)(*it)->getID()};
 for(int i=0; i<(*it)->getNrOfBasisFunctions(); ++i) {
 tempComplexArray[i]=vector[i];
 }
 matrixContainer->ierr_=VecSetValuesBlocked(matrixContainer->RHS_,1,places,&tempComplexArray[0],ADD_VALUES);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 //	}
 }
 matrixContainer->ierr_=MatAssemblyBegin(matrixContainer->S_,MAT_FINAL_ASSEMBLY);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 //elF=&hpGemUIExtentions::initialConditions;
 for(hpGemUIExtentions::ElementIterator it=matrixContainer->elementColBegin(); it!=matrixContainer->elementColEnd(); ++it) {
 //	if(((*it)->getID())/(numberOfElements/TotalAmountOfProcessors+1)==localProcessorNumber) {
 vector.resize((*it)->getNrOfBasisFunctions());
 elIntegral.integrate<LinearAlgebra::NumericalVector>((*it),matrixContainer,vector);
 int places[]= {(int)(*it)->getID()};
 for(int i=0; i<(*it)->getNrOfBasisFunctions(); ++i) {
 tempComplexArray[i]=vector[i];
 }
 matrixContainer->ierr_=VecSetValuesBlocked(matrixContainer->x_,1,places,&tempComplexArray[0],ADD_VALUES);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 //	}
 }
 matrixContainer->ierr_=VecAssemblyBegin(matrixContainer->x_);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 //elF=&hpGemUIExtentions::initialConditionsDeriv;
 for(hpGemUIExtentions::ElementIterator it=matrixContainer->elementColBegin(); it!=matrixContainer->elementColEnd(); ++it) {
 //	if(((*it)->getID())/(numberOfElements/TotalAmountOfProcessors+1)==localProcessorNumber) {
 vector.resize((*it)->getNrOfBasisFunctions());
 elIntegral.integrate((*it),elF,vector,matrixContainer);
 int places[]= {(int)(*it)->getID()};
 for(int i=0; i<(*it)->getNrOfBasisFunctions(); ++i) {
 tempComplexArray[i]=matrix[i];
 }
 matrixContainer->ierr_=VecSetValuesBlocked(matrixContainer->derivative_,1,places,&tempComplexArray[0],ADD_VALUES);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 //	}
 }
 matrixContainer->ierr_=VecAssemblyBegin(matrixContainer->derivative_);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 //faF = &hpGemUIExtentions::sourceTermBoundary;
 for(hpGemUIExtentions::FaceIterator it=matrixContainer->faceColBegin(); it!=matrixContainer->faceColEnd(); ++it){
 //	if(((*it)->getPtrElementLeft()->getID())/(numberOfElements/TotalAmountOfProcessors+1)==localProcessorNumber) {
 if((*it)->isInternal()) {
 //internal faces have no boundary contributions
 } else {
 vector.resize((*it)->getPtrElementLeft()->getNrOfBasisFunctions());
 int places[]= {(int)(*it)->getPtrElementLeft()->getID()};
 faIntegral.integrate<LinearAlgebra::NumericalVector>(*it,matrixContainer,vector);
 for(int i=0; i<(*it)->getPtrElementLeft()->getNrOfBasisFunctions(); ++i) {
 tempComplexArray[i]=vector[i];
 }
 matrixContainer->ierr_=VecSetValuesBlocked(matrixContainer->RHS_,1,places,&tempComplexArray[0],ADD_VALUES);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 }
 //	}
 }
 time(&newTime);
 std::cout<<"filling vectors took "<<difftime(newTime,oldTime)<<" seconds"<<std::endl;
 oldTime=newTime;
 matrixContainer->ierr_=VecAssemblyBegin(matrixContainer->RHS_);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=VecAssemblyEnd(matrixContainer->x_);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=VecAssemblyEnd(matrixContainer->derivative_);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=MatAssemblyEnd(matrixContainer->S_,MAT_FINAL_ASSEMBLY);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=VecAssemblyEnd(matrixContainer->RHS_);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 
 matrixContainer->ierr_=MatDestroy(&DInternal);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 matrixContainer->ierr_=MatDestroy(&DBoundary);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
 time(&newTime);
 std::cout<<"communication took "<<difftime(newTime,oldTime)<<" extra seconds"<<std::endl;
 delete[] tempComplexArray;
 }
 */
