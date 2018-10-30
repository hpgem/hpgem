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

#include "BaseExtended.h"
#include "petscsys.h"
#include "petscksp.h"

#include <cmath>

hpGemUIExtentions::hpGemUIExtentions(Base::GlobalData * const globalConfig, Base::ConfigurationData* elementConfig)
        : HpgemAPIBase(globalConfig, elementConfig)
{
    PetscErrorCode ierr_;
#if PETSC_VERSION_GE(3, 7, 0)
    ierr_ = PetscLogDefaultBegin();
#else
    ierr_ = PetscLogBegin();
#endif
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);

}

hpGemUIExtentions::~hpGemUIExtentions()
{

    PetscErrorCode ierr_;
    PetscViewer log;
    ierr_ = PetscViewerASCIIOpen(PETSC_COMM_WORLD, "maxwell.log", &log);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = PetscLogView(log);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);

    ierr_ = PetscViewerDestroy(&log);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);

}

std::size_t hpGemUIExtentions::addMesh(Base::MeshManipulator<DIM>* mesh)
{
    meshes_.push_back(mesh);
    return meshes_.size() - 1;
}


Base::MeshManipulator<DIM> * hpGemUIExtentions::getMesh(std::size_t meshId)
{
    return meshes_[meshId];
}

const Base::ConfigurationData* hpGemUIExtentions::getConfigData()
{
    return configData_;
}

/* See solveDOS() in baseExtended.h for motivation of commenting out.
void hpGemUIExtentions::LDOSIntegrand(Base::PhysicalElement<DIM>& element, double &ret)
{ //currently LDOS is computed by evaluation at a point so integrand is a bit of a misnomer
    //ElementInfos* info = static_cast<ElementInfos*>(const_cast<ElementT*>(element)->getUserData());
    LinearAlgebra::SmallVector<DIM> Phi, PhiRealI, PhiRealJ, PhiImagI, PhiImagJ;
    //std::vector<LinearAlgebra::NumericalVector> functionValues(element->getNrOfBasisFunctions());
    //info->makeFunctionValuesVector(element,p,functionValues);
    for (int i = 0; i < element->getNrOfBasisFunctions(); ++i)
    {
        
        //PhiRealI+=functionValues[i]*element->getData(0,0,i);
        //PhiImagI+=functionValues[i]*element->getData(1,0,i);
        //PhiRealJ+=functionValues[i]*element->getData(2,0,i);
        //PhiImagJ+=functionValues[i]*element->getData(3,0,i);
        
        element->basisFunction(i, p, Phi);
        
        PhiRealI += Phi * element->getData(0, 0, i);
        PhiImagI += Phi * element->getData(1, 0, i);
        PhiRealJ += Phi * element->getData(2, 0, i);
        PhiImagJ += Phi * element->getData(3, 0, i);
    }
    //current implementation: orientation averaged LDOS; modify the next expression for other orientations
    ret = PhiRealI[0] * PhiRealJ[0] + PhiRealI[1] * PhiRealJ[1] + PhiRealI[2] * PhiRealJ[2] + PhiImagI[0] * PhiImagJ[0] + PhiImagI[1] * PhiImagJ[1] + PhiImagI[2] * PhiImagJ[2];
    ret *= 48; //assume computation is done only in the irreducible brillouin zone
    //cout<<ret<<std::endl;
}

void hpGemUIExtentions::makeFunctionValue(Vec eigenVector, LinearAlgebra::MiddleSizeVector& result)
{
    Vec scaledVec;
    double partialResult;
    int index(0);
    ierr_ = VecDuplicate(eigenVector, &scaledVec);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = MatMult(M_, eigenVector, scaledVec);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    result.resize(4 * meshes_[0]->getNumberOfElements());
    for (ElementIterator it = meshes_[0]->elementColBegin(); it != meshes_[0]->elementColEnd(); ++it)
    {
        ierr_ = VecGetArrayRead(eigenVector, storage_);
        CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        for (int i = 0; i < (*it)->getNrOfBasisFunctions(); ++i)
        {
            (*it)->setData(0, 0, i, (*storage_)[(*it)->getNrOfBasisFunctions() * ((*it)->getID()) + i].real());
            (*it)->setData(1, 0, i, (*storage_)[(*it)->getNrOfBasisFunctions() * ((*it)->getID()) + i].imag());
        }
        ierr_ = VecRestoreArrayRead(eigenVector, storage_);
        CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        ierr_ = VecGetArrayRead(scaledVec, storage_);
        CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        for (int i = 0; i < (*it)->getNrOfBasisFunctions(); ++i)
        {
            (*it)->setData(2, 0, i, (*storage_)[(*it)->getNrOfBasisFunctions() * ((*it)->getID()) + i].real());
            (*it)->setData(3, 0, i, (*storage_)[(*it)->getNrOfBasisFunctions() * ((*it)->getID()) + i].imag());
        }
        ierr_ = VecRestoreArrayRead(scaledVec, storage_);
        CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        for (int i = 0; i < 4; ++i)
        {
            //PointElementReferenceT p(3);
            const PointElementReferenceT& p = (*it)->getReferenceGeometry()->getCenter();
            LDOSIntegrand(*it, partialResult);
            result[index] = partialResult;
            ++index;
        }
    }
    ierr_ = VecDestroy(&scaledVec);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
}
*/

// See baseExtended.h
//void hpGemUIExtentions::exportMatrixes()
//{
//    std::cout << "genereting Matlab scripts to load the matrixes" << std::endl;
//    MHasToBeInverted_ = false;
//    //assembler->fillMatrices(this);
//    PetscViewer viewM, viewS;
//    PetscViewerASCIIOpen(MPI_COMM_WORLD, "M.m", &viewM);
//    PetscViewerASCIIOpen(MPI_COMM_WORLD, "S.m", &viewS);
//    PetscViewerSetFormat(viewM, PETSC_VIEWER_ASCII_MATLAB);
//    PetscViewerSetFormat(viewS, PETSC_VIEWER_ASCII_MATLAB);
//    MatView(M_, viewM);
//    MatView(S_, viewS);
//    PetscViewerDestroy(&viewM);
//    PetscViewerDestroy(&viewS);
//}


/*
void hpGemUIExtentions::solveDOS()
{
    std::cout << "finding the local density of states" << std::endl;
    const MaxwellData* actualdata = getData();
    int degreesOfFreedomPerElement = configData_->numberOfBasisFunctions_;
    std::valarray<PetscScalar> blockvalues(degreesOfFreedomPerElement * degreesOfFreedomPerElement);
    int measureAmount = 0;
    std::vector<double> eigenvalues;
    std::vector<LinearAlgebra::MiddleSizeVector> functionValues;
    std::vector<IS> xboundaryRow, xboundaryCol, yboundaryRow, yboundaryCol, zboundaryRow, zboundaryCol;
    findBoundaryBlocks(xboundaryRow, xboundaryCol, yboundaryRow, yboundaryCol, zboundaryRow, zboundaryCol);
    KspaceData brillouinZone(5);
    
    //For IP-DG solving a general eigenproblem is slightly faster,
    //but the Brezzi formulation needs an inverse of each block in the mass matrix anyway
    //so there the general eigenproblem is just more work
    MHasToBeInverted_ = true;
    assembler->fillMatrices(this);
    
    Utilities::GlobalPetscMatrix M_(meshes_[0], 0, -1), S_(meshes_[0], 1, 0);
    std::cout << "GlobalPetscMatrix initialised" << std::endl;
    Utilities::GlobalPetscVector x_(meshes_[0], 0, -1), derivative_(meshes_[0], 1, -1), RHS_(meshes_[0], 2, 0);
    std::cout << "GlobalPetscVector initialised" << std::endl;
    x_.assemble();
    std::cout << "x_ assembled" << std::endl;
    RHS_.assemble();
    std::cout << "RHS_ assembled" << std::endl;
    derivative_.assemble();
    std::cout << "derivative_ assembled" << std::endl;
    
    //Divide the eigenvalues by pi^2 (because some people cant do that by heart and the useful info shouldn't be hidden)
    //ierr_=MatScale(M_,1/M_PI/M_PI);CHKERRABORT(PETSC_COMM_WORLD,ierr_);
    
    //The eigensolver doesn't like block structures
    //ierr_=MatConvert(M_,"aij",MAT_REUSE_MATRIX,&M_);CHKERRABORT(PETSC_COMM_WORLD,ierr_);
    //ierr_=MatConvert(S_,"aij",MAT_REUSE_MATRIX,&S_);CHKERRABORT(PETSC_COMM_WORLD,ierr_);
    Mat product, dummy;
    ierr_ = MatMatMult(M_, S_, MAT_INITIAL_MATRIX, 1.0, &product);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    
    //    ierr_=EPSSetProblemType(eigenSolver_,EPS_HEP);
    ierr_ = EPSSetOperators(eigenSolver_, product, NULL);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = EPSSetUp(eigenSolver_);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = EPSSetDimensions(eigenSolver_, 24, PETSC_DECIDE, PETSC_DECIDE);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    // 	ierr_=EPSSetType(eigenSolver_,"jd");CHKERRABORT(PETSC_COMM_WORLD,ierr_);
    // 	ST transformation;
    // 	KSP spectralProblem;
    // 	PC preconditioner;
    // 	ierr_=EPSGetST(eigenSolver_,&transformation);CHKERRABORT(PETSC_COMM_WORLD,ierr_);
    // 	ierr_=STSetType(transformation,"precond");CHKERRABORT(PETSC_COMM_WORLD,ierr_);
    // 	ierr_=STGetKSP(transformation,&spectralProblem);CHKERRABORT(PETSC_COMM_WORLD,ierr_);
    // 	ierr_=KSPGetPC(spectralProblem,&preconditioner);CHKERRABORT(PETSC_COMM_WORLD,ierr_);
    // 	ierr_=PCSetType(preconditioner,"ilu");CHKERRABORT(PETSC_COMM_WORLD,ierr_);
    // 	ierr_=EPSSetST(eigenSolver_,transformation);CHKERRABORT(PETSC_COMM_WORLD,ierr_);
    
    //everything that is set in the code, but before this line is overridden by command-line options
    ierr_ = EPSSetFromOptions(eigenSolver_);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    //everything that is set in the code, but after this line overrides the comand-line options
    
    PetscScalar target;
    ierr_ = EPSGetTarget(eigenSolver_, &target);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    
    ierr_ = EPSSolve(eigenSolver_);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    
    PetscScalar neededOnlyForRealPetsc, eigenvalue;
    Vec example, *eigenVectors;
    eigenVectors = new Vec[40]; //a few extra in case SLEPc finds more than the requested amount of eigenvalues
    ierr_ = VecCreate(PETSC_COMM_WORLD, &example);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecSetSizes(example, PETSC_DECIDE, 61);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecSetUp(example);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecDuplicateVecs(x_, 40, &eigenVectors);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecDestroy(&example);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    int converged;
    ierr_ = EPSGetConverged(eigenSolver_, &converged);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    
    ierr_ = EPSGetInvariantSubspace(eigenSolver_, eigenVectors);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    for (int i = 0; i < converged && i < 24; ++i)
    {
        ierr_ = EPSGetEigenvalue(eigenSolver_, i, &eigenvalue, &neededOnlyForRealPetsc);
        CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        //SLEPSc sorts the eigenvalues by closest to the target first, this reorders them for the better smallest first
        //if SLEPSc manages to find some of the zero eigenvalues they should be skipped
        if (eigenvalue.real() > 1e-6)
        {
            LinearAlgebra::MiddleSizeVector functionvalue;
            //makeFunctionValue(eigenVectors[i],functionvalue);
            if (eigenvalue.real() < target.real())
            {
                eigenvalues.insert(eigenvalues.begin(), sqrt(eigenvalue.real()));
                functionValues.insert(functionValues.begin(), functionvalue);
            }
            else
            {
                eigenvalues.push_back(sqrt(eigenvalue.real()));
                functionValues.push_back(functionvalue);
            }
        }
    }
    //eigenvalues.insert(eigenvalues.begin(),2,0);//at the infinite wavelength limit the two constant functions also contribute to the DOS
    brillouinZone.setOmega(eigenvalues);
    //brillouinZone.setFunctionValues(functionValues);
    
    Vec waveVec, waveVecConjugate;
    const int *rows, *columns;
    ierr_ = VecCreate(PETSC_COMM_WORLD, &waveVec);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecSetType(waveVec, "mpi");
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecSetBlockSize(waveVec, configData_->numberOfBasisFunctions_);
    ierr_ = VecSetSizes(waveVec, PETSC_DECIDE, globalData_->numberOfUnknowns_);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecSetUp(waveVec);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    LinearAlgebra::NumericalVector k(3);
    k[0] = M_PI / 60;
    k[1] = 0;
    k[2] = 0;
    // See the Eigenvalue problem.
    makeShiftMatrix(k, waveVec);
    ierr_ = VecAssemblyBegin(waveVec);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecAssemblyEnd(waveVec);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecDuplicate(waveVec, &waveVecConjugate);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecCopy(waveVec, waveVecConjugate);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecConjugate(waveVecConjugate);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    
    PetscLayout distribution;
    int localProcessorNumber, blockProcessorNumber, localsize;
    MPI_Comm_rank(PETSC_COMM_WORLD, &localProcessorNumber);
    ierr_ = PetscLayoutCreate(PETSC_COMM_WORLD, &distribution);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = MatGetLocalSize(product, &localsize, NULL);
    ierr_ = PetscLayoutSetLocalSize(distribution, localsize);
    ierr_ = PetscLayoutSetUp(distribution);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    
    while (brillouinZone.hasNextPoint())
    {
        k = brillouinZone.nextPoint();
        ierr_ = MatDiagonalScale(product, waveVec, waveVecConjugate);
        CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        
        //find a better strategy -- at least there is no communication in this set-up so MatAssemblyBegin should be reasonably fast
        //a cleaner set-up can be made using two iterators
        for (int j = 0; j < xboundaryCol.size(); ++j)
        {
            ierr_ = ISGetIndices(xboundaryRow[j], &rows);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            ierr_ = ISGetIndices(xboundaryCol[j], &columns);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            PetscLayoutFindOwner(distribution, rows[0], &blockProcessorNumber);
            if (blockProcessorNumber == localProcessorNumber)
            {
                ierr_ = MatGetValues(product, degreesOfFreedomPerElement, rows, degreesOfFreedomPerElement, columns, &blockvalues[0]);
                CHKERRABORT(PETSC_COMM_WORLD, ierr_);
                blockvalues *= exp(std::complex<double>(0, k[0] * (j % 2 == 0 ? -1 : 1)));
                ierr_ = MatSetValues(product, degreesOfFreedomPerElement, rows, degreesOfFreedomPerElement, columns, &blockvalues[0], INSERT_VALUES);
                CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            }
            ierr_ = MatAssemblyBegin(product, MAT_FINAL_ASSEMBLY);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            ierr_ = MatAssemblyEnd(product, MAT_FINAL_ASSEMBLY);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            ierr_ = ISRestoreIndices(xboundaryRow[j], &rows);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            ierr_ = ISRestoreIndices(xboundaryCol[j], &columns);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        }
        for (int j = 0; j < yboundaryCol.size(); ++j)
        {
            ierr_ = ISGetIndices(yboundaryRow[j], &rows);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            ierr_ = ISGetIndices(yboundaryCol[j], &columns);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            PetscLayoutFindOwner(distribution, rows[0], &blockProcessorNumber);
            if (blockProcessorNumber == localProcessorNumber)
            {
                ierr_ = MatGetValues(product, degreesOfFreedomPerElement, rows, degreesOfFreedomPerElement, columns, &blockvalues[0]);
                CHKERRABORT(PETSC_COMM_WORLD, ierr_);
                blockvalues *= exp(std::complex<double>(0, k[1] * (j % 2 == 0 ? -1 : 1)));
                ierr_ = MatSetValues(product, degreesOfFreedomPerElement, rows, degreesOfFreedomPerElement, columns, &blockvalues[0], INSERT_VALUES);
                CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            }
            ierr_ = MatAssemblyBegin(product, MAT_FINAL_ASSEMBLY);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            ierr_ = MatAssemblyEnd(product, MAT_FINAL_ASSEMBLY);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            ierr_ = ISRestoreIndices(yboundaryRow[j], &rows);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            ierr_ = ISRestoreIndices(yboundaryCol[j], &columns);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        }
        for (int j = 0; j < zboundaryCol.size(); ++j)
        {
            ierr_ = ISGetIndices(zboundaryRow[j], &rows);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            ierr_ = ISGetIndices(zboundaryCol[j], &columns);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            PetscLayoutFindOwner(distribution, rows[0], &blockProcessorNumber);
            if (blockProcessorNumber == localProcessorNumber)
            {
                ierr_ = MatGetValues(product, degreesOfFreedomPerElement, rows, degreesOfFreedomPerElement, columns, &blockvalues[0]);
                CHKERRABORT(PETSC_COMM_WORLD, ierr_);
                blockvalues *= exp(std::complex<double>(0, k[2] * (j % 2 == 0 ? -1 : 1)));
                ierr_ = MatSetValues(product, degreesOfFreedomPerElement, rows, degreesOfFreedomPerElement, columns, &blockvalues[0], INSERT_VALUES);
                CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            }
            ierr_ = MatAssemblyBegin(product, MAT_FINAL_ASSEMBLY);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            ierr_ = MatAssemblyEnd(product, MAT_FINAL_ASSEMBLY);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            ierr_ = ISRestoreIndices(zboundaryRow[j], &rows);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            ierr_ = ISRestoreIndices(zboundaryCol[j], &columns);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        }
        
        ierr_ = EPSSetOperators(eigenSolver_, product, NULL);
        CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        ierr_ = EPSSetInitialSpace(eigenSolver_, converged, eigenVectors);
        CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        ierr_ = EPSSetUp(eigenSolver_);
        CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        ierr_ = EPSSolve(eigenSolver_);
        CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        ierr_ = EPSGetInvariantSubspace(eigenSolver_, eigenVectors);
        CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        ierr_ = EPSGetConverged(eigenSolver_, &converged);
        CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        
        eigenvalues.clear();
        functionValues.clear();
        for (int j = 0; j < converged && j < 24; ++j)
        {
            ierr_ = EPSGetEigenvalue(eigenSolver_, j, &eigenvalue, &neededOnlyForRealPetsc);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            if (eigenvalue.real() > 1e-6)
            {
                LinearAlgebra::NumericalVector functionvalue;
                //makeFunctionValue(eigenVectors[j],functionvalue);
                if (eigenvalue.real() < target.real())
                {
                    eigenvalues.insert(eigenvalues.begin(), sqrt(eigenvalue.real()));
                    functionValues.insert(functionValues.begin(), functionvalue);
                }
                else
                {
                    eigenvalues.push_back(sqrt(eigenvalue.real()));
                    functionValues.push_back(functionvalue);
                }
            }
        }
        brillouinZone.setOmega(eigenvalues);
        //brillouinZone.setFunctionValues(functionValues);
    }
    for (int i = 0; i < xboundaryCol.size(); ++i)
    {
        ISDestroy(&xboundaryCol[i]);
        ISDestroy(&xboundaryRow[i]);
    }
    for (int i = 0; i < yboundaryCol.size(); ++i)
    {
        ISDestroy(&yboundaryCol[i]);
        ISDestroy(&yboundaryRow[i]);
    }
    for (int i = 0; i < zboundaryCol.size(); ++i)
    {
        ISDestroy(&zboundaryCol[i]);
        ISDestroy(&zboundaryRow[i]);
    }
    
    ierr_ = VecDestroy(&waveVec);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecDestroy(&waveVecConjugate);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    
    //always clean up after you are done
    ierr_ = MatDestroy(&product);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    
    LinearAlgebra::NumericalVector result(meshes_[0]->getNumberOfElements());
    for (int i = 0; i < 1001; ++i)
    {
        brillouinZone.getIntegral(double(i) / 100., result);
        std::cout << result[0] << " " << result[6] << " " << result[10] << std::endl;
        result[0] = 0;
        //result[6]=0;
        //result[10]=0;
    }
}
*/
