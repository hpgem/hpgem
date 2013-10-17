#include "fillMatrices.hpp"
#include "BaseExtended.hpp"

void matrixFillerIP::fillMatrixes(hpGemUIExtentions< 3 >* matrixContainer){
        time_t oldTime,newTime;
	time(&oldTime);
        cout<<"using an IP-DG method with penalty parameter "<<matrixContainer->getData()->StabCoeff_<<endl;
        //learn your own processor number and the total amount of processors so there is no double work done while filling the matrices
        int TotalAmountOfProcessors,localProcessorNumber,numberOfElements;
        MPI_Comm_size(PETSC_COMM_WORLD,&TotalAmountOfProcessors);
        MPI_Comm_rank(PETSC_COMM_WORLD,&localProcessorNumber);
        numberOfElements=matrixContainer->getNumberOfElements(0);
	
// 	if(matrixContainer->getNumberOfElements(0)%TotalAmountOfProcessors!=0){
// 	    cout<<"WARNING: the case where the number of elements is not a multiple of the number of processors has not been toroughly tested; gliches may occur!";
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

        //TODO make this consist of a bit less code repetitions...
        //TODO check if it is faster to do all the integrations in one element loop and one face loop
        //(but requires a lot of changing of integrands, also connot communicate data that is not yet needed)
        LinearAlgebra::Matrix matrix(1,1);
        PetscScalar tempComplexArray[4*matrixContainer->getConfigData()->numberOfBasisFunctions_*matrixContainer->getConfigData()->numberOfBasisFunctions_];
        hpGemUIExtentions<3>::ElementFunction elF = &hpGemUIExtentions<3>::elementMassIntegrand;
        Integration::ElementIntegral<3> elIntegral(false);
	time(&newTime);
	cout<<"filling the matrices; preparation took "<<difftime(newTime,oldTime)<<" seconds"<<endl;
	oldTime=newTime;
        for(hpGemUIExtentions<3>::ElementIterator it=matrixContainer->elementColBegin(); it!=matrixContainer->elementColEnd(); ++it) {
            if(((*it)->getID()-1)/(numberOfElements/TotalAmountOfProcessors+1)==localProcessorNumber) {
                matrix.resize((*it)->getNrOfUnknows(),(*it)->getNrOfUnknows());
                elIntegral.integrate((*it),elF,matrix,matrixContainer);
                if(matrixContainer->MHasToBeInverted_) {
                    matrix.inverse(matrix);
                }
                int places[]= {(*it)->getID()-1};
                for(int i=0; i<(*it)->getNrOfUnknows()*(*it)->getNrOfUnknows(); ++i) {
                    tempComplexArray[i]=matrix[i];//TODO Make complex PETSc read reals
                }
                matrixContainer->ierr_=MatSetValuesBlocked(matrixContainer->M_,1,places,1,places,&tempComplexArray[0],INSERT_VALUES);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
            }//the & and the [0] are still here for easy conversion back to &matrix[0]
        }
        //PETSc can start assembling M if it wants to; we wont set any extra antries anymore
        matrixContainer->ierr_=MatAssemblyBegin(matrixContainer->M_,MAT_FINAL_ASSEMBLY);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);  
	time(&newTime);
	cout<<"filling the mass matrix took "<<difftime(newTime,oldTime)<<" seconds"<<endl;
	oldTime=newTime;                                              
        elF=&hpGemUIExtentions<3>::elementStiffnessIntegrand;
        for(hpGemUIExtentions<3>::ElementIterator it=matrixContainer->elementColBegin(); it!=matrixContainer->elementColEnd(); ++it) {
            if(((*it)->getID()-1)/(numberOfElements/TotalAmountOfProcessors+1)==localProcessorNumber) {
                matrix.resize((*it)->getNrOfUnknows(),(*it)->getNrOfUnknows());
                elIntegral.integrate((*it),elF,matrix,matrixContainer);
                int places[]= {(*it)->getID()-1};
                for(int i=0; i<(*it)->getNrOfUnknows()*(*it)->getNrOfUnknows(); ++i) {
                    tempComplexArray[i]=matrix[i];//TODO Make complex PETSc read reals
                }
                matrixContainer->ierr_=MatSetValuesBlocked(matrixContainer->S_,1,places,1,places,&tempComplexArray[0],ADD_VALUES);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
            }
        }
        hpGemUIExtentions<3>::FaceFunction faF = &hpGemUIExtentions<3>::faceIntegrand;
        Integration::FaceIntegral<3> faIntegral(false);
        for(hpGemUIExtentions<3>::FaceIterator it=matrixContainer->faceColBegin(); it!=matrixContainer->faceColEnd(); ++it) {
            //faces dont have an ID; pick an arbitrary processor to do the work
            if(((*it)->getPtrElementLeft()->getID()-1)/(numberOfElements/TotalAmountOfProcessors+1)==localProcessorNumber) {
                if((*it)->isInternal()) {
                    matrix.resize((*it)->getPtrElementLeft()->getNrOfUnknows()+(*it)->getPtrElementRight()->getNrOfUnknows(),(*it)->getPtrElementLeft()->getNrOfUnknows()+(*it)->getPtrElementRight()->getNrOfUnknows());
                    int places[]= {(*it)->getPtrElementLeft()->getID()-1,(*it)->getPtrElementRight()->getID()-1};
                    faIntegral.integrate(*it,faF,matrix,matrixContainer);
                    for(int i=0; i<4*(*it)->getPtrElementLeft()->getNrOfUnknows()*(*it)->getPtrElementLeft()->getNrOfUnknows(); ++i) {
                        tempComplexArray[i]=matrix[i];//TODO Make complex PETSc read reals
                    }
                    matrixContainer->ierr_=MatSetValuesBlocked(matrixContainer->S_,2,places,2,places,&tempComplexArray[0],ADD_VALUES);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
                } else {
                    matrix.resize((*it)->getPtrElementLeft()->getNrOfUnknows(),(*it)->getPtrElementLeft()->getNrOfUnknows());
                    int places[]= {(*it)->getPtrElementLeft()->getID()-1};
                    faIntegral.integrate(*it,faF,matrix,matrixContainer);
                    for(int i=0; i<(*it)->getPtrElementLeft()->getNrOfUnknows()*(*it)->getPtrElementLeft()->getNrOfUnknows(); ++i) {
                        tempComplexArray[i]=matrix[i];//TODO Make complex PETSc read reals
                    }
                    matrixContainer->ierr_=MatSetValuesBlocked(matrixContainer->S_,1,places,1,places,&tempComplexArray[0],ADD_VALUES);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
                }
            }
        }
        faF = &hpGemUIExtentions<3>::faceIntegrandIPPart;
        for(hpGemUIExtentions<3>::FaceIterator it=matrixContainer->faceColBegin(); it!=matrixContainer->faceColEnd(); ++it) {
            if(((*it)->getPtrElementLeft()->getID()-1)/(numberOfElements/TotalAmountOfProcessors+1)==localProcessorNumber) {
                if((*it)->isInternal()) {
                    matrix.resize((*it)->getPtrElementLeft()->getNrOfUnknows()+(*it)->getPtrElementRight()->getNrOfUnknows(),(*it)->getPtrElementLeft()->getNrOfUnknows()+(*it)->getPtrElementRight()->getNrOfUnknows());
                    int places[]= {(*it)->getPtrElementLeft()->getID()-1,(*it)->getPtrElementRight()->getID()-1};
                    faIntegral.integrate(*it,faF,matrix,matrixContainer);
                    for(int i=0; i<4*(*it)->getPtrElementLeft()->getNrOfUnknows()*(*it)->getPtrElementLeft()->getNrOfUnknows(); ++i) {
                        tempComplexArray[i]=matrix[i];//TODO Make complex PETSc read reals
                    }
                    matrixContainer->ierr_=MatSetValuesBlocked(matrixContainer->S_,2,places,2,places,&tempComplexArray[0],ADD_VALUES);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
                } else {
                    matrix.resize((*it)->getPtrElementLeft()->getNrOfUnknows(),(*it)->getPtrElementLeft()->getNrOfUnknows());
                    int places[]= {(*it)->getPtrElementLeft()->getID()-1};
                    faIntegral.integrate(*it,faF,matrix,matrixContainer);
                    for(int i=0; i<(*it)->getPtrElementLeft()->getNrOfUnknows()*(*it)->getPtrElementLeft()->getNrOfUnknows(); ++i) {
                        tempComplexArray[i]=matrix[i];//TODO Make complex PETSc read reals
                    }
                    matrixContainer->ierr_=MatSetValuesBlocked(matrixContainer->S_,1,places,1,places,&tempComplexArray[0],ADD_VALUES);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
                }
            }
        }
	time(&newTime);
	cout<<"filling the stiffness matrix took "<<difftime(newTime,oldTime)<<" seconds"<<endl;
	oldTime=newTime;
        matrixContainer->ierr_=MatAssemblyBegin(matrixContainer->S_,MAT_FINAL_ASSEMBLY);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
        elF=&hpGemUIExtentions<3>::sourceTerm;
        for(hpGemUIExtentions<3>::ElementIterator it=matrixContainer->elementColBegin(); it!=matrixContainer->elementColEnd(); ++it) {
            if(((*it)->getID()-1)/(numberOfElements/TotalAmountOfProcessors+1)==localProcessorNumber) {
                matrix.resize((*it)->getNrOfUnknows(),1);
                elIntegral.integrate((*it),elF,matrix,matrixContainer);
                int places[]= {(*it)->getID()-1};
                for(int i=0; i<(*it)->getNrOfUnknows(); ++i) {
                    tempComplexArray[i]=matrix[i];//TODO Make complex PETSc read reals
                }
                matrixContainer->ierr_=VecSetValuesBlocked(matrixContainer->RHS_,1,places,&tempComplexArray[0],ADD_VALUES);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
            }
        }
        elF=&hpGemUIExtentions<3>::initialConditions;
        for(hpGemUIExtentions<3>::ElementIterator it=matrixContainer->elementColBegin(); it!=matrixContainer->elementColEnd(); ++it) {
            if(((*it)->getID()-1)/(numberOfElements/TotalAmountOfProcessors+1)==localProcessorNumber) {
                matrix.resize((*it)->getNrOfUnknows(),1);
                elIntegral.integrate((*it),elF,matrix,matrixContainer);
                int places[]= {(*it)->getID()-1};
                for(int i=0; i<(*it)->getNrOfUnknows(); ++i) {
                    tempComplexArray[i]=matrix[i];//TODO Make complex PETSc read reals
                }
                matrixContainer->ierr_=VecSetValuesBlocked(matrixContainer->x_,1,places,&tempComplexArray[0],ADD_VALUES);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
            }
        }
        matrixContainer->ierr_=VecAssemblyBegin(matrixContainer->x_);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
        elF=&hpGemUIExtentions<3>::initialConditionsDeriv;
        for(hpGemUIExtentions<3>::ElementIterator it=matrixContainer->elementColBegin(); it!=matrixContainer->elementColEnd(); ++it) {
            if(((*it)->getID()-1)/(numberOfElements/TotalAmountOfProcessors+1)==localProcessorNumber) {
                matrix.resize((*it)->getNrOfUnknows(),1);
                elIntegral.integrate((*it),elF,matrix,matrixContainer);
                int places[]= {(*it)->getID()-1};
                for(int i=0; i<(*it)->getNrOfUnknows(); ++i) {
                    tempComplexArray[i]=matrix[i];//TODO Make complex PETSc read reals
                }
                matrixContainer->ierr_=VecSetValuesBlocked(matrixContainer->derivative_,1,places,&tempComplexArray[0],ADD_VALUES);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
            }
        }
        matrixContainer->ierr_=VecAssemblyBegin(matrixContainer->derivative_);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
        faF = &hpGemUIExtentions<3>::sourceTermBoundaryIP;
        for(hpGemUIExtentions<3>::FaceIterator it=matrixContainer->faceColBegin(); it!=matrixContainer->faceColEnd(); ++it) {
            if(((*it)->getPtrElementLeft()->getID()-1)/(numberOfElements/TotalAmountOfProcessors+1)==localProcessorNumber) {
                if((*it)->isInternal()) {
                    //internal faces dont produce boundary contributions to the RHS
                } else {
                    matrix.resize((*it)->getPtrElementLeft()->getNrOfUnknows(),1);
                    int places[]= {(*it)->getPtrElementLeft()->getID()-1};
                    faIntegral.integrate(*it,faF,matrix,matrixContainer);
                    for(int i=0; i<(*it)->getPtrElementLeft()->getNrOfUnknows(); ++i) {
                        tempComplexArray[i]=matrix[i];//TODO Make complex PETSc read reals
                    }
                    matrixContainer->ierr_=VecSetValuesBlocked(matrixContainer->RHS_,1,places,&tempComplexArray[0],ADD_VALUES);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
                }
            }
        }
	time(&newTime);
	cout<<"filling the vectors took "<<difftime(newTime,oldTime)<<" seconds"<<endl;
	oldTime=newTime;
        matrixContainer->ierr_=VecAssemblyBegin(matrixContainer->RHS_);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);

        //this functions guarantees assembled matrices so lock untill assembly is finished
        matrixContainer->ierr_=VecAssemblyEnd(matrixContainer->x_);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
        matrixContainer->ierr_=VecAssemblyEnd(matrixContainer->derivative_);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
        matrixContainer->ierr_=MatAssemblyEnd(matrixContainer->M_,MAT_FINAL_ASSEMBLY);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
        matrixContainer->ierr_=MatAssemblyEnd(matrixContainer->S_,MAT_FINAL_ASSEMBLY);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
        matrixContainer->ierr_=VecAssemblyEnd(matrixContainer->RHS_);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
	time(&newTime);
	cout<<"communication took "<<difftime(newTime,oldTime)<<" extra seconds"<<endl;

}

void matrixFillerBR::fillMatrixes(hpGemUIExtentions< 3 >* matrixContainer){
    time_t oldTime,newTime;
    time(&oldTime);
    cout<<"using a Brezzi-flux with penalty parameter "<<matrixContainer->getData()->StabCoeff_<<endl;
    int TotalAmountOfProcessors,localProcessorNumber,numberOfElements;
    MPI_Comm_size(PETSC_COMM_WORLD,&TotalAmountOfProcessors);
    MPI_Comm_rank(PETSC_COMM_WORLD,&localProcessorNumber);
    numberOfElements=matrixContainer->getNumberOfElements(0);

//     if(matrixContainer->getNumberOfElements(0)%TotalAmountOfProcessors!=0){
// 	cout<<"WARNING: the case where the number of elements is not a multiple of the number of processors has not been toroughly tested; gliches may occur!";
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
    hpGemUIExtentions<3>::ElementFunction elF = &hpGemUIExtentions<3>::elementMassIntegrand;
    PetscScalar tempComplexArray[4*matrixContainer->getConfigData()->numberOfBasisFunctions_*matrixContainer->getConfigData()->numberOfBasisFunctions_];
    Integration::ElementIntegral<3> elIntegral(false);
    time(&newTime);
    cout<<"preparation took "<<difftime(newTime,oldTime)<<" seconds"<<endl;
    oldTime=newTime;
    for(hpGemUIExtentions<3>::ElementIterator it=matrixContainer->elementColBegin(); it!=matrixContainer->elementColEnd(); ++it) {
	if(((*it)->getID()-1)/(numberOfElements/TotalAmountOfProcessors+1)==localProcessorNumber) {
	    matrix.resize((*it)->getNrOfUnknows(),(*it)->getNrOfUnknows());
	    elIntegral.integrate((*it),elF,matrix,matrixContainer);
	    if(matrixContainer->MHasToBeInverted_) {
		matrix.inverse(matrix);
	    }
	    int places[]= {(*it)->getID()-1};
	    for(int i=0; i<(*it)->getNrOfUnknows()*(*it)->getNrOfUnknows(); ++i) {
		tempComplexArray[i]=matrix[i];//TODO Make complex PETSc read reals
	    }
	    matrixContainer->ierr_=MatSetValuesBlocked(matrixContainer->M_,1,places,1,places,&tempComplexArray[0],ADD_VALUES);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
	}
    }
    time(&newTime);
    cout<<"filling the mass matrix took "<<difftime(newTime,oldTime)<<" seconds"<<endl;
    oldTime=newTime;
    //we are done filling M. PETSc may start assembling M in the background
    matrixContainer->ierr_=MatAssemblyBegin(matrixContainer->M_,MAT_FINAL_ASSEMBLY);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
    elF=&hpGemUIExtentions<3>::elementStiffnessIntegrand;
    for(hpGemUIExtentions<3>::ElementIterator it=matrixContainer->elementColBegin(); it!=matrixContainer->elementColEnd(); ++it) {
	if(((*it)->getID()-1)/(numberOfElements/TotalAmountOfProcessors+1)==localProcessorNumber) {
	    matrix.resize((*it)->getNrOfUnknows(),(*it)->getNrOfUnknows());
	    elIntegral.integrate((*it),elF,matrix,matrixContainer);
	    int places[]= {(*it)->getID()-1};
	    for(int i=0; i<(*it)->getNrOfUnknows()*(*it)->getNrOfUnknows(); ++i) {
		tempComplexArray[i]=matrix[i];//TODO Make complex PETSc read reals
	    }
	    matrixContainer->ierr_=MatSetValuesBlocked(matrixContainer->S_,1,places,1,places,&tempComplexArray[0],ADD_VALUES);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
	}
    }
    hpGemUIExtentions<3>::FaceFunction faF = &hpGemUIExtentions<3>::faceIntegrand;
    Integration::FaceIntegral<3> faIntegral(false);
    for(hpGemUIExtentions<3>::FaceIterator it=matrixContainer->faceColBegin(); it!=matrixContainer->faceColEnd(); ++it) {
	if(((*it)->getPtrElementLeft()->getID()-1)/(numberOfElements/TotalAmountOfProcessors+1)==localProcessorNumber) {
	    if((*it)->isInternal()) {
		matrix.resize((*it)->getPtrElementLeft()->getNrOfUnknows()+(*it)->getPtrElementRight()->getNrOfUnknows(),(*it)->getPtrElementLeft()->getNrOfUnknows()+(*it)->getPtrElementRight()->getNrOfUnknows());
		int places[]= {(*it)->getPtrElementLeft()->getID()-1,(*it)->getPtrElementRight()->getID()-1};
		faIntegral.integrate(*it,faF,matrix,matrixContainer);
		for(int i=0; i<4*(*it)->getPtrElementLeft()->getNrOfUnknows()*(*it)->getPtrElementLeft()->getNrOfUnknows(); ++i) {
		    tempComplexArray[i]=matrix[i];//TODO Make complex PETSc read reals
		}
		matrixContainer->ierr_=MatSetValuesBlocked(matrixContainer->S_,2,places,2,places,&tempComplexArray[0],ADD_VALUES);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
	    } else {
		matrix.resize((*it)->getPtrElementLeft()->getNrOfUnknows(),(*it)->getPtrElementLeft()->getNrOfUnknows());
		int places[]= {(*it)->getPtrElementLeft()->getID()-1};
		faIntegral.integrate(*it,faF,matrix,matrixContainer);
		for(int i=0; i<(*it)->getPtrElementLeft()->getNrOfUnknows()*(*it)->getPtrElementLeft()->getNrOfUnknows(); ++i) {
		    tempComplexArray[i]=matrix[i];//TODO Make complex PETSc read reals
		}
		matrixContainer->ierr_=MatSetValuesBlocked(matrixContainer->S_,1,places,1,places,&tempComplexArray[0],ADD_VALUES);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
	    }
	}
    }
    //we will be needing a few entries of M so assembly should be done by now
    matrixContainer->ierr_=MatAssemblyEnd(matrixContainer->M_,MAT_FINAL_ASSEMBLY);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
    
    //make some assistant matrices //TODO look if both are neccecary
    Mat DInternal,DBoundary,MLocal,dummy,dummy2;
    matrixContainer->ierr_=MatCreateDense(PETSC_COMM_SELF,PETSC_DECIDE,PETSC_DECIDE,2*matrixContainer->getConfigData()->numberOfBasisFunctions_,2*matrixContainer->getConfigData()->numberOfBasisFunctions_,NULL,&DInternal);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
    matrixContainer->ierr_=MatCreateDense(PETSC_COMM_SELF,PETSC_DECIDE,PETSC_DECIDE,matrixContainer->getConfigData()->numberOfBasisFunctions_,matrixContainer->getConfigData()->numberOfBasisFunctions_,NULL,&DBoundary);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
    matrixContainer->ierr_=MatSetUp(DInternal);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
    matrixContainer->ierr_=MatSetUp(DBoundary);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
    int localnElements;
    int dimension,skip(0);
    faF = &hpGemUIExtentions<3>::faceIntegrandBRPart;
    for(hpGemUIExtentions<3>::FaceIterator it=matrixContainer->faceColBegin(); it!=matrixContainer->faceColEnd(); ++it) {
	int placesBlocked[]= {(*it)->getPtrElementLeft()->getID()-1,-1};
	if((*it)->isInternal()) {
	    localnElements=2;
	    dimension=(*it)->getPtrElementLeft()->getNrOfUnknows()+(*it)->getPtrElementRight()->getNrOfUnknows();
	    placesBlocked[1]=(*it)->getPtrElementRight()->getID()-1;
	} else {
	    localnElements=1;
	    dimension=(*it)->getPtrElementLeft()->getNrOfUnknows();
	}

	//PETSc likes to have only unique values in an IS when loacing from a submatrix so create some unique bogus
	//in the parts that are not needed
	if(((*it)->getPtrElementLeft()->getID()-1)/(numberOfElements/TotalAmountOfProcessors+1)!=localProcessorNumber) {
	    int temp1=placesBlocked[0];
	    placesBlocked[0]=4*localProcessorNumber;
	    while(temp1==placesBlocked[0]||placesBlocked[1]==placesBlocked[0]) {
		placesBlocked[0]++;
		placesBlocked[0]%=numberOfElements;
	    }
	    int temp2=placesBlocked[1];
	    placesBlocked[1]=placesBlocked[0]+1;
	    while(temp1==placesBlocked[1]||temp2==placesBlocked[1]) {
		placesBlocked[1]++;
		placesBlocked[1]%=numberOfElements;
	    }
	}
	IS ISplaces;
	matrixContainer->ierr_=ISCreateBlock(PETSC_COMM_WORLD,(*it)->getPtrElementLeft()->getNrOfUnknows(),localnElements,placesBlocked,PETSC_COPY_VALUES,&ISplaces);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
	matrixContainer->ierr_=MatGetSubMatrix(matrixContainer->M_,ISplaces,ISplaces,MAT_INITIAL_MATRIX,&MLocal);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
	if(((*it)->getPtrElementLeft()->getID()-1)/(numberOfElements/TotalAmountOfProcessors+1)==localProcessorNumber) {
	    matrix.resize(dimension,dimension);
	    faIntegral.integrate(*it,faF,matrix,matrixContainer);
	    for(int i=0; i<dimension*dimension; ++i) {
		tempComplexArray[i]=matrix[i];//TODO Make complex PETSc read reals
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
		faF = &hpGemUIExtentions<3>::sourceTermBoundaryBR;
		matrix.resize(dimension,1);
		faIntegral.integrate(*it,faF,matrix,matrixContainer);
		for(int i=0; i<dimension; ++i) {
		    tempComplexArray[i]=matrix[i];//TODO Make complex PETSc read reals
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
		faF = &hpGemUIExtentions<3>::faceIntegrandBRPart;
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
	}
    }
    time(&newTime);
    cout<<"filling the stiffness matrix took "<<difftime(newTime,oldTime)<<" seconds"<<endl;
    oldTime=newTime;
    elF=&hpGemUIExtentions<3>::sourceTerm;
    for(hpGemUIExtentions<3>::ElementIterator it=matrixContainer->elementColBegin(); it!=matrixContainer->elementColEnd(); ++it) {
	if(((*it)->getID()-1)/(numberOfElements/TotalAmountOfProcessors+1)==localProcessorNumber) {
	    matrix.resize((*it)->getNrOfUnknows(),1);
	    elIntegral.integrate((*it),elF,matrix,matrixContainer);
	    int places[]= {(*it)->getID()-1};
	    for(int i=0; i<(*it)->getNrOfUnknows(); ++i) {
		tempComplexArray[i]=matrix[i];//TODO Make complex PETSc read reals
	    }
	    matrixContainer->ierr_=VecSetValuesBlocked(matrixContainer->RHS_,1,places,&tempComplexArray[0],ADD_VALUES);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
	}
    }
    matrixContainer->ierr_=MatAssemblyBegin(matrixContainer->S_,MAT_FINAL_ASSEMBLY);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
    elF=&hpGemUIExtentions<3>::initialConditions;
    for(hpGemUIExtentions<3>::ElementIterator it=matrixContainer->elementColBegin(); it!=matrixContainer->elementColEnd(); ++it) {
	if(((*it)->getID()-1)/(numberOfElements/TotalAmountOfProcessors+1)==localProcessorNumber) {
	    matrix.resize((*it)->getNrOfUnknows(),1);
	    elIntegral.integrate((*it),elF,matrix,matrixContainer);
	    int places[]= {(*it)->getID()-1};
	    for(int i=0; i<(*it)->getNrOfUnknows(); ++i) {
		tempComplexArray[i]=matrix[i];//TODO Make complex PETSc read reals
	    }
	    matrixContainer->ierr_=VecSetValuesBlocked(matrixContainer->x_,1,places,&tempComplexArray[0],ADD_VALUES);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
	}
    }
    matrixContainer->ierr_=VecAssemblyBegin(matrixContainer->x_);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
    elF=&hpGemUIExtentions<3>::initialConditionsDeriv;
    for(hpGemUIExtentions<3>::ElementIterator it=matrixContainer->elementColBegin(); it!=matrixContainer->elementColEnd(); ++it) {
	if(((*it)->getID()-1)/(numberOfElements/TotalAmountOfProcessors+1)==localProcessorNumber) {
	    matrix.resize((*it)->getNrOfUnknows(),1);
	    elIntegral.integrate((*it),elF,matrix,matrixContainer);
	    int places[]= {(*it)->getID()-1};
	    for(int i=0; i<(*it)->getNrOfUnknows(); ++i) {
		tempComplexArray[i]=matrix[i];//TODO Make complex PETSc read reals
	    }
	    matrixContainer->ierr_=VecSetValuesBlocked(matrixContainer->derivative_,1,places,&tempComplexArray[0],ADD_VALUES);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
	}
    }
    matrixContainer->ierr_=VecAssemblyBegin(matrixContainer->derivative_);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
    faF = &hpGemUIExtentions<3>::sourceTermBoundary;
    for(hpGemUIExtentions<3>::FaceIterator it=matrixContainer->faceColBegin(); it!=matrixContainer->faceColEnd(); ++it){
	if(((*it)->getPtrElementLeft()->getID()-1)/(numberOfElements/TotalAmountOfProcessors+1)==localProcessorNumber) {
	    if((*it)->isInternal()) {
		//internal faces have no boundary contributions
	    } else {
		matrix.resize((*it)->getPtrElementLeft()->getNrOfUnknows(),1);
		int places[]= {(*it)->getPtrElementLeft()->getID()-1};
		faIntegral.integrate(*it,faF,matrix,matrixContainer);
		for(int i=0; i<(*it)->getPtrElementLeft()->getNrOfUnknows(); ++i) {
		    tempComplexArray[i]=matrix[i];//TODO Make complex PETSc read reals
		}
		matrixContainer->ierr_=VecSetValuesBlocked(matrixContainer->RHS_,1,places,&tempComplexArray[0],ADD_VALUES);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
	    }
	}
    }
    time(&newTime);
    cout<<"filling vectors took "<<difftime(newTime,oldTime)<<" seconds"<<endl;
    oldTime=newTime;
    matrixContainer->ierr_=VecAssemblyBegin(matrixContainer->RHS_);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
    matrixContainer->ierr_=VecAssemblyEnd(matrixContainer->x_);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
    matrixContainer->ierr_=VecAssemblyEnd(matrixContainer->derivative_);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
    matrixContainer->ierr_=MatAssemblyEnd(matrixContainer->S_,MAT_FINAL_ASSEMBLY);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
    matrixContainer->ierr_=VecAssemblyEnd(matrixContainer->RHS_);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);

    matrixContainer->ierr_=MatDestroy(&DInternal);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
    matrixContainer->ierr_=MatDestroy(&DBoundary);CHKERRABORT(PETSC_COMM_WORLD,matrixContainer->ierr_);
    time(&newTime);
    cout<<"communication took "<<difftime(newTime,oldTime)<<" extra seconds"<<endl;
}
