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

#include "Base/MpiContainer.h"
#include "GlobalVector.h"
#include <vector>
#include "Base/MeshManipulatorBase.h"
#include "Base/Edge.h"
#include "Base/Face.h"
#include "Base/Element.h"
#include "Base/ElementCacheData.h"
#include "LinearAlgebra/MiddleSizeVector.h"
#include "Base/FaceCacheData.h"
#include "Geometry/PointPhysical.h"
#include "Geometry/PhysicalGeometry.h"
#include "Geometry/ReferenceGeometry.h"
#include "Geometry/PointReference.h"
#include <Logger.h>
#include <numeric>
#if defined(HPGEM_USE_PETSC) || defined(HPGEM_USE_COMPLEX_PETSC)
    #include "petscis.h"
#endif
#if defined(HPGEM_USE_SUNDIALS)
	#include "nvector/nvector_serial.h"
#endif

namespace Utilities
{
    
    GlobalVector::GlobalVector(Base::MeshManipulatorBase* theMesh, int elementVectorID, int faceVectorID)
            : meshLevel_(-2), elementVectorID_(elementVectorID), faceVectorID_(faceVectorID), startPositionsOfElementsInTheVector_(), theMesh_(theMesh)
    {
        logger.assert_debug(theMesh != nullptr, "Invalid mesh passed");
    }
    
#if defined(HPGEM_USE_PETSC) || defined(HPGEM_USE_COMPLEX_PETSC)
    
    GlobalPetscVector::GlobalPetscVector(Base::MeshManipulatorBase* theMesh, int elementVectorID, int faceVectorID)
            : GlobalVector(theMesh, elementVectorID, faceVectorID)
            , indexing_() // Will be initialized by reset().
    {
        logger.assert_debug(theMesh != nullptr, "Invalid mesh passed");
        PetscBool petscRuns;
        PetscInitialized(&petscRuns);
        logger.assert_debug(petscRuns == PETSC_TRUE, "Early call, firstly the command line arguments should be parsed");
        VecCreateSeq(PETSC_COMM_SELF, 1, &b_);
        
        reset();
    }
    
    GlobalPetscVector::~GlobalPetscVector()
    {
        int ierr = VecDestroy(&b_);
        CHKERRV(ierr);
    }
    
    GlobalPetscVector::operator Vec()
    {
        if(HPGEM_LOGLEVEL>=Log::DEBUG)
        {
            VecChop(b_, 1e-13);
            VecScale(b_, 9.);
            VecView(b_, PETSC_VIEWER_STDOUT_WORLD);
            VecScale(b_, 1. / 9.);
        }
        return b_;
    }

    //debug note: GlobalPetscMatrix 'independently' chooses an ordering for the degrees of freedom, but hpGEM assumes both orderings to be the same
    void GlobalPetscVector::reset()
    {
        int ierr = VecDestroy(&b_);
        CHKERRV(ierr);
        indexing_.reset(theMesh_, GlobalIndexing::BLOCKED_PROCESSOR);
        ierr = VecCreateMPI(PETSC_COMM_WORLD, indexing_.getNumberOfLocalBasisFunctions(), PETSC_DETERMINE, &b_);
        CHKERRV(ierr);
    }
    
    void GlobalPetscVector::assemble()
    {
        reset();
        
        LinearAlgebra::MiddleSizeVector elementVector;
        std::vector<PetscInt> elementToGlobal (0);
        
        if (elementVectorID_ >= 0)
        {
            for (Base::Element* element : theMesh_->getElementsList())
            {
                indexing_.getGlobalIndices(element, elementToGlobal);
                elementVector = element->getElementVector(elementVectorID_);
                int ierr = VecSetValues(b_, elementToGlobal.size(), elementToGlobal.data(), elementVector.data(), ADD_VALUES);
                CHKERRV(ierr);
                
            }
        }
        
        LinearAlgebra::MiddleSizeVector faceVector;
        std::vector<PetscInt> faceToGlobal (0);
        if (faceVectorID_ >= 0)
        {
            for (Base::Face* face : theMesh_->getFacesList())
            {
                if (!face->isOwnedByCurrentProcessor())
                    continue;

                faceToGlobal.clear();
                indexing_.getGlobalIndices(face, faceToGlobal);
                faceVector = face->getFaceVector(faceVectorID_);
                int ierr = VecSetValues(b_, faceToGlobal.size(), faceToGlobal.data(), faceVector.data(), ADD_VALUES);
                CHKERRV(ierr);
            }
        }
        
        int ierr = VecAssemblyBegin(b_);
        ierr = VecAssemblyEnd(b_);
        CHKERRV(ierr);
    }
    
    void GlobalPetscVector::constructFromTimeIntegrationVector(std::size_t timeIntegrationVectorId, std::size_t solutionVar)
    {
        reset();
        
        LinearAlgebra::MiddleSizeVector elementData;
        std::vector<PetscInt> localToGlobal;
        for (Base::Element* element : theMesh_->getElementsList())
        {
            indexing_.getGlobalIndices(element, localToGlobal);
            elementData.resize(element->getTotalNumberOfBasisFunctions());
            for(std::size_t j = 0; j < element->getNumberOfUnknowns(); ++j)
            {
                for (std::size_t i = 0; i < element->getNumberOfBasisFunctions(j); ++i)
                {

                    if(j == solutionVar)
                    {
                        elementData[element->convertToSingleIndex(i, solutionVar)] = element->getTimeIntegrationData(timeIntegrationVectorId, solutionVar, i);
                    }
                    else
                    {
                        localToGlobal[element->convertToSingleIndex(i, j)] = -1;
                    }

                }
            }
            int ierr = VecSetValues(b_, localToGlobal.size(), localToGlobal.data(), elementData.data(), INSERT_VALUES);
            CHKERRV(ierr);
        }
        
        int ierr = VecAssemblyBegin(b_);
        ierr = VecAssemblyEnd(b_);
        CHKERRV(ierr);
    }
    
    void GlobalPetscVector::constructFromTimeIntegrationVector(std::size_t timeIntegrationVectorId)
    {
        reset();

        LinearAlgebra::MiddleSizeVector elementData;
        std::vector<PetscInt> localToGlobal;
        for (Base::Element* element : theMesh_->getElementsList())
        {
            std::size_t numberOfBasisFunctions = element->getTotalNumberOfBasisFunctions();
            indexing_.getGlobalIndices(element, localToGlobal);
            int ierr = VecSetValues(b_, numberOfBasisFunctions, localToGlobal.data(), element->getTimeIntegrationVector(timeIntegrationVectorId).data(), INSERT_VALUES);
            CHKERRV(ierr);
        }

        int ierr = VecAssemblyBegin(b_);
        ierr = VecAssemblyEnd(b_);
        CHKERRV(ierr);
    }

    void GlobalPetscVector::writeTimeIntegrationVector(std::size_t timeIntegrationVectorId)
    {
        PetscScalar *data;
        
        VecScatter scatter;
        Vec localB;
        //create a local vector...
        VecScatterCreateToAll(b_, &scatter, &localB);
        //we dont need the scatter context, just the vector
        VecScatterDestroy(&scatter);
        
        std::vector<PetscInt> positions;
        std::size_t totalPositions = 0;
        for(Base::Element* element : theMesh_->getElementsList())
        {
            totalPositions += element->getTotalNumberOfBasisFunctions();
        }
        positions.reserve(totalPositions);
        std::vector<PetscInt> newPositions;
        for (Base::Element* element : theMesh_->getElementsList())
        {
            indexing_.getGlobalIndices(element, newPositions);
            for (auto& a : newPositions)
            {
                positions.push_back(a);
            }
        }
        
        IS scatterIS;
        ISCreateGeneral(PETSC_COMM_SELF, positions.size(), positions.data(), PETSC_COPY_VALUES, &scatterIS);
        ISSortRemoveDups(scatterIS);
        VecScatterCreate(b_, scatterIS, localB, scatterIS, &scatter);
        VecScatterBegin(scatter, b_, localB, INSERT_VALUES, SCATTER_FORWARD);
        VecScatterEnd(scatter, b_, localB, INSERT_VALUES, SCATTER_FORWARD);
        VecScatterDestroy(&scatter);
        ISDestroy(&scatterIS);
        
        int ierr = VecGetArray(localB, &data);
        CHKERRV(ierr);
        //Question: wat doet deze lijn
        for (Base::MeshManipulatorBase::ElementIterator it = theMesh_->elementColBegin(); it != theMesh_->elementColEnd(); ++it)
        {
        	//Create vector for the local data that has to be written
            LinearAlgebra::MiddleSizeVector localData((*it)->getTotalNumberOfBasisFunctions());
            std::size_t runningTotal = 0;
            for(std::size_t index = 0; index < (*it)->getNumberOfUnknowns(); ++index) // for every variable iV
            {
                std::size_t nElementBasis = (*it)->getLocalNumberOfBasisFunctions(index);
                int elementBasis0 = indexing_.getGlobalIndex((*it), index);
                for(std::size_t i = 0; i < nElementBasis; ++i) // Get the local basis functions of the element
                {
                	//Copy the values from data to the localData and update the running number
                    localData[runningTotal] = std::real(data[elementBasis0 + i]);
                    ++runningTotal;
                }

                for (std::size_t i = 0; i < (*it)->getPhysicalGeometry()->getNumberOfFaces(); ++i) // for all faces of the element
                {
                    std::size_t nFaceBasis = (*it)->getFace(i)->getLocalNumberOfBasisFunctions(index);
                    int faceBasis0 = indexing_.getGlobalIndex((*it)->getFace(i), index);
                    for (std::size_t j = 0; j < nFaceBasis; ++j) // get local basis functions of a face
                    {
                    	//Copy the values from data to the localData and update the running number
                        localData[runningTotal] = std::real(data[faceBasis0 + j]);
                        ++runningTotal;
                    }
                }
                for (std::size_t i = 0; i < (*it)->getNumberOfEdges(); ++i) // For all edges of the element
                {
                    std::size_t nEdgeBasis = (*it)->getEdge(i)->getLocalNumberOfBasisFunctions(index);
                    int edgeBasis0 = indexing_.getGlobalIndex((*it)->getEdge(i), index);
                    for (std::size_t j = 0; j < nEdgeBasis; ++j) //Get the local basis function of an edge
                    {
                    	//Copy the values from data to the localData and update the running number

                        localData[runningTotal] = std::real(data[edgeBasis0 + j]);
                        ++runningTotal;
                    }
                }
                if (theMesh_->dimension() > 1) // There are no nodes in a 1D problem
                {
                    for (std::size_t i = 0; i < (*it)->getNumberOfNodes(); ++i) //For all nodes
                    {
                        std::size_t nNodeBasis = (*it)->getNode(i)->getLocalNumberOfBasisFunctions(index);
                        int nodeBasis0 = indexing_.getGlobalIndex((*it)->getNode(i), index);
                        for (std::size_t j = 0; j < nNodeBasis; ++j) //Get the local number of basis function of a node
                        {
                        	//Copy the values from data to the localData and update the running number
                            localData[runningTotal] = std::real(data[nodeBasis0 + j]);
                            ++runningTotal;
                        }
                    }
                }
            }
            //Put the localData in the element
            logger.assert_debug(localData.size() == runningTotal, "not enough info to fill the vector");
            (*it)->setTimeIntegrationVector(timeIntegrationVectorId, localData);
        }
        ierr = VecRestoreArray(localB, &data);
        VecDestroy(&localB);
        CHKERRV(ierr);
    }

    void GlobalPetscVector::writeTimeIntegrationVector(std::size_t timeIntegrationVectorId, std::size_t unknown)
    {
        PetscScalar *data;

        VecScatter scatter;
        Vec localB;
        //create a local vector...
        VecScatterCreateToAll(b_, &scatter, &localB);
        //we dont need the scatter context, just the vector
        VecScatterDestroy(&scatter);

        std::vector<PetscInt> positions;
        std::size_t totalPositions = 0;
        for(Base::Element* element : theMesh_->getElementsList())
        {
            totalPositions += element->getTotalNumberOfBasisFunctions();
        }
        positions.reserve(totalPositions);
        std::vector<PetscInt> localPositions;
        for (Base::Element* element : theMesh_->getElementsList())
        {
            indexing_.getGlobalIndices(element, localPositions);
            for (auto& a : localPositions)
            {
                positions.push_back(a);
            }
        }

        IS scatterIS;
        ISCreateGeneral(PETSC_COMM_SELF, positions.size(), positions.data(), PETSC_COPY_VALUES, &scatterIS);
        ISSortRemoveDups(scatterIS);
        VecScatterCreate(b_, scatterIS, localB, scatterIS, &scatter);
        VecScatterBegin(scatter, b_, localB, INSERT_VALUES, SCATTER_FORWARD);
        VecScatterEnd(scatter, b_, localB, INSERT_VALUES, SCATTER_FORWARD);
        VecScatterDestroy(&scatter);
        ISDestroy(&scatterIS);

        int ierr = VecGetArray(localB, &data);
        CHKERRV(ierr);
        for (Base::MeshManipulatorBase::ElementIterator it = theMesh_->elementColBegin(); it != theMesh_->elementColEnd(); ++it)
        {
            LinearAlgebra::MiddleSizeVector localData((*it)->getTotalNumberOfBasisFunctions());
            std::size_t runningTotal = 0;
            for(std::size_t index = 0; index < (*it)->getNumberOfUnknowns(); ++index)
            {
                std::size_t nElementBasis = (*it)->getLocalNumberOfBasisFunctions(index);
                int elementBasis0 = indexing_.getGlobalIndex((*it), index);
                for(std::size_t i = 0; i < nElementBasis; ++i)
                {
                    localData[runningTotal] = std::real(data[elementBasis0 + i]);
                    ++runningTotal;
                }

                for (std::size_t i = 0; i < (*it)->getPhysicalGeometry()->getNumberOfFaces(); ++i)
                {
                    std::size_t nFaceBasis = (*it)->getFace(i)->getLocalNumberOfBasisFunctions(index);
                    int faceBasis0 = indexing_.getGlobalIndex((*it)->getFace(i), index);
                    for (std::size_t j = 0; j < nFaceBasis; ++j)
                    {
                        localData[runningTotal] = std::real(data[faceBasis0 + j]);
                        ++runningTotal;
                    }
                }
                for (std::size_t i = 0; i < (*it)->getNumberOfEdges(); ++i)
                {
                    std::size_t nEdgeBasis = (*it)->getEdge(i)->getLocalNumberOfBasisFunctions(index);
                    int edgeBasis0 = indexing_.getGlobalIndex((*it)->getEdge(i), index);
                    for (std::size_t j = 0; j < nEdgeBasis; ++j)
                    {
                        localData[runningTotal] = std::real(data[edgeBasis0 + j]);
                        ++runningTotal;
                    }
                }
                if (theMesh_->dimension() > 1)
                {
                    for (std::size_t i = 0; i < (*it)->getNumberOfNodes(); ++i)
                    {
                        std::size_t nNodeBasis = (*it)->getNode(i)->getLocalNumberOfBasisFunctions(index);
                        int nodeBasis0 = indexing_.getGlobalIndex((*it)->getNode(i), index);
                        for (std::size_t j = 0; j < nNodeBasis; ++j)
                        {
                            localData[runningTotal] = std::real(data[nodeBasis0 + j]);
                            ++runningTotal;
                        }
                    }
                }
            }
            logger.assert_debug(localData.size() == runningTotal, "not enough info to fill the vector");
            LinearAlgebra::MiddleSizeVector singleUnknownData((*it)->getNumberOfBasisFunctions());
            for(std::size_t i = 0; i < (*it)->getNumberOfBasisFunctions(); ++i)
            {
                singleUnknownData[i] = localData[(*it)->convertToSingleIndex(i, unknown)];
            }
            (*it)->setTimeIntegrationSubvector(timeIntegrationVectorId, unknown, singleUnknownData);
        }
        ierr = VecRestoreArray(localB, &data);
        VecDestroy(&localB);
        CHKERRV(ierr);
    }
#endif


#if defined(HPGEM_USE_SUNDIALS)

    //Constructor that will initialise the vector with zero entries but the correct length of the problem, nDOF
    GlobalSundialsVector::GlobalSundialsVector(Base::MeshManipulatorBase* theMesh, int elementVectorID, int faceVectorID)
            : GlobalVector(theMesh, elementVectorID, faceVectorID)
    {
        logger.assert_debug(theMesh!=nullptr, "Invalid mesh passed");
        //b_ = N_VNew_Serial(1); //Create a temporary placeholder vector
        reset();
    }

    GlobalSundialsVector::~GlobalSundialsVector()
    {
    	//N_VDestroy_Serial(b_);
    }

    GlobalSundialsVector::operator N_Vector()
    {
        return b_;
    }

    /// \brief Resets the N_Vector
    void GlobalSundialsVector::reset()
    {
        //Destroy old vector
    	//N_VDestroy_Serial(b_);

    	//Assign private variables to the global vector
    	//Initialise totalNrOfDOF to zero, this will be calculated as we go
        std::size_t totalNumberOfDOF(0), DIM(theMesh_->dimension());
        //Resize the startPosition vectors
        startPositionsOfElementsInTheVector_.resize(theMesh_->getNumberOfElements(Base::IteratorType::GLOBAL));
        startPositionsOfFacesInTheVector_.resize(theMesh_->getNumberOfFaces(Base::IteratorType::GLOBAL));
        startPositionsOfEdgesInTheVector_.resize(theMesh_->getNumberOfEdges(Base::IteratorType::GLOBAL));
        startPositionsOfNodesInTheVector_.resize(theMesh_->getNumberOfNodes(Base::IteratorType::GLOBAL));

        for (Base::Element* element : theMesh_->getElementsList())
        {
        	//Append the start of a new element to the vector, Then update the new total number of DOF
            startPositionsOfElementsInTheVector_[element->getID()] = totalNumberOfDOF; //New element position is appended to the vector
            totalNumberOfDOF += element->getLocalNumberOfBasisFunctions() * element->getNumberOfUnknowns(); //Calculate the new total DOF

            for (std::size_t i = 0; i < element->getNumberOfFaces(); ++i)
            {
                //faces at the boundary of the subdomain should also be added only once, so add them here is the left element of the face is in the subdomain
                if ((element->getFace(i)->getFaceType() == Geometry::FaceType::SUBDOMAIN_BOUNDARY || element->getFace(i)->getFaceType() == Geometry::FaceType::PERIODIC_SUBDOMAIN_BC)
                        && element->getFace(i)->getPtrElementLeft() == element)
                {
                    startPositionsOfFacesInTheVector_[element->getFace(i)->getID()] = totalNumberOfDOF;
                    totalNumberOfDOF += element->getFace(i)->getLocalNumberOfBasisFunctions() * element->getNumberOfUnknowns();
                    std::cout << "face local basis functions: " << element->getFace(i)->getLocalNumberOfBasisFunctions() << std::endl;
                }

            }
        }
        for (Base::Face* face : theMesh_->getFacesList())
        {
            //skip faces at the subdomain boundary because we already treated them
            if (face->getFaceType() != Geometry::FaceType::SUBDOMAIN_BOUNDARY && face->getFaceType() != Geometry::FaceType::PERIODIC_SUBDOMAIN_BC)
            {
                startPositionsOfFacesInTheVector_[face->getID()] = totalNumberOfDOF;
                totalNumberOfDOF += face->getLocalNumberOfBasisFunctions() * face->getPtrElementLeft()->getNumberOfUnknowns();
            }
        }
        for (Base::Edge* edge : theMesh_->getEdgesList())
        {
            startPositionsOfEdgesInTheVector_[edge->getID()] = totalNumberOfDOF;
            totalNumberOfDOF += edge->getLocalNumberOfBasisFunctions() * edge->getElement(0)->getNumberOfUnknowns();
        }
        //DIM == 1 faces and nodes are the same entities, skip one of them
        if (DIM > 1)
        {
            for (Base::Node* node : theMesh_->getNodesList())
            {
                startPositionsOfNodesInTheVector_[node->getID()] = totalNumberOfDOF;
                totalNumberOfDOF += node->getLocalNumberOfBasisFunctions() * node->getElement(0)->getNumberOfUnknowns();
            }
        }

        //b_ = N_VNew_Serial(totalNumberOfDOF);
        totalNumberOfDOF_ = totalNumberOfDOF;
    }

    void GlobalSundialsVector::constructFromTimeIntegrationVector(std::size_t timeIntegrationVectorId)
   {
       	//Extract data pointer from the input vector
    	double *data = NV_DATA_S(b_);

    		//For all elements
            for (Base::Element* element : theMesh_->getElementsList())
            {
            	//Create vector for the local data that has to be read
                LinearAlgebra::MiddleSizeVector localData = element->getTimeIntegrationVector(timeIntegrationVectorId);
                std::size_t runningTotal = 0;
                for(std::size_t index = 0; index < element->getNumberOfUnknowns(); ++index) // for every variable iV
                {
                    for(std::size_t i = 0; i < element->getLocalNumberOfBasisFunctions(); ++i) // Get the local basis functions of the element
                    {
                    	//Copy the values from localData to data and update the running number
                        data[startPositionsOfElementsInTheVector_[element->getID()] + i + index * element->getLocalNumberOfBasisFunctions()] = localData[runningTotal];
                        ++runningTotal;

                    }
                    // The code below does not return anything in case of DG, however in case of conforming FEM they will return data
                    for (std::size_t i = 0; i < element->getPhysicalGeometry()->getNumberOfFaces(); ++i) // for all faces of the element
                    {
                        for (std::size_t j = 0; j < element->getFace(i)->getLocalNumberOfBasisFunctions(); ++j) // get local basis functions of a face
                        {
                            data[startPositionsOfFacesInTheVector_[element->getFace(i)->getID()] + j + index * element->getFace(i)->getLocalNumberOfBasisFunctions()] = localData[runningTotal];
                            ++runningTotal;
                        }
                    }
                    for (std::size_t i = 0; i < element->getNumberOfEdges(); ++i) // For all edges of the element
                    {
                        for (std::size_t j = 0; j < element->getEdge(i)->getLocalNumberOfBasisFunctions(); ++j) //Get the local basis function of an edge
                        {
                        	//Copy the values from localData to data and update the running number
                        	data[startPositionsOfEdgesInTheVector_[element->getEdge(i)->getID()] + j + index * element->getEdge(i)->getLocalNumberOfBasisFunctions()] = localData[runningTotal];
                            ++runningTotal;
                        }
                    }
                    if (theMesh_->dimension() > 1) // There are no nodes in a 1D problem
                    {
                        for (std::size_t i = 0; i < element->getNumberOfNodes(); ++i) //For all nodes
                        {
                            for (std::size_t j = 0; j < element->getNode(i)->getLocalNumberOfBasisFunctions(); ++j) //Get the local number of basis function of a node
                            {
                            	//Copy the values from localData to data and update the running number
                            	data[startPositionsOfNodesInTheVector_[element->getNode(i)->getID()] + j + index * element->getNode(i)->getLocalNumberOfBasisFunctions()] = localData[runningTotal];
                                ++runningTotal;
                            }
                        }
                    }
                }
                logger.assert_debug(localData.size() == runningTotal, "not enough info to fill the vector");
            }
    }

    void GlobalSundialsVector::setScale(LinearAlgebra::MiddleSizeVector scaleFactor)
    {
       	//Extract data pointer from the input vector
    	double *data = NV_DATA_S(b_);

    	//For all elements
    	for (Base::Element* element : theMesh_->getElementsList())
    	{
    		std::size_t runningTotal = 0;
    		for(std::size_t index = 0; index < element->getNumberOfUnknowns(); ++index) // for every variable iV
    		{
    			for(std::size_t i = 0; i < element->getLocalNumberOfBasisFunctions(); ++i) // Get the local basis functions of the element
    			{
    				//Copy the scale value
    				data[startPositionsOfElementsInTheVector_[element->getID()] + i + index * element->getLocalNumberOfBasisFunctions()] = scaleFactor(index);
    				++runningTotal;

    			}
    			// The code below does not return anything in case of DG, however in case of conforming FEM they will return data
    			for (std::size_t i = 0; i < element->getPhysicalGeometry()->getNumberOfFaces(); ++i) // for all faces of the element
    			{
    				for (std::size_t j = 0; j < element->getFace(i)->getLocalNumberOfBasisFunctions(); ++j) // get local basis functions of a face
    				{
    					data[startPositionsOfFacesInTheVector_[element->getFace(i)->getID()] + j + index * element->getFace(i)->getLocalNumberOfBasisFunctions()] = scaleFactor(index);
    					++runningTotal;
    				}
    			}
    			for (std::size_t i = 0; i < element->getNumberOfEdges(); ++i) // For all edges of the element
    			{
    				for (std::size_t j = 0; j < element->getEdge(i)->getLocalNumberOfBasisFunctions(); ++j) //Get the local basis function of an edge
    				{
    					//Copy the scale value
    					data[startPositionsOfEdgesInTheVector_[element->getEdge(i)->getID()] + j + index * element->getEdge(i)->getLocalNumberOfBasisFunctions()] = scaleFactor(index);
    					++runningTotal;
    				}
    			}
    			if (theMesh_->dimension() > 1) // There are no nodes in a 1D problem
    			{
    				for (std::size_t i = 0; i < element->getNumberOfNodes(); ++i) //For all nodes
    				{
    					for (std::size_t j = 0; j < element->getNode(i)->getLocalNumberOfBasisFunctions(); ++j) //Get the local number of basis function of a node
    					{
    						//Copy the scale value
    						data[startPositionsOfNodesInTheVector_[element->getNode(i)->getID()] + j + index * element->getNode(i)->getLocalNumberOfBasisFunctions()] = scaleFactor(index);
    						++runningTotal;
    					}
    				}
    			}
    		}
    	}
    }

    void GlobalSundialsVector::writeTimeIntegrationVector(std::size_t timeIntegrationVectorId)
    {
    	//do something here
    	double *data = NV_DATA_S(b_);

        for (Base::Element* element : theMesh_->getElementsList())
        {
        	//Create vector for the local data that has to be written
            LinearAlgebra::MiddleSizeVector localData(element->getNumberOfBasisFunctions() * element->getNumberOfUnknowns());
            std::size_t runningTotal = 0;
            for(std::size_t index = 0; index < element->getNumberOfUnknowns(); ++index) // for every variable iV
            {
                for(std::size_t i = 0; i < element->getLocalNumberOfBasisFunctions(); ++i) // Get the local basis functions of the element
                {
                	//Copy the values from data to the localData and update the running number
                    localData[runningTotal] = std::real(data[startPositionsOfElementsInTheVector_[element->getID()] + i + index * element->getLocalNumberOfBasisFunctions()]);
                    ++runningTotal;
                }
                // The code below does not return anything in case of DG, however in case of conforming FEM they will return data
                for (std::size_t i = 0; i < element->getPhysicalGeometry()->getNumberOfFaces(); ++i) // for all faces of the element
                {
                    for (std::size_t j = 0; j < element->getFace(i)->getLocalNumberOfBasisFunctions(); ++j) // get local basis functions of a face
                    {
                    	//Copy the values from data to the localData and update the running number
                        localData[runningTotal] = std::real(data[startPositionsOfFacesInTheVector_[element->getFace(i)->getID()] + j + index * element->getFace(i)->getLocalNumberOfBasisFunctions()]);
                        ++runningTotal;
                    }
                }
                for (std::size_t i = 0; i < element->getNumberOfEdges(); ++i) // For all edges of the element
                {
                    for (std::size_t j = 0; j < element->getEdge(i)->getLocalNumberOfBasisFunctions(); ++j) //Get the local basis function of an edge
                    {
                    	//Copy the values from data to the localData and update the running number
                        localData[runningTotal] = std::real(data[startPositionsOfEdgesInTheVector_[element->getEdge(i)->getID()] + j + index * element->getEdge(i)->getLocalNumberOfBasisFunctions()]);
                        ++runningTotal;
                    }
                }
                if (theMesh_->dimension() > 1) // There are no nodes in a 1D problem
                {
                    for (std::size_t i = 0; i < element->getNumberOfNodes(); ++i) //For all nodes
                    {
                        for (std::size_t j = 0; j < element->getNode(i)->getLocalNumberOfBasisFunctions(); ++j) //Get the local number of basis function of a node
                        {
                        	//Copy the values from data to the localData and update the running number
                            localData[runningTotal] = std::real(data[startPositionsOfNodesInTheVector_[element->getNode(i)->getID()] + j + index * element->getNode(i)->getLocalNumberOfBasisFunctions()]);
                            ++runningTotal;
                        }
                    }
                }
            }
            //Put the localData in the element
            logger.assert_debug(localData.size() == runningTotal, "not enough info to fill the vector");
            element->setTimeIntegrationVector(timeIntegrationVectorId, localData);
        }
    }


    void GlobalSundialsVector::print()
    {
    	std::size_t length = NV_LENGTH_S(b_);
    	double *data = NV_DATA_S(b_);

		std::cout << "[" << data[0];
    	for (std::size_t i = 1; i < length; i++)
    	{
    		std::cout << ", " << data[i];
    	}
    	std::cout << ']' << std::endl;
    }

    //Possibly deprecated
/*    LinearAlgebra::MiddleSizeVector GlobalSundialsVector::getLocalVector(const Base::Element* ptrElement)
    {
    	//Create the local vector
    	std::size_t numberOfDOF = ptrElement->getNumberOfBasisFunctions() * ptrElement->getNumberOfUnknowns();
    	LinearAlgebra::MiddleSizeVector localVector(numberOfDOF);

    	//Fill the local vector
    	for (std::size_t i = 0; i < numberOfDOF; i++)
    	{
    		localVector(i) = startPositionsOfElementsInTheVector_[ptrElement->getID()] + i;
    	}

    	return localVector;
    }*/

    void GlobalSundialsVector::setVector(N_Vector b)
    {
    	b_ = b;
    }

#endif

}

