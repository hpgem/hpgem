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

//WARNING: THIS HEADER PARTLY DEPENDS ON C LIBRARIES - INCLUDING THIS HEADER BEFORE OTHER HEADERS GREATLY INCLEASES THE RISK OF NAMING CONFLICTS OR OTHER COMPILE TIME ERRORS
#ifndef GLOBALMATRIX_HPP_
#define GLOBALMATRIX_HPP_

#if defined(HPGEM_USE_PETSC) || defined(HPGEM_USE_COMPLEX_PETSC)
#include <petscmat.h>
#endif
#include <vector>

namespace Base
{
    class MeshManipulator;
    class Face;
    class Element;
}

namespace Utilities
{
    
    /**
     *  General global assembly class for matrices. Defines the general routines that specializations will always need
     *  Subclasses should also derive from the external package specific matrix type they try to implement
     */
    class GlobalMatrix
    {
        
    public:
        
        ///use the deconstructor of the subclass in case this is needed
        virtual ~GlobalMatrix()
        {
        }
        
        ///constructs the global matrix and performs element assembly
        GlobalMatrix(Base::MeshManipulator* theMesh, int elementMatrixID, int faceMatrixID);

        ///signals the matrix that (some) element matrixes have changed and that the matrix can be assembled again
        ///also used after mesh refinement to make sure the global matrix reflects the changes in the mesh
        virtual void reAssemble()=0;

        ///cleans out all the entries, putting them back to 0 - keeps the non-zero strucure of the matrix however.
        ///After this it assembles the matrix to put it back in a legal state
        virtual void reset()=0;

        void getMatrixBCEntries(const Base::Face* face, int& numberOfEntries, std::vector<int>& entries);

    protected:
        
        int meshLevel_, elementMatrixID_, faceMatrixID_;
        std::vector<int> startPositionsOfElementsInTheMatrix_;
        std::vector<int> startPositionsOfFacesInTheMatrix_;
        std::vector<int> startPositionsOfEdgesInTheMatrix_;
        std::vector<int> startPositionsOfVerticesInTheMatrix_;
        Base::MeshManipulator *theMesh_;
        
    };
#if defined(HPGEM_USE_PETSC) || defined(HPGEM_USE_COMPLEX_PETSC)
    ///\bug this class depends on PETSc and is likely to cause naming conflicts between the c and c++ standard libraries (workaround: make sure to include all other needed hpGEM headers before including this header)
    class GlobalPetscMatrix : public GlobalMatrix
    {
        
    public:
        ///for now provides implicit conversion to Mat (the PETSc matrix type)
        ///this needs a special function because deriving from Mat appears to be illegal
        ///\bug need a better way to provide an interface to the supported Mat routines AND to other routines that need a Mat (like KSPSetOperators()) (but not stuff like MatDestroy())
        operator Mat();

        GlobalPetscMatrix(Base::MeshManipulator* theMesh, int elementMatrixID, int faceMatrixID = -1);
        ~GlobalPetscMatrix();

        void reset();

        void reAssemble();

    private:
        
        std::vector<PetscInt> makePositionsInMatrix(const Base::Element*);

    private:
        
        Mat A_;
    };
#endif
}

#endif /* GLOBALMATRIX_HPP_ */
