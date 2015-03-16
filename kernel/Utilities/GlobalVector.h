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
#ifndef GLOBALVECTOR_HPP_
#define GLOBALVECTOR_HPP_

#if defined(HPGEM_USE_PETSC) || defined(HPGEM_USE_COMPLEX_PETSC)
#include "petscvec.h"
#endif
#include <vector>
namespace Base
{
    class MeshManipulator;
    class Element;
}

namespace Utilities
{
    
    /**
     *  General global assembly class for vectors. Defines the general routines that specializations will always need
     *  Subclasses should also derive from the external package specific vector type they try to implement
     */
    class GlobalVector
    {
        
    public:
        ///use the destructor of the subclass in case this is needed
        virtual ~GlobalVector()
        {
        }
        
        ///construct the global vector does not do assembly by default
        ///because some vectors (like the solution of the linear problem) are filled by external means
        GlobalVector(Base::MeshManipulator* theMesh, int elementVectorID = 0, int faceVectorID = 0);

        ///for post-processing: puts the solution in the time-level data of the elements
        virtual void writeTimeLevelData(int timeLevel, int variable = 0)=0;
        virtual void constructFromTimeLevelData(int timelevel, int variable = 0)=0;

        ///(re-)collects element vectors and boundary information into this vector
        virtual void assemble()=0;

        ///cleans out all the entries, putting them back to 0
        virtual void reset()=0;

    protected:
        
        int meshLevel_, elementVectorID_, faceVectorID_;
        std::vector<int> startPositionsOfElementsInTheVector_;
        std::vector<int> startPositionsOfFacesInTheVector_;
        std::vector<int> startPositionsOfEdgesInTheVector_;
        std::vector<int> startPositionsOfVerticesInTheVector_;
        Base::MeshManipulator *theMesh_;
        
    };
    
#if defined(HPGEM_USE_PETSC) || defined(HPGEM_USE_COMPLEX_PETSC)
    ///\bug this class depends on PETSc and is likely to cause naming conflicts between the c and c++ standard libraries (workaround: make sure to include all other needed hpGEM headers before including this header)
    class GlobalPetscVector : public GlobalVector
    {
        
    public:
        ///for now provides implicit conversion to Vec (the PETSc vector type)
        ///this needs a special function because deriving from Mat appears to be illegal
        ///this class handles the data management itself, please DON'T pass it to functions like MatDestroy or MatCreate
        ///\bug need a better way to provide an interface to the supported Mat routines AND to other routines that need a Mat (like KSPSolve())
        operator Vec();

        GlobalPetscVector(Base::MeshManipulator* theMesh, int elementVectorID = 0, int faceVectorID = 0);
        ~GlobalPetscVector();

        void writeTimeLevelData(int timeLevel, int variable = 0);
        void constructFromTimeLevelData(int timelevel, int variable = 0);

        void reset();

        void assemble();

    private:
        
        std::vector<PetscInt> makePositionsInVector(const Base::Element*);

    private:
        
        Vec b_;
    };
#endif

}

#endif /* GLOBALVECTOR_HPP_ */
