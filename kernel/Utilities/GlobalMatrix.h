/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2014, University of Twente
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

 1. Redistributions of source code must retain the above copyright notice, this
 list of conditions and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.

 3. Neither the name of the copyright holder nor the names of its contributors
 may be used to endorse or promote products derived from this software without
 specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

// WARNING: THIS HEADER PARTLY DEPENDS ON C LIBRARIES - INCLUDING THIS HEADER
// BEFORE OR AFTER OTHER HEADERS GREATLY INCLEASES THE RISK OF NAMING CONFLICTS
// OR OTHER COMPILE TIME ERRORS
#ifndef HPGEM_KERNEL_GLOBALMATRIX_H
#define HPGEM_KERNEL_GLOBALMATRIX_H

#if defined(HPGEM_USE_ANY_PETSC)
#include <petscmat.h>
#endif
#include "GlobalIndexing.h"
#include <vector>
#include <map>

namespace Base {
class MeshManipulatorBase;
class Face;
class Element;
}  // namespace Base

namespace Utilities {

/**
 *  General global assembly class for matrices. Defines the general routines
 * that specializations will always need Subclasses should also derive from the
 * external package specific matrix type they try to implement
 */
class GlobalMatrix {

   public:
    /// use the destructor of the subclass in case this is needed
    virtual ~GlobalMatrix() = default;

    GlobalMatrix(const GlobalMatrix& other) = delete;

    /// constructs the global matrix and performs element assembly
    GlobalMatrix(const GlobalIndexing& rowIndexing,
                 const GlobalIndexing& columnIndexing, int elementMatrixID,
                 int faceMatrixID);

    /// \brief Reinitialize the matrix
    ///
    /// Reinitialize the matrix to match the current state of the local
    /// matrices and the GlobalIndex.
    virtual void reinit() = 0;

    /// \brief Assemble the matrix from elements and faces
    virtual void assemble() = 0;

    // Retrieve the global indices of all basis functions owned by the face, its
    // edges and its nodes. Indices are with respect to the rows of the matrix.
    void getMatrixBCEntries(const Base::Face* face,
                            std::size_t& numberOfEntries,
                            std::vector<int>& entries);

   protected:
    int meshLevel_, elementMatrixID_, faceMatrixID_;
    const GlobalIndexing& rowIndexing_;
    const GlobalIndexing& columnIndexing_;
    /// Whether the row and column index are the physically same object
    const bool symmetricIndexing_;
};
#if defined(HPGEM_USE_ANY_PETSC)
///\bug this class depends on PETSc and is likely to cause naming conflicts
/// between the c and c++ standard libraries (workaround: make sure to include
/// all other needed hpGEM headers before including this header)
class GlobalPetscMatrix : public GlobalMatrix {

   public:
    /// for now provides implicit conversion to Mat (the PETSc matrix type)
    /// this needs a special function because deriving from Mat appears to be
    /// illegal \bug need a better way to provide an interface to the supported
    /// Mat routines AND to other routines that need a Mat (like
    /// KSPSetOperators()) (but not stuff like MatDestroy())
    operator Mat();

    GlobalPetscMatrix(const GlobalIndexing& indexing, int elementMatrixID,
                      int faceMatrixID = -1)
        : GlobalPetscMatrix(indexing, indexing, elementMatrixID, faceMatrixID) {
    }

    GlobalPetscMatrix(const GlobalIndexing& rowIndexing,
                      const GlobalIndexing& columnIndexing, int elementMatrixID,
                      int faceMatrixID = -1);

    ~GlobalPetscMatrix() override;

    void assemble() override;

    void reinit() override;

    void printMatInfo(MatInfoType type, std::ostream& stream);
    void writeMatlab(const std::string& fileName);

    const GlobalIndexing& getGlobalIndex() const { return rowIndexing_; }

   private:
    void createMat();

    Mat A_;
};
#endif
}  // namespace Utilities

#endif  // HPGEM_KERNEL_GLOBALMATRIX_H
