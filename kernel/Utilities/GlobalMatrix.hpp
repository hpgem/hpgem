/*
 * GlobalMatrix.hpp
 *
 *  Created on: Mar 20, 2014
 *      Author: brinkf
 */

#ifndef GLOBALMATRIX_HPP_
#define GLOBALMATRIX_HPP_

#ifdef hpGEM_INCLUDE_PETSC_SUPPORT
#include "petscmat.h"
#endif
#include <vector>

namespace Base{
class MeshManipulator;
class Face;
class Element;
}

namespace Utilities{

	/**
	 *  General global assembly class for matrices. Defines the general routines that specializations will always need
	 *  Subclasses should also derive from the external package specific matrix type they try to implement
	 */
	class GlobalMatrix{

	public:

		///use the deconstructor of the subclass in case this is needed
		virtual ~GlobalMatrix(){}

		///constructs the global matrix and performs element assembly
		GlobalMatrix(Base::MeshManipulator* theMesh, int elementMatrixID, int faceMatrixID);

		///signals the matrix that (some) element matrixes have changed and that the matrix can be assembled again
		///also used after mesh refinement to make sure the global matrix reflects the changes in the mesh
		virtual void reAssemble()=0;

		///cleans out all the entries, putting them back to 0 - keeps the non-zero strucure of the matrix however.
		///After this it assembles the matrix to put it back in a legal state
		virtual void reset()=0;

		void getMatrixBCEntries(Base::Face* face, int& numberOfEntries, std::vector<int>& entries);

	protected:

		int meshLevel_,elementMatrixID_,faceMatrixID_;
		std::vector<int> startPositionsOfElementsInTheMatrix_;
		std::vector<int> startPositionsOfFacesInTheMatrix_;
		std::vector<int> startPositionsOfEdgesInTheMatrix_;
		std::vector<int> startPositionsOfVerticesInTheMatrix_;
		Base::MeshManipulator *theMesh_;

	};
#ifdef hpGEM_INCLUDE_PETSC_SUPPORT
	class GlobalPetscMatrix:public GlobalMatrix{

	public:
		///for now provides implicit conversion to Mat (the PETSc matrix type)
		///this needs a special function because deriving from Mat appears to be illegal
		///\bug need a better way to provide an interface to the supported Mat routines AND to other routines that need a Mat (like KSPSetOperators()) (but not stuff like MatDestroy())
		operator Mat();

		GlobalPetscMatrix(Base::MeshManipulator* theMesh, int elementMatrixID, int faceMatrixID=-1);
		~GlobalPetscMatrix();

		void reset();

		void reAssemble();

	private:

		void makePositionsInMatrix(int,const Base::Element*,int[]);

	private:

		Mat A_;
	};
#endif
}


#endif /* GLOBALMATRIX_HPP_ */
