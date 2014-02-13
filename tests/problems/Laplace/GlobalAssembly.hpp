/*
 * GlobalAssembly.hpp
 *
 *  Created on: Jan 17, 2014
 *      Author: brinkf
 */

#ifndef GLOBALASSEMBLY_HPP_
#define GLOBALASSEMBLY_HPP_

///\TODO split into multiple files & refactor into proper folder

#include "petscmat.h"
#include "Base/MeshManipulator.hpp"

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
		GlobalMatrix(Base::MeshManipulator* theMesh,unsigned int elementMatrixID,unsigned int faceMatrixID);

		///signals the matrix that (some) element matrixes have changed and that the matrix can be assembled again
		///also used after mesh refinement to make sure the global matrix reflects the changes in the mesh
		virtual void reAssemble()=0;

		///cleans out all the entries, putting them back to 0 - keeps the non-zero strucure of the matrix however.
		///After this it assembles the matrix to put it back in a legal state
		virtual void reset()=0;

	protected:

		int meshLevel_,elementMatrixID_,faceMatrixID_;
		std::vector<int> startPositionsOfElementsInTheMatrix_;
		std::vector<int> startPositionsOfFacesInTheMatrix_;
		std::vector<int> startPositionsOfEdgesInTheMatrix_;
		std::vector<int> startPositionsOfVerticesInTheMatrix_;
		Base::MeshManipulator *theMesh_;

	};

	class GlobalPetscMatrix:public GlobalMatrix{

	public:
		///for now provides implicit conversion to Mat (the PETSc matrix type)
		///this needs a special function because deriving from Mat appears to be illegal
		///this class handles the data management itself, please DON'T pass it to functions like MatDestroy or MatCreate
		///\bug need a better way to provide an interface to the supported Mat routines AND to other routines that need a Mat (like KSPSetOperators())
		operator Mat();

		GlobalPetscMatrix(Base::MeshManipulator* theMesh,unsigned int elementMatrixID,unsigned int faceMatrixID);
		~GlobalPetscMatrix();

		void reset();

		void reAssemble();

	private:

		Mat A_;
	};

	/**
	 *  General global assembly class for vectors. Defines the general routines that specializations will always need
	 *  Subclasses should also derive from the external package specific vector type they try to implement
	 */
	class GlobalVector{

	public:
		///use the deconstructor of the subclass in case this is needed
		virtual ~GlobalVector(){}

		///construct the global vector does not do assembly by default
		///because some vectors (like the solution of the linear problem) are filled by external means
		GlobalVector(Base::MeshManipulator* theMesh,int elementVectorID=0,int faceVectorID=0);

		///for post-processing: puts the solution in the time-level data of the elements
		virtual void writeTimeLevelData(int timeLevel)=0;

		///(re-)collects element vectors and boundary information into this vector
		virtual void assemble()=0;

		///cleans out all the entries, putting them back to 0
		virtual void reset()=0;

	protected:

		int meshLevel_,elementVectorID_,faceVectorID_;
		std::vector<int> startPositionsOfElementsInTheVector_;
		std::vector<int> startPositionsOfFacesInTheVector_;
		std::vector<int> startPositionsOfEdgesInTheVector_;
		std::vector<int> startPositionsOfVerticesInTheVector_;
		Base::MeshManipulator *theMesh_;

	};

	class GlobalPetscVector:public GlobalVector{

	public:
		///for now provides implicit conversion to Vec (the PETSc vector type)
		///this needs a special function because deriving from Mat appears to be illegal
		///this class handles the data management itself, please DON'T pass it to functions like MatDestroy or MatCreate
		///\bug need a better way to provide an interface to the supported Mat routines AND to other routines that need a Mat (like KSPSolve())
		operator Vec();

		GlobalPetscVector(Base::MeshManipulator* theMesh,int elementVectorID=0,int faceVectorID=0);
		~GlobalPetscVector();

		void writeTimeLevelData(int timeLevel);

		void reset();

		void assemble();

	private:

		Vec b_;
	};

};


#endif //GLOBALASSEMBLY_HPP_
