/*
 * GlobalVector.hpp
 *
 *  Created on: Mar 20, 2014
 *      Author: brinkf
 */

#ifndef GLOBALVECTOR_HPP_
#define GLOBALVECTOR_HPP_

#ifdef hpGEM_INCLUDE_PETSC_SUPPORT
#include "petscvec.h"
#endif
#include <vector>
namespace Base{
class MeshManipulator;
class Element;
}

namespace Utilities{

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
		virtual void constructFromTimeLevelData(int timelevel)=0;

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

#ifdef hpGEM_INCLUDE_PETSC_SUPPORT
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
		void constructFromTimeLevelData(int timelevel);

		void reset();

		void assemble();

	private:

		void makePositionsInVector(int,const Base::Element*,int[]);

	private:

		Vec b_;
	};
#endif

}





#endif /* GLOBALVECTOR_HPP_ */
