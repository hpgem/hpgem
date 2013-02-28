#ifndef MeshRefiner_hpp
#define MeshRefiner_hpp

#include "TypedefsInBase.hpp"

namespace Base
{
    // to refine an element, refiner requires access to:
    //   - global nodes
    //   - physical nodes 
    //   - list of Elements and levelTree  (two different entities?)
    //   - list of Faces and levelTree  (two different entities?)
    class MeshRefiner
    {
      public:
	  MeshRefiner();
	  
      private:
	  typedef Vector<unsigned int> VectorOfLocalIndices;
	  
	  LocalIndexT	createNewNodes(const ElementIteratorT el, const RefinementT refType);
	  
	  LocalIndexT	createSubElements(const ElementIteratorT el, const RefinementT refType);
	  bool          isRefineThisFace(const FaceIteratorT fa) const;
	  LocalIndexT	createSubFaces(const FaceIteratorT fa);

	  LocalIndexT	subElementsOnFace(const RefinementT refType, const LocalIndexT iFace,
					  VectorOfLocalIndexT&) const;

	  LocalIndexT	getLocalSubFacesNr(const ElementIteratorT el, const int refType,
					   VectorOfLocalIndexT& localSubFacesNr);

	  void pairingCheck(const ElementIteratorT elL, const DimType locFaceNrL, 
			    const ElementIteratorT elR, const DimType locFaceNrR, 
			    int& pairingValue, bool& sizeOrder);

	  bool isPeriodicMatch(const DimType periodicDim, const PhysSpacePoint<dim>& P, 
			       const PhysSpacePoint<dim>& Q, const CoordType tolerance);

	  bool isPeriodicMatch(const FaceIteratorT fa,
			       const ElementIteratorT el1, const DimType localFaceNr1,
			       const ElementIteratorT el2, const DimType localFaceNr2,
			       const CoordType tolerance);

      private:
	  VectorOfNodes     nodes_;	// parent's and childrens' nodes during refining an element
	  VectorOfElements  subElements_;	// sub-elements of an element
	  VectorOfFaces     faces_;		
    };
    
} // close namespace Base
#endif //  MeshRefiner_hpp
