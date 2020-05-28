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

#ifndef MeshRefiner_h
#define MeshRefiner_h

//#include "TypedefsInBase.h"

namespace Base {
// to refine an element, refiner requires access to:
//   - global nodes
//   - physical nodes
//   - list of Elements and levelTree  (two different entities?)
//   - list of Faces and levelTree  (two different entities?)
///\deprecated This should not even compile, please do not try to use it in it's
///current state.
class MeshRefiner {
   public:
    MeshRefiner();

   private:
    using VectorOfLocalIndices = Vector<unsigned int>;

    LocalIndexT createNewNodes(const ElementIteratorT el,
                               const RefinementT refType);

    LocalIndexT createSubElements(const ElementIteratorT el,
                                  const RefinementT refType);
    bool isRefineThisFace(const FaceIteratorT fa) const;
    LocalIndexT createSubFaces(const FaceIteratorT fa);

    LocalIndexT subElementsOnFace(const RefinementT refType,
                                  const LocalIndexT iFace,
                                  VectorOfLocalIndexT&) const;

    LocalIndexT getLocalSubFacesNr(const ElementIteratorT el, const int refType,
                                   VectorOfLocalIndexT& localSubFacesNr);

    void pairingCheck(const ElementIteratorT elL, const DimType locFaceNrL,
                      const ElementIteratorT elR, const DimType locFaceNrR,
                      int& pairingValue, bool& sizeOrder);

    bool isPeriodicMatch(const DimType periodicDim,
                         const PhysSpacePoint<dim>& P,
                         const PhysSpacePoint<dim>& Q,
                         const CoordType tolerance);

    bool isPeriodicMatch(const FaceIteratorT fa, const ElementIteratorT el1,
                         const DimType localFaceNr1, const ElementIteratorT el2,
                         const DimType localFaceNr2, const CoordType tolerance);

   private:
    VectorOfNodes nodes_;  // parent's and childrens' nodes during refining an
                           // element
    VectorOfElements subElements_;  // sub-elements of an element
    VectorOfFaces faces_;
};

}  // namespace Base
#endif  //  MeshRefiner_h
