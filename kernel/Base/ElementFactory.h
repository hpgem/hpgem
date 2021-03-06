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

#ifndef HPGEM_KERNEL_ELEMENTFACTORY_H
#define HPGEM_KERNEL_ELEMENTFACTORY_H

#include <vector>

#include "Element.h"
#include "GlobalUniqueIndex.h"
#include "Zone.h"

namespace hpgem {

namespace Geometry {
template <std::size_t DIM>
class PointPhysical;
}

namespace Base {
class Element;
class BasisFunctionSet;

//! Element constructors need a lot of information that is the same for each
//! element, this information goes here
class ElementFactory {
   public:
    using CollectionOfBasisFunctionSets =
        Element::CollectionOfBasisFunctionSets;
    static ElementFactory& instance() {
        static ElementFactory theInstance;
        return theInstance;
    }

    ElementFactory(const ElementFactory& orig) = delete;
    ElementFactory& operator=(const ElementFactory& other) = delete;

    //! provide the non-constant information and get an Element!
    template <std::size_t DIM>
    Element* makeElement(const std::vector<std::size_t>& globalNodeIndexes,
                         std::vector<Geometry::PointPhysical<DIM> >& points,
                         Zone& zone, std::size_t owner, bool owning);

    //! mesh creation routines can use this to set their desired defaults
    void setCollectionOfBasisFunctionSets(
        const CollectionOfBasisFunctionSets* functions);

    //! mesh creation routines can use this to set their desired defaults
    void setNumberOfUnknowns(std::size_t unknowns);

    //! mesh creation routines can use this to set their desired defaults
    void setNumberOfTimeLevels(std::size_t timeLevels);

    //! mesh creation routines can use this to set their desired defaults
    void setNumberOfMatrices(std::size_t matrices);

    //! mesh creation routines can use this to set their desired defaults
    void setNumberOfVectors(std::size_t vectors);

   private:
    ElementFactory();

    std::size_t unknowns_;
    const CollectionOfBasisFunctionSets* basisFunctionSets_;
    std::size_t timeLevels_;
    std::size_t numberOfElementMatrices_;
    std::size_t numberOfElementVectors_;
};

//! provide the non-constant information and get an Element!
template <std::size_t DIM>
Element* ElementFactory::makeElement(
    const std::vector<std::size_t>& globalNodeIndexes,
    std::vector<Geometry::PointPhysical<DIM> >& points, Zone& zone,
    std::size_t owner, bool owning) {
    return new Element(
        globalNodeIndexes, basisFunctionSets_, points, unknowns_, timeLevels_,
        GlobalUniqueIndex::instance().getElementIndex(), zone, owner, owning,
        numberOfElementMatrices_, numberOfElementVectors_);
}

}  // namespace Base

}  // namespace hpgem

#endif  // HPGEM_KERNEL_ELEMENTFACTORY_H
