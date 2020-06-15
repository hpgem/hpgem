/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2019, University of Twente
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

#ifndef HPGEM_ELEMENTBASISFUNCTIONS_H
#define HPGEM_ELEMENTBASISFUNCTIONS_H

#include <memory>
#include <vector>

#include "BasisFunctionSet.h"

namespace Base {

/// \brief Association information of BasisFunctions with an Element
///
/// Association between the collection of basisFunctionSet's (with all basis
/// functions used on any element for any unknown) to those used for each
/// unknown. Each unknown can be associated with several basisFunctionSet's
/// from the global one, which results in a linear index of the basis
/// functions for each unknown.
///
/// Note that this order is important for GlobalIndexing, as it expects the
/// specific layout of first the basis functions of the element, than those
/// of the faces, edges and nodes each following their local index. Which is
/// needed for the correct assembly with conforming basis functions.
class ElementBasisFunctions {
   public:
    /// \brief Apply legacy behaviour with respect to systems with multiple
    /// unknowns.
    ///
    /// Before the introduction of different basis functions for different
    /// unknowns, the code expected that each unknown had the same basis
    /// functions and no 'unknown' parameter was present for most functions.
    /// To prevent accidental mistakes one can still use some functions
    /// without unknown argument, which will then check if all unknowns have
    /// the same number of the relevant quantities like the same number of
    /// basis functions. This is done by passing this value as default to
    /// the 'unknown' argument of those functions.
    const static std::size_t LEGACY_BEHAVIOUR =
        std::numeric_limits<std::size_t>::max();
    using CollectionOfBasisFunctionSets =
        std::vector<std::shared_ptr<const BasisFunctionSet>>;

    /// Create an ElementBasisFunction with no basis functions and no unknowns.
    ElementBasisFunctions() : sets_(nullptr), setPositions_() {}

    /// \brief Constructor that sets the collection of basisFunctionSets and
    ///   number of unknowns without assigning basis functions.
    ///
    /// Note unlike the legacy behavior there is no association made between
    /// basis functions in the set and any of the unknowns.
    /// \param sets The set of basisFunctionSet's
    /// \param numberOfUnknowns The number of unknowns.
    ElementBasisFunctions(const CollectionOfBasisFunctionSets *sets,
                          std::size_t numberOfUnknowns)
        : sets_(sets), setPositions_(numberOfUnknowns) {}

    // Ensure all standard functions are present
    ElementBasisFunctions(const ElementBasisFunctions &other) = default;
    ElementBasisFunctions(ElementBasisFunctions &&other) = default;

    ~ElementBasisFunctions() = default;

    ElementBasisFunctions &operator=(const ElementBasisFunctions &other) =
        default;
    ElementBasisFunctions &operator=(ElementBasisFunctions &&other) = default;

    /// \brief The number of basis functions for the unknown that have
    /// support on the element.
    ///
    /// \param unknown The unknown
    /// \return The number of basis functions on the element
    std::size_t getNumberOfBasisFunctions(
        std::size_t unknown = LEGACY_BEHAVIOUR) const;

    /// \brief The number of basis functions with only support on the element
    ///
    /// Count the number of basis functions that have only support on the
    /// element and its boundary but not on other elements.
    /// \param unknown The unknown
    /// \return The amount
    std::size_t getNumberOfLocalBasisFunctions(
        std::size_t unknown = LEGACY_BEHAVIOUR) const;

    /// \brief Total count of basis functions with only support on the element
    ///
    /// \return The amount
    std::size_t getTotalLocalNumberOfBasisFunctions() const;

    /// \brief Maximum polynomial order of the basis functions used on this
    /// element \return The order
    std::size_t getMaximumOrder() const;

    /// \brief Convert from element basis function index to the
    ///   basisFunctionSet to which it belongs and the index in that set.
    ///
    /// Convert the linear index on the basis functions on an element for a
    /// single unknown, to the BasisFunctionSet of which this basis function
    /// is part and the index in this set.
    ///
    /// \param index The index of the basis function for the element
    /// \param unknown The unknown for which to covert this index
    /// \return The basisFunctionSet to which this basis function belongs
    ///   and the index in that set.
    std::pair<const BasisFunctionSet *, std::size_t>
        getBasisFunctionSetAndIndex(
            size_t index, std::size_t unknown = LEGACY_BEHAVIOUR) const;

    /// Validate the internal consistency
    void validatePositions() const;

    /// \brief Remove all basis functions for an unknown.
    void clearBasisFunctionPosition(std::size_t unknown);
    /// \brief Associate a basisFunctionSet with an unknown
    ///
    /// Associate a basisFunctionSet with an unknown. To allow the ordering
    /// of the basisFunctionSet's for a single unknown a 'place' can be
    /// specified. Sets with lower values for place get a lower basis
    /// function index that those with higher place. Only one set can be
    /// registered for each place, subsequent registrations will override
    /// previous ones.
    ///
    /// \param unknown The unknown to associate it with
    /// \param place The place (not necessarily consecutive) in which this
    ///   set should appear. Should be non negative
    /// \param position The position in of the basisFunctionSet in the (global)
    ///   collection of basisFunctionSet.
    void registerBasisFunctionPosition(std::size_t unknown, std::size_t place,
                                       std::size_t position);

    /// \brief Get the offset for the basis functions at a certain place.
    ///
    /// Get the number of basis functions for the unknown that are registered
    /// before 'place'. The basis functions corresponding registered at place
    /// are thus preceded by this amount of basis functions in the local
    /// ordering.
    ///
    /// \param unknown The unknown for the basis functions
    /// \param place The place of the basis functions
    /// \return The offset of the basis functions at place. Undefined if none
    /// are registered.
    std::size_t getBasisFunctionOffset(std::size_t unknown,
                                       std::size_t place) const;

   private:
    std::size_t getNumberOfUnknowns() const { return setPositions_.size(); }

    /// Check if the unknown is valid.
    /// \param unknown The unknown to test
    /// \param allowAbsent Whether to allow -1 for the unknown being absent
    /// (legacy handling)
    void assertValidUnknown(long unknown, bool allowAbsent) const {
        logger.assert_debug(
            unknown < static_cast<long>(getNumberOfUnknowns()) ||
                (allowAbsent && unknown == -1),
            "Invalid unknown %, there are only % unknowns", unknown,
            getNumberOfUnknowns());
    }

    /// \brief Collection of all basisfunctions used on the mesh.
    const CollectionOfBasisFunctionSets *sets_;

    /// \brief Indices of the basis functions used on this element
    ///
    /// Double vector with indices to what BasisFunctionSet's are used for
    /// each unknown. This is a double vector as there may be multiple sets
    /// for a single unknown (e.g. in the case of conforming basis
    /// functions).
    ///
    /// The indices are into sets_, with -1 indicating that the place in
    /// this vector is not in use.
    std::vector<std::vector<int>> setPositions_;
};

// Local implementation for inlining

inline std::size_t ElementBasisFunctions::getNumberOfBasisFunctions(
    std::size_t unknown) const {
    assertValidUnknown(unknown, true);
    if (unknown != LEGACY_BEHAVIOUR) {
        std::size_t total = 0;
        for (int pos : setPositions_[unknown]) {
            if (pos != -1) {
                total += sets_->at(pos)->size();
            }
        }
        return total;
    } else if (getNumberOfUnknowns() != 0) {
        // Legacy
        std::size_t dofs = getNumberOfBasisFunctions(0);
        for (std::size_t i = 1; i < getNumberOfUnknowns(); ++i) {
            logger.assert_debug(dofs == getNumberOfBasisFunctions(i),
                                "Unequal number of basis functions");
        }
        return dofs;
    } else {
        // Seems to be a good default
        logger(WARN, "No unknowns");
        return 0;
    }
}

inline std::size_t ElementBasisFunctions::getNumberOfLocalBasisFunctions(
    std::size_t unknown) const {
    assertValidUnknown(unknown, true);
    if (unknown != LEGACY_BEHAVIOUR) {
        if (setPositions_[unknown].empty() || setPositions_[unknown][0] == -1) {
            return 0;
        } else {
            return sets_->at(setPositions_[unknown][0])->size();
        }

    } else if (getNumberOfUnknowns() != 0) {
        // Legacy behaviour
        std::size_t result = getNumberOfLocalBasisFunctions(0);
        for (std::size_t i = 1; i < getNumberOfUnknowns(); ++i) {
            logger.assert_debug(result == getNumberOfLocalBasisFunctions(i),
                                "Unequal number of local basis functions");
        }
        return result;
    } else {
        logger(WARN, "No unknowns, so assuming no basis functions");
        return 0;
    }
}

inline std::size_t ElementBasisFunctions::getTotalLocalNumberOfBasisFunctions()
    const {
    std::size_t total = 0;
    for (std::size_t i = 0; i < getNumberOfUnknowns(); ++i) {
        total += getNumberOfLocalBasisFunctions(i);
    }
    return total;
}
}  // namespace Base

#endif  // HPGEM_ELEMENTBASISFUNCTIONS_H
