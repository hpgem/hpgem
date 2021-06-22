/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2021, University of Twente
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
#ifndef HPGEM_ABSTRACTDIMENSIONLESSBASE_H
#define HPGEM_ABSTRACTDIMENSIONLESSBASE_H

#include "Logger.h"

namespace hpgem {

/// \brief Abstraction around dimension erasure
///
/// In hpgem there are many classes that have a fixed  dimension, for example:
/// Points and Quadrature rules. For efficiency reasons this dimension is
/// usually a template parameter reducing the required dynamic allocation.
/// However, the downside of having such a template parameter is that the
/// code becomes more verbose (template parameters everywhere) and slower to
/// compile.
///
/// To reduce the impact of having dimension template parameters on some classes
/// a common pattern is to introduce a base class that is independent of the
/// dimension thus
/// \code
/// class TBase;
/// template<std::size_t> class T : public TBase;
/// \endcode
///
/// This base class then serves two purposes:
///  1. To allow pointing to an instance from a class which is not templated by
///     dimension.
///  2. As place for constants and methods which are (in their signature) not
///     dependent on the dimension template parameter. For example, the number
///     of nodes in a reference geometry, or the L2 norm of a vector.
///
/// The reference from purpose 1 allows calling dimension-template independent
/// methods. However, there are many places where we actually know (from
/// context) the dimension of the actual instance, and need a method that
/// depends in signature on the dimensional template. Thus we need a way to
/// upcast from TBase to T with some known dimension.
///
/// To simplify this common upcast we use that this conversion is only used for
/// some low dimensions (0-4). To simplify this we allow for implicit conversion
/// from TBase to T for dimension 0-4. (and const variants). Thus allowing us to
/// write \code T<3>& t = object->getTBase(); \endcode (where we know from
/// context that it is the 3 dimensional instance).
///
/// However writing these conversions is rather tedious, on needs to implement
/// 10 methods (dimensions 0-4, both const and non-const) for every base class.
/// This abstraction class removes this writing of boilerplate code. It should
/// be used as
/// \code
/// template<std::size_t DIM>
/// class T;
///
/// class TBase : public AbstractDimensionlessBase<TBase, T> {};
///
/// template<std::size_t>
/// class T : TBase {};
/// \endcode
///
/// \tparam BaseType The base type which does not have a dimension template.
/// \tparam TypeWithDim The actual type with a dimension template.
template <typename BaseType, template <std::size_t> class TypeWithDim>
class AbstractDimensionlessBase {
   public:
    // Implicit conversion operator to allow easy conversion from TBase to
    // T<0> - T<4>. As this implicit conversion is intended we need to silence
    // clang-tidy warnings about the implicit conversions.

    // NOLINTNEXTLINE(google-explicit-constructor)
    operator TypeWithDim<0>&() { return castDimension<0>(); }
    // NOLINTNEXTLINE(google-explicit-constructor)
    operator TypeWithDim<1>&() { return castDimension<1>(); }
    // NOLINTNEXTLINE(google-explicit-constructor)
    operator TypeWithDim<2>&() { return castDimension<2>(); }
    // NOLINTNEXTLINE(google-explicit-constructor)
    operator TypeWithDim<3>&() { return castDimension<3>(); }
    // NOLINTNEXTLINE(google-explicit-constructor)
    operator TypeWithDim<4>&() { return castDimension<4>(); }

    // NOLINTNEXTLINE(google-explicit-constructor)
    operator const TypeWithDim<0>&() const { return castDimension<0>(); }
    // NOLINTNEXTLINE(google-explicit-constructor)
    operator const TypeWithDim<1>&() const { return castDimension<1>(); }
    // NOLINTNEXTLINE(google-explicit-constructor)
    operator const TypeWithDim<2>&() const { return castDimension<2>(); }
    // NOLINTNEXTLINE(google-explicit-constructor)
    operator const TypeWithDim<3>&() const { return castDimension<3>(); }
    // NOLINTNEXTLINE(google-explicit-constructor)
    operator const TypeWithDim<4>&() const { return castDimension<4>(); }

    /// \brief Fast conversion to implementation type with dimension.
    ///
    /// Cast the current dimensionless instance to the instance with the right
    /// dimension. Casting using the wrong dimension will result in an assertion
    /// failure when hpgem is compiled with asserts, otherwise it results in
    /// undefined behaviour.
    ///
    /// \tparam DIM The actual dimension
    /// \return Reference to this entity with the given dimension
    template <std::size_t DIM>
    inline TypeWithDim<DIM>& castDimension() {
        checkType<DIM>();
        // This static cast is allowed because the template instantiation is
        // done in TBase, hence we know that this is actually of type BaseType*.
        BaseType* baseThis = static_cast<BaseType*>(this);
        // Now we use a static cast for the conversion. The validity is checked
        // (when asserts are enabled) by checkType<>();
        return static_cast<TypeWithDim<DIM>&>(*baseThis);
    }

    template <std::size_t DIM>
    inline const TypeWithDim<DIM>& castDimension() const {
        checkType<DIM>();
        return static_cast<const TypeWithDim<DIM>&>(
            *static_cast<const BaseType*>(this));
    }

   private:
    /// When asserts are enabled: Check that casting to TypeWithDim<DIM> is
    /// valid. Without assertions this check is not done and thus a no-op.
    ///
    /// \tparam DIM The dimension of the target.
    template <std::size_t DIM>
    inline void checkType() const {
#if HPGEM_ASSERTS
        // Use dynamic cast to allow deriving from TypeWithDim.
        const TypeWithDim<DIM>* test = dynamic_cast<const TypeWithDim<DIM>*>(
            static_cast<const BaseType*>(this));
        logger.assert_debug(test != nullptr,
                            "Invalid conversion to dimension %", DIM);
#endif
    }
};

}  // namespace hpgem

#endif  // HPGEM_ABSTRACTDIMENSIONLESSBASE_H
