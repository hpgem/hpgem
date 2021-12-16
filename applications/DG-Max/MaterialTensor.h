/*
This file forms part of hpGEM. This package has been developed over a number of
years by various people at the University of Twente and a full list of
contributors can be found at http://hpgem.org/about-the-code/team

This code is distributed using BSD 3-Clause License. A copy of which can found
below.


Copyright (c) 2021, University of Twente
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef HPGEM_MATERIALCOEFFICIENT_H
#define HPGEM_MATERIALCOEFFICIENT_H

#include <complex>
#include <limits>
#include "LinearAlgebra/SmallVector.h"

namespace DGMax {

class MaterialTensor {
   public:
    using DiagonalTensor = hpgem::LinearAlgebra::SmallVectorC<3>;
    template <std::size_t d>
    using VecC = hpgem::LinearAlgebra::SmallVectorC<d>;

   private:
    enum class Type { SCALAR, DIAGONAL_TENSOR } type;

    union Value {
        std::complex<double> scalar;
        DiagonalTensor diag;

        Value() : scalar(std::numeric_limits<double>::signaling_NaN()){};
        Value(std::complex<double> scalar) : scalar(scalar){};
        Value(const DiagonalTensor& tensor) : diag(tensor){};

        Value(const Value& other, Type type) { set(other, type); }

        void set(const Value& other, Type type) {
            switch (type) {
                case Type::SCALAR:
                    scalar = other.scalar;
                    break;
                case Type::DIAGONAL_TENSOR:
                    diag = other.diag;
                    break;
                default:
                    hpgem::logger.fail("Unknown material tensor type");
            }
        }
    } value;

   public:
    MaterialTensor() = default;
    MaterialTensor(const MaterialTensor& other)
        : type(other.type), value(other.value, other.type) {}

    MaterialTensor& operator=(const MaterialTensor& other) {
        type = other.type;
        value.set(other.value, type);
        return *this;
    }

    explicit MaterialTensor(std::complex<double> val) : type(Type::SCALAR), value(val){};

    explicit MaterialTensor(const DiagonalTensor& val)
        : type(Type::DIAGONAL_TENSOR), value(val){};

    MaterialTensor adjoint() const {
        switch (type) {
            case Type::SCALAR:
                return MaterialTensor(std::conj(value.scalar));
            case Type::DIAGONAL_TENSOR:
                return MaterialTensor(value.diag.conj());
            default:
                hpgem::logger.fail("Unknown material tensor type");
        }
    }

    template<std::size_t d>
    VecC<d> applyDiv(const hpgem::LinearAlgebra::SmallVector<d>& vec) const {
        return applyDiv(VecC<d>(vec));
    }

    /**
     * Apply the material tensor to regular vector like quantities, similar to
     * how it is used in the divergence constraint.
     * @tparam d The dimension of the vector
     * @param vec The vector to apply it to
     * @return The resulting vector
     */
    template <std::size_t d>
    VecC<d> applyDiv(VecC<d> vec) const {
        switch (type) {
            case Type::SCALAR:
                vec *= value.scalar;
                break;
            case Type::DIAGONAL_TENSOR:
                for (std::size_t i = 0; i < d; ++i) {
                    vec[i] *= value.diag[i];
                }
                break;
            default:
                hpgem::logger.fail("Unknown material tensor type");
        }
        return vec;
    }
    /**
     * Apply the material tensor to curl-vector quantities, similar how it is
     * used in 'Curl muinv Curl E'.
     *
     * This will include the inverse that is needed in these cases.
     *
     * Note for 2D this applies the z-component of the material tensor.
     * @param vec The vector to apply it to
     * @return The resulting vector
     */
    VecC<2> applyCurl(VecC<2> vec) const {
        switch (type) {
            case Type::SCALAR:
                vec[0] /= value.scalar;
                break;
            case Type::DIAGONAL_TENSOR:
                vec[0] /= value.diag[2];
                break;
            default:
                hpgem::logger.fail("Unknown material tensor type");
        }
        return vec;
    }
    VecC<3> applyCurl(VecC<3> vec) const {
        switch (type) {
            case Type::SCALAR:
                vec /= value.scalar;
                break;
            case Type::DIAGONAL_TENSOR:
                for (std::size_t i = 0; i < 3; ++i) {
                    vec[i] /= value.diag[i];
                }
                break;
            default:
                hpgem::logger.fail("Unknown material tensor type");
        }
        return vec;
    }
};

}  // namespace DGMax

#endif  // HPGEM_MATERIALCOEFFICIENT_H
