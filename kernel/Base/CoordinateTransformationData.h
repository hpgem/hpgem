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
#ifndef HPGEM_COORDINATETRANSFORMATIONDATA_H
#define HPGEM_COORDINATETRANSFORMATIONDATA_H

#include "Geometry/Jacobian.h"

namespace hpgem {
namespace Base {

/**
 * Interface for CoordinateTransformation providing the information about the
 * underlying coordinate transformation.
 *
 * Note CoordinateTransformation is a misnomer, it is a transformation of
 * function values based on a coordinate transformation.
 *
 * @tparam DIM The dimension of the coordinate change.
 */
template <std::size_t DIM>
class CoordinateTransformationData {
    // Design note: This class is used to transform function values. The usual
    // application is to transform multiple functions in a row. For example
    // transforming all the basis functions at a fixed quadrature point. Hence,
    // the values of this function will be requested multiple times.
    //
    // The functions are thus designed to give a const reference, as the results
    // should be cached by the implementing class.

   public:
    /**
     * Get the Jacobian of the coordinate transformation
     */
    virtual const Geometry::Jacobian<DIM, DIM>& getJacobian() const = 0;
    /**
     * Get the transpose of the Jacobian.
     */
    virtual const Geometry::Jacobian<DIM, DIM>& getTransposeJacobian()
        const = 0;

    /**
     * Get the determinant of the Jacobian.
     */
    virtual double getJacobianDet() const = 0;
};

/**
 * Simple implementation of CoordinateTransformationData that holds the required
 * values.
 */
template <std::size_t DIM>
class ValueCoordinateTransformationData final
    : public CoordinateTransformationData<DIM> {
   public:
    ValueCoordinateTransformationData()
        : jacobian(), transposeJacobian(), jacobianDet(0.0){};

    const Geometry::Jacobian<DIM, DIM>& getJacobian() const final {
        return jacobian;
    }
    const Geometry::Jacobian<DIM, DIM>& getTransposeJacobian() const final {
        return transposeJacobian;
    }
    double getJacobianDet() const final { return jacobianDet; }

    void setJacobian(const Geometry::Jacobian<DIM, DIM>& newJacobian) {
        jacobian = newJacobian;
        transposeJacobian = jacobian.transpose();
        jacobianDet = jacobian.determinant();
    }

   private:
    Geometry::Jacobian<DIM, DIM> jacobian, transposeJacobian;
    double jacobianDet;
};

}  // namespace Base
}  // namespace hpgem

#endif  // HPGEM_COORDINATETRANSFORMATIONDATA_H
