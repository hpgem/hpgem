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
//------------------------------------------------------------------------------
#ifndef HPGEM_KERNEL_GAUSSQUADRATURERULESFORTRIANGULARPRISM_H
#define HPGEM_KERNEL_GAUSSQUADRATURERULESFORTRIANGULARPRISM_H
//---------------------------------------------------------------------------
#include "Integration/QuadratureRules/GaussQuadratureRule.h"
#include "Geometry/PointReference.h"
#include <vector>

namespace hpgem {

//---------------------------------------------------------------------------
namespace QuadratureRules {

//---------------------------------------------------------------------------
class TriPrism_1_1 : public GaussQuadratureRule {
   public:
    static TriPrism_1_1& Instance() {
        static TriPrism_1_1 theInstance;
        return theInstance;
    }
    TriPrism_1_1(const TriPrism_1_1&) = delete;

    std::string getName() const final;
    std::size_t order() const final;
    std::size_t dimension() const final;
    std::size_t getNumberOfPoints() const final;
    double weight(std::size_t i) const final;
    const Geometry::PointReferenceBase& getPoint(std::size_t i) const final;
    Geometry::ReferenceGeometry* forReferenceGeometry() const final;

   private:
    TriPrism_1_1();

    const std::string name_;
    double weight_[1];
    Geometry::ReferenceGeometry* const refGeoPtr_;
    std::vector<Geometry::PointReference<3>> gp_;
};

//---------------------------------------------------------------------------
class TriPrism_5_21 : public GaussQuadratureRule {
   public:
    static TriPrism_5_21& Instance() {
        static TriPrism_5_21 theInstance;
        return theInstance;
    }
    TriPrism_5_21(const TriPrism_5_21&) = delete;

    std::string getName() const final;
    std::size_t order() const final;
    std::size_t dimension() const final;
    std::size_t getNumberOfPoints() const final;
    double weight(std::size_t i) const final;
    const Geometry::PointReferenceBase& getPoint(std::size_t i) const final;
    Geometry::ReferenceGeometry* forReferenceGeometry() const final;

   private:
    TriPrism_5_21();

    const std::string name_;
    double weight_[21];
    Geometry::ReferenceGeometry* const refGeoPtr_;
    std::vector<Geometry::PointReference<3>> gp_;
};

//---------------------------------------------------------------------------
class TriPrism_7_64 : public GaussQuadratureRule {
   public:
    static TriPrism_7_64& Instance() {
        static TriPrism_7_64 theInstance;
        return theInstance;
    }
    TriPrism_7_64(const TriPrism_7_64&) = delete;

    std::string getName() const final;
    std::size_t order() const final;
    std::size_t dimension() const final;
    std::size_t getNumberOfPoints() const final;
    double weight(std::size_t i) const final;
    const Geometry::PointReferenceBase& getPoint(std::size_t i) const final;
    Geometry::ReferenceGeometry* forReferenceGeometry() const final;

   private:
    TriPrism_7_64();

    const std::string name_;
    double weight_[64];
    Geometry::ReferenceGeometry* const refGeoPtr_;
    std::vector<Geometry::PointReference<3>> gp_;
};

//---------------------------------------------------------------------------
}  // namespace QuadratureRules
}  // namespace hpgem

#endif  // HPGEM_KERNEL_GAUSSQUADRATURERULESFORTRIANGULARPRISM_H
