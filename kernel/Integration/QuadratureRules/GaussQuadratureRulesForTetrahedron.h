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
#ifndef HPGEM_KERNEL_GAUSSQUADRATURERULESFORTETRAHEDRON_H
#define HPGEM_KERNEL_GAUSSQUADRATURERULESFORTETRAHEDRON_H
//---------------------------------------------------------------------------
#include "Integration/QuadratureRules/GaussQuadratureRule.h"
#include "Geometry/PointReference.h"
#include <vector>
#include "StandardQaussQuadratureRule.h"

namespace hpgem {

//---------------------------------------------------------------------------
namespace QuadratureRules {

//---------------------------------------------------------------------------
class Tn3_1_1 : public GaussQuadratureRule {
   public:
    static Tn3_1_1& Instance() {
        static Tn3_1_1 theInstance;
        return theInstance;
    }
    Tn3_1_1(const Tn3_1_1&) = delete;

    std::string getName() const final;
    std::size_t order() const final;
    std::size_t dimension() const final;
    std::size_t getNumberOfPoints() const final;
    double weight(std::size_t i) const final;
    const Geometry::PointReferenceBase& getPoint(std::size_t i) const final;
    Geometry::ReferenceGeometry* forReferenceGeometry() const final;

   private:
    Tn3_1_1();

    const std::string name_;
    double weight_[1];
    Geometry::ReferenceGeometry* const refGeoPtr_;
    std::vector<Geometry::PointReference<3>> gp_;
};

//---------------------------------------------------------------------------
class Tn3_2_4 : public GaussQuadratureRule {
   public:
    static Tn3_2_4& Instance() {
        static Tn3_2_4 theInstance;
        return theInstance;
    }
    Tn3_2_4(const Tn3_2_4&) = delete;

    std::string getName() const final;
    std::size_t order() const final;
    std::size_t dimension() const final;
    std::size_t getNumberOfPoints() const final;
    double weight(std::size_t i) const final;
    const Geometry::PointReferenceBase& getPoint(std::size_t i) const final;
    Geometry::ReferenceGeometry* forReferenceGeometry() const final;

   private:
    Tn3_2_4();

    const std::string name_;
    double weight_[4];
    Geometry::ReferenceGeometry* const refGeoPtr_;
    std::vector<Geometry::PointReference<3>> gp_;
};

//---------------------------------------------------------------------------
class Tn3_3_10 : public StandardGaussQuadratureRule<3> {
   public:
    static Tn3_3_10& Instance() {
        static Tn3_3_10 theInstance;
        return theInstance;
    }
    Tn3_3_10(const Tn3_3_10&) = delete;

   private:
    Tn3_3_10();
};

//---------------------------------------------------------------------------
class Tn3_4_14 : public StandardGaussQuadratureRule<3> {
   public:
    static Tn3_4_14& Instance() {
        static Tn3_4_14 theInstance;
        return theInstance;
    }

   private:
    Tn3_4_14();
};

//---------------------------------------------------------------------------
class T3_5_14 : public GaussQuadratureRule {
   public:
    static T3_5_14& Instance() {
        static T3_5_14 theInstance;
        return theInstance;
    }
    T3_5_14(const T3_5_14&) = delete;

    std::string getName() const final;
    std::size_t order() const final;
    std::size_t dimension() const final;
    std::size_t getNumberOfPoints() const final;
    double weight(std::size_t i) const final;
    const Geometry::PointReferenceBase& getPoint(std::size_t i) const final;
    Geometry::ReferenceGeometry* forReferenceGeometry() const final;

   private:
    T3_5_14();

    const std::string name_;
    double weight_[14];
    Geometry::ReferenceGeometry* const refGeoPtr_;
    std::vector<Geometry::PointReference<3>> gp_;
};

//---------------------------------------------------------------------------
class T3_6_24 : public GaussQuadratureRule {
   public:
    static T3_6_24& Instance() {
        static T3_6_24 theInstance;
        return theInstance;
    }
    T3_6_24(const T3_6_24&) = delete;

    std::string getName() const final;
    std::size_t order() const final;
    std::size_t dimension() const final;
    std::size_t getNumberOfPoints() const final;
    double weight(std::size_t i) const final;
    const Geometry::PointReferenceBase& getPoint(std::size_t i) const final;
    Geometry::ReferenceGeometry* forReferenceGeometry() const final;

   private:
    T3_6_24();

    const std::string name_;
    double weight_[24];
    Geometry::ReferenceGeometry* const refGeoPtr_;
    std::vector<Geometry::PointReference<3>> gp_;
};

//---------------------------------------------------------------------------
class T3_7_31 : public StandardGaussQuadratureRule<3> {
   public:
    static T3_7_31& Instance() {
        static T3_7_31 theInstance;
        return theInstance;
    }
    T3_7_31(const T3_7_31&) = delete;

   private:
    T3_7_31();
};

//---------------------------------------------------------------------------
class T3_8_43 : public GaussQuadratureRule {
   public:
    static T3_8_43& Instance() {
        static T3_8_43 theInstance;
        return theInstance;
    }
    T3_8_43(const T3_8_43&) = delete;

    std::string getName() const final;
    std::size_t order() const final;
    std::size_t dimension() const final;
    std::size_t getNumberOfPoints() const final;
    double weight(std::size_t i) const final;
    const Geometry::PointReferenceBase& getPoint(std::size_t i) const final;
    Geometry::ReferenceGeometry* forReferenceGeometry() const final;

   private:
    T3_8_43();

    const std::string name_;
    double weight_[43];
    Geometry::ReferenceGeometry* const refGeoPtr_;
    std::vector<Geometry::PointReference<3>> gp_;
};

//---------------------------------------------------------------------------
class T3_9_53 : public StandardGaussQuadratureRule<3> {
   public:
    static T3_9_53& Instance() {
        static T3_9_53 theInstance;
        return theInstance;
    }
    T3_9_53(const T3_9_53&) = delete;

   private:
    T3_9_53();
};

//---------------------------------------------------------------------------
// the class name lies, this is actually an 11th order rule
class T3_10_126 : public GaussQuadratureRule {
   public:
    static T3_10_126& Instance() {
        static T3_10_126 theInstance;
        return theInstance;
    }
    T3_10_126(const T3_10_126&) = delete;

    std::string getName() const final;
    std::size_t order() const final;
    std::size_t dimension() const final;
    std::size_t getNumberOfPoints() const final;
    double weight(std::size_t i) const final;
    const Geometry::PointReferenceBase& getPoint(std::size_t i) const final;
    Geometry::ReferenceGeometry* forReferenceGeometry() const final;

   private:
    T3_10_126();

    const std::string name_;
    double weight_[126];
    Geometry::ReferenceGeometry* const refGeoPtr_;
    std::vector<Geometry::PointReference<3>> gp_;
};

//---------------------------------------------------------------------------
}  // namespace QuadratureRules
}  // namespace hpgem

#endif  // HPGEM_KERNEL_GAUSSQUADRATURERULESFORTETRAHEDRON_H
