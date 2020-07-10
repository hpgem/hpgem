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
#ifndef HPGEM_KERNEL_GAUSSQUADRATURERULESFORTRIANGLE_H
#define HPGEM_KERNEL_GAUSSQUADRATURERULESFORTRIANGLE_H
//---------------------------------------------------------------------------
#include "Integration/QuadratureRules/GaussQuadratureRule.h"
#include "Geometry/PointReference.h"
#include <vector>

//---------------------------------------------------------------------------
namespace QuadratureRules {

//---------------------------------------------------------------------------
class Tn2_1_1 : public GaussQuadratureRule {
   public:
    static Tn2_1_1& Instance() {
        static Tn2_1_1 theInstance;
        return theInstance;
    }
    Tn2_1_1(const Tn2_1_1&) = delete;

    std::string getName() const final;
    std::size_t order() const final;
    std::size_t dimension() const final;
    std::size_t getNumberOfPoints() const final;
    double weight(std::size_t i) const final;
    const Geometry::PointReferenceBase& getPoint(std::size_t i) const final;
    Geometry::ReferenceGeometry* forReferenceGeometry() const final;

   private:
    Tn2_1_1();

    const std::string name_;
    double weight_[1];
    Geometry::ReferenceGeometry* const refGeoPtr_;
    std::vector<Geometry::PointReference<2>> gp_;
};

//---------------------------------------------------------------------------
class Tn2_2_3 : public GaussQuadratureRule {
   public:
    static Tn2_2_3& Instance() {
        static Tn2_2_3 theInstance;
        return theInstance;
    }
    Tn2_2_3(const Tn2_2_3&) = delete;

    std::string getName() const final;
    std::size_t order() const final;
    std::size_t dimension() const final;
    std::size_t getNumberOfPoints() const final;
    double weight(std::size_t i) const final;
    const Geometry::PointReferenceBase& getPoint(std::size_t i) const final;
    Geometry::ReferenceGeometry* forReferenceGeometry() const final;

   private:
    Tn2_2_3();

    const std::string name_;
    double weight_[3];
    Geometry::ReferenceGeometry* const refGeoPtr_;
    std::vector<Geometry::PointReference<2>> gp_;
};

//---------------------------------------------------------------------------
class Tn2_3_4 : public GaussQuadratureRule {
   public:
    static Tn2_3_4& Instance() {
        static Tn2_3_4 theInstance;
        return theInstance;
    }
    Tn2_3_4(const Tn2_3_4&) = delete;

    std::string getName() const final;
    std::size_t order() const final;
    std::size_t dimension() const final;
    std::size_t getNumberOfPoints() const final;
    double weight(std::size_t i) const final;
    const Geometry::PointReferenceBase& getPoint(std::size_t i) const final;
    Geometry::ReferenceGeometry* forReferenceGeometry() const final;

   private:
    Tn2_3_4();

    const std::string name_;
    double weight_[4];
    Geometry::ReferenceGeometry* const refGeoPtr_;
    std::vector<Geometry::PointReference<2>> gp_;
};

//---------------------------------------------------------------------------
class Tn2_4_6 : public GaussQuadratureRule {
   public:
    static Tn2_4_6& Instance() {
        static Tn2_4_6 theInstance;
        return theInstance;
    }
    Tn2_4_6(const Tn2_4_6&) = delete;

    std::string getName() const final;
    std::size_t order() const final;
    std::size_t dimension() const final;
    std::size_t getNumberOfPoints() const final;
    double weight(std::size_t i) const final;
    const Geometry::PointReferenceBase& getPoint(std::size_t i) const final;
    Geometry::ReferenceGeometry* forReferenceGeometry() const final;

   private:
    Tn2_4_6();

    const std::string name_;
    double weight_[6];
    Geometry::ReferenceGeometry* const refGeoPtr_;
    std::vector<Geometry::PointReference<2>> gp_;
};

//---------------------------------------------------------------------------
class T2_5_7 : public GaussQuadratureRule {
   public:
    static T2_5_7& Instance() {
        static T2_5_7 theInstance;
        return theInstance;
    }
    T2_5_7(const T2_5_7&) = delete;

    std::string getName() const final;
    std::size_t order() const final;
    std::size_t dimension() const final;
    std::size_t getNumberOfPoints() const final;
    double weight(std::size_t i) const final;
    const Geometry::PointReferenceBase& getPoint(std::size_t i) const final;
    Geometry::ReferenceGeometry* forReferenceGeometry() const final;

   private:
    T2_5_7();

    const std::string name_;
    double weight_[7];
    Geometry::ReferenceGeometry* const refGeoPtr_;
    std::vector<Geometry::PointReference<2>> gp_;
};

//---------------------------------------------------------------------------
class T2_6_12 : public GaussQuadratureRule {
   public:
    static T2_6_12& Instance() {
        static T2_6_12 theInstance;
        return theInstance;
    }
    T2_6_12(const T2_6_12&) = delete;

    std::string getName() const final;
    std::size_t order() const final;
    std::size_t dimension() const final;
    std::size_t getNumberOfPoints() const final;
    double weight(std::size_t i) const final;
    const Geometry::PointReferenceBase& getPoint(std::size_t i) const final;
    Geometry::ReferenceGeometry* forReferenceGeometry() const final;

   private:
    T2_6_12();

    const std::string name_;
    double weight_[12];
    Geometry::ReferenceGeometry* const refGeoPtr_;
    std::vector<Geometry::PointReference<2>> gp_;
};

//---------------------------------------------------------------------------
class T2_7_13 : public GaussQuadratureRule {
   public:
    static T2_7_13& Instance() {
        static T2_7_13 theInstance;
        return theInstance;
    }
    T2_7_13(const T2_7_13&) = delete;

    std::string getName() const final;
    std::size_t order() const final;
    std::size_t dimension() const final;
    std::size_t getNumberOfPoints() const final;
    double weight(std::size_t i) const final;
    const Geometry::PointReferenceBase& getPoint(std::size_t i) const final;
    Geometry::ReferenceGeometry* forReferenceGeometry() const final;

   private:
    T2_7_13();

    const std::string name_;
    double weight_[13];
    Geometry::ReferenceGeometry* const refGeoPtr_;
    std::vector<Geometry::PointReference<2>> gp_;
};

//---------------------------------------------------------------------------
class T2_8_16 : public GaussQuadratureRule {
   public:
    static T2_8_16& Instance() {
        static T2_8_16 theInstance;
        return theInstance;
    }
    T2_8_16(const T2_8_16&) = delete;

    std::string getName() const final;
    std::size_t order() const final;
    std::size_t dimension() const final;
    std::size_t getNumberOfPoints() const final;
    double weight(std::size_t i) const final;
    const Geometry::PointReferenceBase& getPoint(std::size_t i) const final;
    Geometry::ReferenceGeometry* forReferenceGeometry() const final;

   private:
    T2_8_16();

    const std::string name_;
    double weight_[16];
    Geometry::ReferenceGeometry* const refGeoPtr_;
    std::vector<Geometry::PointReference<2>> gp_;
};

//---------------------------------------------------------------------------
class T2_9_19 : public GaussQuadratureRule {
   public:
    static T2_9_19& Instance() {
        static T2_9_19 theInstance;
        return theInstance;
    }
    T2_9_19(const T2_9_19&) = delete;

    std::string getName() const final;
    std::size_t order() const final;
    std::size_t dimension() const final;
    std::size_t getNumberOfPoints() const final;
    double weight(std::size_t i) const final;
    const Geometry::PointReferenceBase& getPoint(std::size_t i) const final;
    Geometry::ReferenceGeometry* forReferenceGeometry() const final;

   private:
    T2_9_19();

    const std::string name_;
    double weight_[19];
    Geometry::ReferenceGeometry* const refGeoPtr_;
    std::vector<Geometry::PointReference<2>> gp_;
};

//---------------------------------------------------------------------------
class T2_10_25 : public GaussQuadratureRule {
   public:
    static T2_10_25& Instance() {
        static T2_10_25 theInstance;
        return theInstance;
    }
    T2_10_25(const T2_10_25&) = delete;

    std::string getName() const final;
    std::size_t order() const final;
    std::size_t dimension() const final;
    std::size_t getNumberOfPoints() const final;
    double weight(std::size_t i) const final;
    const Geometry::PointReferenceBase& getPoint(std::size_t i) const final;
    Geometry::ReferenceGeometry* forReferenceGeometry() const final;

   private:
    T2_10_25();

    const std::string name_;
    double weight_[25];
    Geometry::ReferenceGeometry* const refGeoPtr_;
    std::vector<Geometry::PointReference<2>> gp_;
};

//---------------------------------------------------------------------------
class T2_11_28 : public GaussQuadratureRule {
   public:
    static T2_11_28& Instance() {
        static T2_11_28 theInstance;
        return theInstance;
    }
    T2_11_28(const T2_11_28&) = delete;

    std::string getName() const final;
    std::size_t order() const final;
    std::size_t dimension() const final;
    std::size_t getNumberOfPoints() const final;
    double weight(std::size_t i) const final;
    const Geometry::PointReferenceBase& getPoint(std::size_t i) const final;
    Geometry::ReferenceGeometry* forReferenceGeometry() const final;

   private:
    T2_11_28();

    const std::string name_;
    double weight_[28];
    Geometry::ReferenceGeometry* const refGeoPtr_;
    std::vector<Geometry::PointReference<2>> gp_;
};

//---------------------------------------------------------------------------
}  // namespace QuadratureRules
#endif  // HPGEM_KERNEL_GAUSSQUADRATURERULESFORTRIANGLE_H
