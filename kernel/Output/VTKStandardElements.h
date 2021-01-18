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
#ifndef HPGEM_VTKSTANDARDELEMENTS_H
#define HPGEM_VTKSTANDARDELEMENTS_H

#include "VTKElement.h"

namespace hpgem {
namespace Output {

/// Simple hardcoded (linear) elements

class VTKPoint final : public VTKElement<0> {
   public:
    std::uint8_t vtkId() const final { return 1; }

    const std::vector<Geometry::PointReference<0>>& getPoints() const final {
        return points_;
    }

   private:
    static std::vector<Geometry::PointReference<0>> points_;
};

class VTKLine final : public VTKElement<1> {
   public:
    std::uint8_t vtkId() const final { return 3; }

    const std::vector<Geometry::PointReference<1>>& getPoints() const final {
        return points_;
    }

   private:
    static std::vector<Geometry::PointReference<1>> points_;
};

class VTKTriangle final : public VTKElement<2> {
   public:
    std::uint8_t vtkId() const final { return 5; }

    const std::vector<Geometry::PointReference<2>>& getPoints() const final {
        return points_;
    }

   private:
    static std::vector<Geometry::PointReference<2>> points_;
};

class VTKQuad final : public VTKElement<2> {
   public:
    std::uint8_t vtkId() const final { return 9; }

    const std::vector<Geometry::PointReference<2>>& getPoints() const final {
        return points_;
    }

   private:
    static std::vector<Geometry::PointReference<2>> points_;
};

class VTKTetra final : public VTKElement<3> {
   public:
    std::uint8_t vtkId() const final { return 10; }

    const std::vector<Geometry::PointReference<3>>& getPoints() const final {
        return points_;
    }

   private:
    static std::vector<Geometry::PointReference<3>> points_;
};

class VTKHexahedron final : public VTKElement<3> {
   public:
    std::uint8_t vtkId() const final { return 12; }

    const std::vector<Geometry::PointReference<3>>& getPoints() const final {
        return points_;
    }

   private:
    static std::vector<Geometry::PointReference<3>> points_;
};

class VTKWedge final : public VTKElement<3> {
   public:
    std::uint8_t vtkId() const final { return 13; }

    const std::vector<Geometry::PointReference<3>>& getPoints() const final {
        return points_;
    }

   private:
    static std::vector<Geometry::PointReference<3>> points_;
};

class VTKPyramid final : public VTKElement<3> {
   public:
    std::uint8_t vtkId() const final { return 14; }

    const std::vector<Geometry::PointReference<3>>& getPoints() const final {
        return points_;
    }

   private:
    static std::vector<Geometry::PointReference<3>> points_;
};
}  // namespace Output
}  // namespace hpgem

#endif  // HPGEM_VTKSTANDARDELEMENTS_H
