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
#include "VTKStandardElements.h"

namespace hpgem {
namespace Output {

// See also http://www.vtk.org/VTK/img/file-formats.pdf
// for graphical versions

std::vector<Geometry::PointReference<0>> VTKPoint::points_(
    {Geometry::PointReference<0>()});

// Same as reference
std::vector<Geometry::PointReference<1>> VTKLine::points_(
    {Geometry::PointReference<1>({-1.0}), Geometry::PointReference<1>({1.0})});

// Same as reference
std::vector<Geometry::PointReference<2>> VTKTriangle::points_(
    {Geometry::PointReference<2>({0.0, 0.0}),
     Geometry::PointReference<2>({1.0, 0.0}),
     Geometry::PointReference<2>({0.0, 1.0})});

// Swapped points 2/3
std::vector<Geometry::PointReference<2>> VTKQuad::points_(
    {Geometry::PointReference<2>({-1.0, -1.0}),
     Geometry::PointReference<2>({+1.0, -1.0}),
     Geometry::PointReference<2>({+1.0, +1.0}),
     Geometry::PointReference<2>({-1.0, +1.0})});

// Same as reference
std::vector<Geometry::PointReference<3>> VTKTetra::points_(
    {Geometry::PointReference<3>({0.0, 0.0, 0.0}),
     Geometry::PointReference<3>({1.0, 0.0, 0.0}),
     Geometry::PointReference<3>({0.0, 1.0, 0.0}),
     Geometry::PointReference<3>({0.0, 0.0, 1.0})});

// Swapped points 2/3 and 6/7
std::vector<Geometry::PointReference<3>> VTKHexahedron::points_(
    {Geometry::PointReference<3>({-1.0, -1.0, -1.0}),
     Geometry::PointReference<3>({+1.0, -1.0, -1.0}),
     Geometry::PointReference<3>({+1.0, +1.0, -1.0}),
     Geometry::PointReference<3>({-1.0, +1.0, -1.0}),
     Geometry::PointReference<3>({-1.0, -1.0, +1.0}),
     Geometry::PointReference<3>({+1.0, -1.0, +1.0}),
     Geometry::PointReference<3>({+1.0, +1.0, +1.0}),
     Geometry::PointReference<3>({-1.0, +1.0, +1.0})});

// Same as reference
std::vector<Geometry::PointReference<3>> VTKWedge::points_(
    {Geometry::PointReference<3>({0.0, 0.0, -1.0}),
     Geometry::PointReference<3>({1.0, 0.0, -1.0}),
     Geometry::PointReference<3>({0.0, 1.0, -1.0}),
     Geometry::PointReference<3>({0.0, 0.0, +1.0}),
     Geometry::PointReference<3>({1.0, 0.0, +1.0}),
     Geometry::PointReference<3>({0.0, 1.0, +1.0})});

// Different order, top is last instead of first point
std::vector<Geometry::PointReference<3>> VTKPyramid::points_(
    {Geometry::PointReference<3>({-1.0, -1.0, 0.0}),
     Geometry::PointReference<3>({+1.0, -1.0, 0.0}),
     Geometry::PointReference<3>({+1.0, +1.0, 0.0}),
     Geometry::PointReference<3>({-1.0, +1.0, 0.0}),
     Geometry::PointReference<3>({0.0, 0.0, 1.0})});

}  // namespace Output
}  // namespace hpgem
