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
#include "PredefinedStructure.h"
#include "Base/Element.h"

namespace DGMax {

ElementInfos* PredefinedStructureDescription::createElementInfo(
    const Base::Element* element) {
    switch (dimension_) {
        case 2:
            return createElementInfoDim<2>(element);
        case 3:
            return createElementInfoDim<3>(element);
        default:
            logger.assert_always(false, "Not implemented for dimension %",
                                 dimension_);
            return nullptr;
    }
}

template <std::size_t DIM>
ElementInfos* PredefinedStructureDescription::createElementInfoDim(
    const Base::Element* element) const {
    logger.assert_debug(element->getReferenceGeometry()->getDimension() == DIM,
                        "Incorrect dimension");
    const Geometry::PointReference<DIM>& p =
        element->getReferenceGeometry()->getCenter();
    Geometry::PointPhysical<DIM> pPhys = element->referenceToPhysical(p);

    return new ElementInfos(jelmerStructure(pPhys, structure_));
}

PredefinedStructure structureFromInt(std::size_t value) {
    switch (value) {
        case 0:
            return PredefinedStructure::VACUUM;
        case 1:
            return PredefinedStructure::BRAGG_STACK;
        case 2:
            return PredefinedStructure::CYLINDER;
        case 3:
            return PredefinedStructure::SQUARE_HOLE;
        case 4:
            return PredefinedStructure::INVERSE_WOODPILE_OLD;
        case 5:
            return PredefinedStructure::INVERSE_WOODPILE_NEW;
        default:
            logger.assert_always(false, "Undefined structure %", value);
            return PredefinedStructure::VACUUM;  // To keep the compiler happy
    }
}

/**
 * Check if a coordinate is inside the pores in one of the directions of an
 * inverse woodpile structure. The unit-cell of the pore layout is a
 * rectangle of size c (xz-coordinate) by a (y-coordinate), with one corner at
 * the origin. The pores are placed at the four corners and at the center.
 *
 * For a 3D structure one needs to give the pores in the z direction an offset
 * of a/4 compared to the pores in x direction.
 * @param xz The coordinate in the x or z plane
 * @param y The coordinate in the y direction
 * @param a The lattice constant in the y-direction
 * @param c The lattice constant in the x-direction
 * @param r The radius of the pores
 * @return Whether the coordinate falls in th pore.
 */
bool insideIWPore(double xz, double y, double a, double c, double r) {
    // Code is written with the assumption that only the pores at the four
    // corners and the centre pore are overlapping with the unit cell. This
    // gives an upper bound for the radius.
    logger.assert_debug(r <= a / 2, "Too large a");
    logger.assert_debug(r <= (c / 2), "Too large c");
    // Compute coordinate within a unit cell
    xz -= std::floor(xz / c) * c;
    y -= std::floor(y / a) * a;

    double r2 = r * r;

    // Try corner pores
    double y2 = y * y;
    double xz2 = xz * xz;
    double y2c = (a - y) * (a - y);
    double xz2c = (c - xz) * (c - xz);
    if (y2 + xz2 < r2 || y2 + xz2c < r2 || y2c + xz2 < r2 || y2c + xz2c < r2) {
        return true;
    }
    // Try the center pore
    double ymid = y - a / 2;
    double xzmid = xz - c / 2;
    return ymid * ymid + xzmid * xzmid < r2;
}

// Actually programmed structures
template <std::size_t DIM>
double jelmerStructure(const Geometry::PointPhysical<DIM>& pPhys,
                       PredefinedStructure structureType) {
    if (structureType == PredefinedStructure::VACUUM) {  // Vacuum Case
        return 1;
    }
    if (structureType == PredefinedStructure::BRAGG_STACK) {  // Bragg Stack
        if (pPhys[0] < 0.5) {
            return 13;
        }
        return 1;

    } else if (structureType == PredefinedStructure::CYLINDER) {  // Cylinder
        // Case with
        // radius 0.2a
        if ((pPhys[0] - 0.5) * (pPhys[0] - 0.5) +
                (pPhys[1] - 0.5) * (pPhys[1] - 0.5) <=
            0.2 * 0.2) {
            return 13;
            // std::cout << pPhys[0] << " " << pPhys[1] << " " << pPhys[2]  << "
            // " << epsilon_ << "\n";
        }
        return 1;
        // std::cout << pPhys[0] << " " << pPhys[1] << " " << pPhys[2]  << "
        // " << epsilon_ << "\n";

    } else if (structureType ==
               PredefinedStructure::SQUARE_HOLE) {  // Cube in Cuboid Case with
        // width of
        // pilars of 0.1a
        if (pPhys[0] < 0.1 || pPhys[0] > 0.9 || pPhys[1] < 0.1 ||
            pPhys[1] > 0.9) {
            return 1;
            // std::cout << pPhys[0] << " " << pPhys[1] << " " << pPhys[2]  << "
            // " << epsilon_ << "\n";
        }
        return 13;
        // std::cout << pPhys[0] << " " << pPhys[1] << " " << pPhys[2]  << "
        // " << epsilon_ << "\n";

    }

    else if (structureType ==
             PredefinedStructure::INVERSE_WOODPILE_OLD) {  // Inverse Woodpile
        // here y has length 10, whereas the lenght of x and y are
        // length(y)/sqrt(2) = 7.07.
        // The first three circles, 2 halves and 1 whole circle, are in x,y
        // plane. The second set of circles, 4 quarters and 1 whole, is in the
        // y,z plane. Diameter of cylinders are 0.19a.
        /*
        if((pPhys[0]-0.3535)*(pPhys[0]-0.3535) + (pPhys[1]-0.75)*(pPhys[1]-0.75)
        <= 0.19*0.19    || (pPhys[0]-0.707)*(pPhys[0]-0.707) +
        (pPhys[1]-0.25)*(pPhys[1]-0.25) <= 0.19*0.19      ||
           (pPhys[0])*(pPhys[0]) + (pPhys[1]-0.25)*(pPhys[1]-0.25) <= 0.19*0.19
        ||

           (pPhys[2]-0.3535)*(pPhys[2]-0.3535) + (pPhys[1]-0.5)*(pPhys[1]-0.5)
        <= 0.19*0.19      || (pPhys[2])*(pPhys[2]) + (pPhys[1])*(pPhys[1]) <=
        0.19*0.19                            || (pPhys[2])*(pPhys[2]) +
        (pPhys[1]-1)*(pPhys[1]-1) <= 0.19*0.19                        ||
           (pPhys[2]-0.707)*(pPhys[2]-0.707) + (pPhys[1])*(pPhys[1]) ||
           (pPhys[2]-0.707)*(pPhys[2]-0.707) + (pPhys[1]-1)*(pPhys[1]-1) <=
        0.19*0.19)
        {
            epsilon_ = 13;
            std::cout << pPhys[0] << " " << pPhys[1] << " " << pPhys[2]  << " "
        << epsilon_ << "\n";
        }
        else
        {
            epsilon_ = 1;
            std::cout << pPhys[0] << " " << pPhys[1] << " " << pPhys[2]  << " "
        << epsilon_ << "\n";
        }
         */

        if ((pPhys[2] - 0.3535) * (pPhys[2] - 0.3535) +
                    (pPhys[1] - 0.5) * (pPhys[1] - 0.5) <=
                0.19 * 0.19 ||
            (pPhys[2] - 0.707) * (pPhys[2] - 0.707) + (pPhys[1]) * (pPhys[1]) <=
                0.19 * 0.19 ||
            (pPhys[2] - 0.707) * (pPhys[2] - 0.707) +
                    (pPhys[1] - 1) * (pPhys[1] - 1) <=
                0.19 * 0.19 ||
            (pPhys[2]) * (pPhys[2]) + (pPhys[1]) * (pPhys[1]) <= 0.19 * 0.19 ||
            (pPhys[2]) * (pPhys[2]) + (pPhys[1] - 1) * (pPhys[1] - 1) <=
                0.19 * 0.19 ||

            (pPhys[0] - 0.3535) * (pPhys[0] - 0.3535) +
                    (pPhys[1] - 0.75) * (pPhys[1] - 0.75) <=
                0.19 * 0.19 ||
            (pPhys[0]) * (pPhys[0]) + (pPhys[1] - 0.25) * (pPhys[1] - 0.25) <=
                0.19 * 0.19 ||
            (pPhys[0] - 0.707) * (pPhys[0] - 0.707) +
                    (pPhys[1] - 0.25) * (pPhys[1] - 0.25) <=
                0.19 * 0.19) {
            return 1;
            // std::cout << pPhys[0] << " " << pPhys[1] << " " << pPhys[2]  << "
            // " << epsilon_ << "\n";
        }
        return 13;
        // std::cout << pPhys[0] << " " << pPhys[1] << " " << pPhys[2]  << "
        // " << epsilon_ << "\n";

    } else if (structureType == PredefinedStructure::INVERSE_WOODPILE_NEW) {
        const double a = 1.414;
        const double c = 1.0;
        const double r = 0.3111;
        const double a8 = a / 8;

        // Offsets, for a standard woodpile the pores in X and Z direction are
        // spaced a/4 apart. So it is expected that the multipliers differ by
        // either 2*a8 (=a/4) or 6*a8 (=3a/4)
        double offXPores = +1.0 * a8;
        double offZPores = -1.0 * a8;

        if (insideIWPore(pPhys[0], pPhys[1] + offZPores, a, c, r) ||
            insideIWPore(pPhys[2], pPhys[1] + offXPores, a, c, r)) {
            return 1.0;
        } else {
            // Use by default the same value as in COPS MPB computations.
            return 12.1;
        }

    } else {
        std::cout << "Incorrect value for SetEpsilon"
                  << "\n";
        return 1.0;
    }
}

}  // namespace DGMax