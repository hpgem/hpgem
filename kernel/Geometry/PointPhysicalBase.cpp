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

#include "PointPhysical.h"

namespace Geometry {

PointPhysicalBase::operator PointPhysical<0> &() {
    logger.assert_debug(typeid(*this) == typeid(PointPhysical<0>),
                        "Trying to convert to a point of the wrong dimension");
    return static_cast<PointPhysical<0>&>(*this);
}

PointPhysicalBase::operator PointPhysical<1> &() {
    logger.assert_debug(typeid(*this) == typeid(PointPhysical<1>),
                        "Trying to convert to a point of the wrong dimension");
    return static_cast<PointPhysical<1>&>(*this);
}

PointPhysicalBase::operator PointPhysical<2> &() {
    logger.assert_debug(typeid(*this) == typeid(PointPhysical<2>),
                        "Trying to convert to a point of the wrong dimension");
    return static_cast<PointPhysical<2>&>(*this);
}

PointPhysicalBase::operator PointPhysical<3> &() {
    logger.assert_debug(typeid(*this) == typeid(PointPhysical<3>),
                        "Trying to convert to a point of the wrong dimension");
    return static_cast<PointPhysical<3>&>(*this);
}

PointPhysicalBase::operator PointPhysical<4> &() {
    logger.assert_debug(typeid(*this) == typeid(PointPhysical<4>),
                        "Trying to convert to a point of the wrong dimension");
    return static_cast<PointPhysical<4>&>(*this);
}

PointPhysicalBase::operator const PointPhysical<0> &() const {
    logger.assert_debug(typeid(*this) == typeid(PointPhysical<0>),
                        "Trying to convert to a point of the wrong dimension");
    return static_cast<const PointPhysical<0>&>(*this);
}

PointPhysicalBase::operator const PointPhysical<1> &() const {
    logger.assert_debug(typeid(*this) == typeid(PointPhysical<1>),
                        "Trying to convert to a point of the wrong dimension");
    return static_cast<const PointPhysical<1>&>(*this);
}

PointPhysicalBase::operator const PointPhysical<2> &() const {
    logger.assert_debug(typeid(*this) == typeid(PointPhysical<2>),
                        "Trying to convert to a point of the wrong dimension");
    return static_cast<const PointPhysical<2>&>(*this);
}

PointPhysicalBase::operator const PointPhysical<3> &() const {
    logger.assert_debug(typeid(*this) == typeid(PointPhysical<3>),
                        "Trying to convert to a point of the wrong dimension");
    return static_cast<const PointPhysical<3>&>(*this);
}

PointPhysicalBase::operator const PointPhysical<4> &() const {
    logger.assert_debug(typeid(*this) == typeid(PointPhysical<4>),
                        "Trying to convert to a point of the wrong dimension");
    return static_cast<const PointPhysical<4>&>(*this);
}

// output operator can just pretend it is a point of appropriate dimension
std::ostream& operator<<(std::ostream& out, const PointPhysicalBase& point) {
    switch (point.size()) {
        case 0:
            out << Point<0>(static_cast<const PointPhysical<0>&>(point));
            break;
        case 1:
            out << Point<1>(static_cast<const PointPhysical<1>&>(point));
            break;
        case 2:
            out << Point<2>(static_cast<const PointPhysical<2>&>(point));
            break;
        case 3:
            out << Point<3>(static_cast<const PointPhysical<3>&>(point));
            break;
        case 4:
            out << Point<4>(static_cast<const PointPhysical<4>&>(point));
    }
    return out;
}
}  // namespace Geometry
