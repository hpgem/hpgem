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

// Validates that a moved mesh contains no gaps or overlaps by integrating a
// series of (non-)linear functions over the entire domain.
#include "Base/MeshManipulator.h"
#include "Integration/ElementIntegral.h"
#include "Base/ConfigurationData.h"
#include "Base/MeshMoverBase.h"
#include "Geometry/Mappings/MappingReferenceToPhysical.h"
#include "Base/CommandLineOptions.h"

#include "unordered_set"
#include "Logger.h"
#include <cmath>
#include <CMakeDefinitions.h>

template <std::size_t DIM>
void move(Base::MeshManipulator<DIM>* mesh) {
    class : public Base::MeshMoverBase<DIM> {
       public:
        void movePoint(Geometry::PointPhysical<DIM>& p) const { p *= 2; }
    } mover;
    for (Geometry::PointPhysical<DIM>& node : mesh->getNodeCoordinates()) {
        mover.movePoint(node);
    }
    for (Base::Element* element : mesh->getElementsList()) {
        element->getReferenceToPhysicalMap()->reinit();
    }
}

template <std::size_t DIM>
void testMesh(Base::MeshManipulator<DIM>* test) {
    move(test);
    class : public Integration::ElementIntegrandBase<double, DIM> {
        void elementIntegrand(Base::PhysicalElement<DIM>& element,
                              double& ret) {
            ret = 1;
        }
    } one;

    class : public Integration::ElementIntegrandBase<double, DIM> {
        void elementIntegrand(Base::PhysicalElement<DIM>& element,
                              double& ret) {
            ret = 0;
            const Geometry::PointPhysical<DIM>& pPhys =
                element.getPointPhysical();
            for (std::size_t i = 0; i < pPhys.size(); ++i) {
                ret += pPhys[i];
            }
        }
    } linear;
    class : public Integration::ElementIntegrandBase<double, DIM> {
        void elementIntegrand(Base::PhysicalElement<DIM>& element,
                              double& ret) {
            ret = 1;
            const Geometry::PointPhysical<DIM>& pPhys =
                element.getPointPhysical();
            for (std::size_t i = 0; i < pPhys.size(); ++i) {
                ret *= pPhys[i];
            }
        }
    } trilinear;
    Integration::ElementIntegral<DIM> elIntegral;
    double total = 0;
    double result;
    for (Base::Element* element : test->getElementsList()) {
        result = elIntegral.integrate(element, &one);
        total += result;
    }
    logger.assert_always(
        (std::abs(total - std::pow(2., test->dimension())) < 1e-12),
        "total mesh volume");
    total = 0;
    for (Base::Element* element : test->getElementsList()) {
        result = elIntegral.integrate(element, &linear);
        total += result;
    }
    logger.assert_always(
        (std::abs(total - test->dimension() * std::pow(2., test->dimension())) <
         1e-12),
        "linear function");
    total = 0;
    for (Base::Element* element : test->getElementsList()) {
        result = elIntegral.integrate(element, &trilinear);
        total += result;
    }
    logger.assert_always(
        (std::abs(total - std::pow(2., test->dimension())) < 1e-12),
        "trilinear function");
}

int main(int argc, char** argv) {
    using namespace std::string_literals;
    Base::parse_options(argc, argv);
    // dim 1

    Base::MeshManipulator<1>* test =
        new Base::MeshManipulator<1>(new Base::ConfigurationData(1, 0), 0);
    test->readMesh(Base::getCMAKE_hpGEM_SOURCE_DIR() + "/tests/files/"s +
                   "1Drectangular1mesh.hpgem"s);

    test->useMonomialBasisFunctions(2);
    testMesh(test);

    delete test;
    test = new Base::MeshManipulator<1>(new Base::ConfigurationData(1));
    test->readMesh(Base::getCMAKE_hpGEM_SOURCE_DIR() + "/tests/files/"s +
                   "1Drectangular2mesh.hpgem"s);
    test->useMonomialBasisFunctions(2);
    testMesh(test);

    // dim 2

    delete test;

    Base::MeshManipulator<2>* test2 =
        new Base::MeshManipulator<2>(new Base::ConfigurationData(1));
    test2->readMesh(Base::getCMAKE_hpGEM_SOURCE_DIR() + "/tests/files/"s +
                    "2Dtriangular1mesh.hpgem"s);
    test2->useMonomialBasisFunctions(2);
    testMesh(test2);

    delete test2;
    test2 = new Base::MeshManipulator<2>(new Base::ConfigurationData(1));
    test2->readMesh(Base::getCMAKE_hpGEM_SOURCE_DIR() + "/tests/files/"s +
                    "2Drectangular1mesh.hpgem"s);
    test2->useMonomialBasisFunctions(2);
    testMesh(test2);

    delete test2;

    test2 = new Base::MeshManipulator<2>(new Base::ConfigurationData(1));
    test2->readMesh(Base::getCMAKE_hpGEM_SOURCE_DIR() + "/tests/files/"s +
                    "2Dtriangular2mesh.hpgem"s);
    test2->useMonomialBasisFunctions(2);
    testMesh(test2);

    delete test2;
    test2 = new Base::MeshManipulator<2>(new Base::ConfigurationData(1));
    test2->readMesh(Base::getCMAKE_hpGEM_SOURCE_DIR() + "/tests/files/"s +
                    "2Drectangular2mesh.hpgem"s);
    test2->useMonomialBasisFunctions(2);
    testMesh(test2);

    // dim 3

    delete test2;

    Base::MeshManipulator<3>* test3 =
        new Base::MeshManipulator<3>(new Base::ConfigurationData(1));
    test3->readMesh(Base::getCMAKE_hpGEM_SOURCE_DIR() + "/tests/files/"s +
                    "3Dtriangular1mesh.hpgem"s);
    test3->useMonomialBasisFunctions(2);
    testMesh(test3);

    delete test3;
    test3 = new Base::MeshManipulator<3>(new Base::ConfigurationData(1));
    test3->readMesh(Base::getCMAKE_hpGEM_SOURCE_DIR() + "/tests/files/"s +
                    "3Drectangular1mesh.hpgem"s);
    test3->useMonomialBasisFunctions(2);
    testMesh(test3);

    delete test3;

    test3 = new Base::MeshManipulator<3>(new Base::ConfigurationData(1));
    test3->readMesh(Base::getCMAKE_hpGEM_SOURCE_DIR() + "/tests/files/"s +
                    "3Dtriangular2mesh.hpgem"s);
    test3->useMonomialBasisFunctions(2);
    testMesh(test3);

    delete test3;
    test3 = new Base::MeshManipulator<3>(new Base::ConfigurationData(1));
    test3->readMesh(Base::getCMAKE_hpGEM_SOURCE_DIR() + "/tests/files/"s +
                    "3Drectangular2mesh.hpgem"s);
    test3->useMonomialBasisFunctions(2);
    testMesh(test3);

    delete test3;

    test3 = new Base::MeshManipulator<3>(new Base::ConfigurationData(1));
    test3->readMesh(Base::getCMAKE_hpGEM_SOURCE_DIR() + "/tests/files/"s +
                    "3Dtriangular3mesh.hpgem"s);
    test3->useMonomialBasisFunctions(2);
    testMesh(test3);

    delete test3;
    test3 = new Base::MeshManipulator<3>(new Base::ConfigurationData(1));
    test3->readMesh(Base::getCMAKE_hpGEM_SOURCE_DIR() + "/tests/files/"s +
                    "3Drectangular3mesh.hpgem"s);
    test3->useMonomialBasisFunctions(2);
    testMesh(test3);
    delete test3;
}
