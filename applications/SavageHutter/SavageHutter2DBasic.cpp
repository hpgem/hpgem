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

#include <fstream>
#include <iomanip>
#include <algorithm>

#include "hpgem-cmake.h"

#include "SavageHutter2DBasic.h"
#include "HeightLimiters/PositiveLayerLimiter.h"

using namespace hpgem;

SavageHutter2DBasic::SavageHutter2DBasic(std::size_t polyOrder,
                                         std::string meshName)
    : SavageHutter2DBase(3, polyOrder) {
    chuteAngle_ = M_PI / 180 * 29.6483997880436;
    epsilon_ = .1;
    // const PointPhysicalT &pPhys = createMeshDescription(1).bottomLeft_;
    const PointPhysicalT pPhys = {0., 0.};
    // inflowBC_ = getInitialSolution(pPhys, 0);
    inflowBC_ = LinearAlgebra::MiddleSizeVector(
        {1, 3 * std::sqrt(std::cos(chuteAngle_) * epsilon_ * 1), 0});
    logger(INFO, "inflow BC: %", inflowBC_);

    std::vector<std::string> variableNames = {"h", "hu", "hv"};
    setOutputNames("output2D", "SavageHutter", "SavageHutter", variableNames);

    readMesh(meshName);

    createContraction();
}

///\details Transform the rectangle to a contracting geometry.
void SavageHutter2DBasic::createContraction() {
    // Set up the move of the mesh; note that the mesh mover gets deleted in the
    // mesh manipulator
    MeshMoverContraction* contraction = new MeshMoverContraction;
    initialiseMeshMover(contraction, 0);

    // Set the parameters for the contraction. Not necessary to do here,
    // but it is nice to have everything at one place.
    contraction->xBegin = 6;
    contraction->xMiddle = 11;
    contraction->xEnd = 12;
    contraction->contractionWidth = .6;

    // Actually move the mesh such that a contraction is formed.
    meshes_[0]->move();

    // fix all elements and faces again, for example normals to faces.
    for (Base::Element* element :
         meshes_[0]->getElementsList(Base::IteratorType::GLOBAL)) {
        element->getReferenceToPhysicalMap()->reinit();
    }
}

LinearAlgebra::MiddleSizeVector SavageHutter2DBasic::getInitialSolution(
    const PointPhysicalT& pPhys, const double& startTime,
    const std::size_t orderTimeDerivative) {
    const double x = pPhys[0];
    /*std::map<double, std::array<double, 3>>::iterator positionInMap =
    inputValues_.lower_bound(x); std::array<double, 3> value =
    positionInMap->second; logger(DEBUG, "initial conditions at x = %: %, %", x,
    value[0], value[1]); return LinearAlgebra::MiddleSizeVector({value[0],
    value[1], value[2]});*/
    double h = 1;
    double froudeIn = 3;
    double hu = h * froudeIn * std::sqrt(std::cos(chuteAngle_) * epsilon_ * h);
    if (x > .5) {
        double u = hu / h;
        h = (-h + std::sqrt(h * h + 4 * u * u * h * 2 / epsilon_ /
                                        std::cos(chuteAngle_))) /
            2;
        logger(DEBUG, "h initial: %", h);
    }
    const double hv = 0;
    return MiddleSizeVector({h, hu, hv});
}

void SavageHutter2DBasic::registerVTKWriteFunctions() {
    HpgemAPISimplified::registerVTKWriteFunctions();

    registerVTKWriteFunction(
        [=](Base::Element* element, const Geometry::PointReference<2>& pRef,
            std::size_t timeLevel) -> double {
            if (element->getSolution(timeLevel, pRef)[0] > 1e-5) {
                const double h = element->getSolution(timeLevel, pRef)[0];
                const double u = element->getSolution(timeLevel, pRef)[1] / h;
                const double v = element->getSolution(timeLevel, pRef)[2] / h;
                return Base::L2Norm(LinearAlgebra::MiddleSizeVector({u, v})) /
                       std::sqrt(h * epsilon_ * std::cos(chuteAngle_));
            }
            return 0;
        },
        "F");
}

HeightLimiter* SavageHutter2DBasic::createHeightLimiter() {
    return new EmptyHeightLimiter;
}

LinearAlgebra::MiddleSizeVector SavageHutter2DBasic::computeSourceTerm(
    const LinearAlgebra::MiddleSizeVector& numericalSolution,
    const PointPhysicalT& pPhys, const double time) {
    logger.assert_debug(chuteAngle_ < M_PI / 2,
                        "Angle must be in radians, not degrees!");
    const double h = numericalSolution(0);
    const double hu = numericalSolution(1);
    const double hv = numericalSolution(2);

    double u = 0;
    double v = 0;
    if (h > dryLimit_) {
        u = hu / h;
        v = hv / h;
    }

    double uNormalized = 0;
    double vNormalized = 0;
    if (Base::L2Norm({u, v}) > 1e-16) {
        uNormalized = u / Base::L2Norm({u, v});
        vNormalized = v / Base::L2Norm({u, v});
    }
    const double mu = computeFriction(numericalSolution);
    const double sourceX = h * std::sin(chuteAngle_) -
                           h * mu * uNormalized * std::cos(chuteAngle_);
    const double sourceY = -h * mu * vNormalized * std::cos(chuteAngle_);
    return MiddleSizeVector({0, sourceX, sourceY});
}

LinearAlgebra::MiddleSizeVector SavageHutter2DBasic::computePhysicalFlux(
    const LinearAlgebra::MiddleSizeVector& numericalSolution) {
    const double h = numericalSolution(0);
    logger.assert_debug(h > -1e-16, "Negative height (%)", h);
    double hu = numericalSolution(1);
    double hv = numericalSolution(2);
    double u = 0;
    double v = 0;
    if (h > dryLimit_) {
        u = hu / h;
        v = hv / h;
    }
    MiddleSizeVector flux(6);
    flux(0) = hu;
    flux(1) = hv;
    flux(2) = hu * u + epsilon_ / 2 * std::cos(chuteAngle_) * h * h;
    flux(3) = hu * v;
    flux(4) = hu * v;
    flux(5) = hv * v + epsilon_ / 2 * std::cos(chuteAngle_) * h * h;
    logger(DEBUG, "flux values: %, ", flux);
    return flux;
}

void SavageHutter2DBasic::tasksAfterSolving() {
    if (false) {
        auto widthValues = widthAverage();

        std::ofstream widthFile("widthFile.dat");
        for (auto myPair : widthValues) {
            widthFile << myPair.first;
            for (std::size_t i = 0; i < 3; ++i) {
                widthFile << '\t' << std::setw(10) << myPair.second[i];
            }
            widthFile << '\n';
        }
    }
}

void SavageHutter2DBasic::tasksBeforeSolving() {
    SavageHutterBase::tasksBeforeSolving();
    if (false) {
        // put values of initial conditions from text file in a map
        std::ifstream inFile(
            getCMAKE_hpGEM_SOURCE_DIR() +
            "/applications/SavageHutter/"
            "initialConditionsGeneratedBy1DMatlabWithShock.txt");
        std::string line;
        while (getline(inFile, line)) {
            std::stringstream lineStream(line);
            double x;
            std::array<double, 2> hAndU;
            std::array<double, 3> valuesForVariables;
            lineStream >> x >> hAndU[0] >> hAndU[1];
            valuesForVariables = {hAndU[0], hAndU[0] * hAndU[1], 0};
            inputValues_[x] = valuesForVariables;
        }
        logger(INFO, "size of the map after reading: %", inputValues_.size());
    }
}
