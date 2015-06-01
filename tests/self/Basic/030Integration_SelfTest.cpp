/*
 This file forms part of hpGEM. This package has been developed over a number of years by various people at the University of Twente and a full list of contributors can be found at
 http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found below.


 Copyright (c) 2014, University of Twente
 All rights reserved.

 Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

//Validates that the mesh contains no gaps or overlaps and that integration works properly by integrating a
//series of (non-)linear functions over the entire domain.
#include "Base/MeshManipulator.h"
#include "Base/RectangularMeshDescriptor.h"
#include "Integration/ElementIntegral.h"
#include "Base/ConfigurationData.h"
#include "Base/ElementCacheData.h"

#include "Base/CommandLineOptions.h"
#include "unordered_set"
#include "Logger.h"
#include <cmath>

template<std::size_t DIM>
void testMesh(Base::MeshManipulator<DIM>* test)
{

    class : public Integration::ElementIntegrandBase<double, DIM>
    {
        void elementIntegrand(const Base::Element* el, const Geometry::PointReference<DIM>& p, double& ret)
        {
            ret = 1.;
        }
    } one;
    
    class : public Integration::ElementIntegrandBase<double, DIM>
    {
        void elementIntegrand(const Base::Element* el, const Geometry::PointReference<DIM>& p, double& ret)
        {
            ret = 0;
            Geometry::PointPhysical<DIM> pPhys = el->referenceToPhysical(p);
            for (std::size_t i = 0; i < DIM; ++i)
            {
                ret += pPhys[i];
            }
        }
    } linear;
    
    class : public Integration::ElementIntegrandBase<double, DIM>
    {
        void elementIntegrand(const Base::Element* el, const Geometry::PointReference<DIM>& p, double& ret)
        {
            ret = 1;
            Geometry::PointPhysical<DIM> pPhys = el->referenceToPhysical(p);
            for (std::size_t i = 0; i < DIM; ++i)
            {
                ret *= pPhys[i];
            }
        }
    } trilinear;
    
    Integration::ElementIntegral elIntegral(false);
    double total = 0;
    double result;
    for (Base::Element* element : test->getElementsList())
    {
        result = elIntegral.integrate(element, &one);
        total += result;
    }
    logger.assert_always((std::abs(total - 1.) < 1e-12), "total mesh volume");
    total = 0;
    for (Base::Element* element : test->getElementsList())
    {
        result = elIntegral.integrate(element, &linear);
        total += result;
    }
    logger.assert_always((std::abs(total - .5 * test->dimension()) < 1e-12), "linear function");
    total = 0;
    for (Base::Element* element : test->getElementsList())
    {
        result = elIntegral.integrate(element, &trilinear);
        total += result;
    }
    logger.assert_always((std::abs(total - std::pow(0.5, test->dimension())) < 1e-12), "trilinear function");
}

int main(int argc, char** argv)
{
    Base::parse_options(argc, argv);
    // dim 1
    Base::RectangularMeshDescriptor<1> description1D;
    Base::RectangularMeshDescriptor<2> description2D;
    Base::RectangularMeshDescriptor<3> description3D;
    description1D.bottomLeft_[0] = 0;
    description2D.bottomLeft_[0] = 0;
    description2D.bottomLeft_[1] = 0;
    description3D.bottomLeft_[0] = 0;
    description3D.bottomLeft_[1] = 0;
    description3D.bottomLeft_[2] = 0;
    description1D.topRight_[0] = 1;
    description2D.topRight_[0] = 1;
    description2D.topRight_[1] = 1;
    description3D.topRight_[0] = 1;
    description3D.topRight_[1] = 1;
    description3D.topRight_[2] = 1;
    description1D.boundaryConditions_[0] = Base::BoundaryType::SOLID_WALL;
    description2D.boundaryConditions_[0] = Base::BoundaryType::SOLID_WALL;
    description2D.boundaryConditions_[1] = Base::BoundaryType::SOLID_WALL;
    description3D.boundaryConditions_[0] = Base::BoundaryType::SOLID_WALL;
    description3D.boundaryConditions_[1] = Base::BoundaryType::SOLID_WALL;
    description3D.boundaryConditions_[2] = Base::BoundaryType::SOLID_WALL;
    
    description1D.numElementsInDIM_[0] = 2;
    
    Base::MeshManipulator<1> *test = new Base::MeshManipulator<1>(new Base::ConfigurationData(1, 1, 2, 0), Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, 2, 0);
    test->createTriangularMesh(description1D.bottomLeft_, description1D.topRight_, description1D.numElementsInDIM_);
    
    testMesh(test);
    
    delete test;
    test = new Base::MeshManipulator<1>(new Base::ConfigurationData(1, 1, 2, 0), Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, 2, 0);
    test->createRectangularMesh(description1D.bottomLeft_, description1D.topRight_, description1D.numElementsInDIM_);
    testMesh(test);
    
    delete test;
    description1D.numElementsInDIM_[0] = 3;
    
    test = new Base::MeshManipulator<1>(new Base::ConfigurationData(1, 1, 2, 0), Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, 2, 0);
    test->createTriangularMesh(description1D.bottomLeft_, description1D.topRight_, description1D.numElementsInDIM_);
    
    testMesh(test);
    
    delete test;
    test = new Base::MeshManipulator<1>(new Base::ConfigurationData(1, 1, 2, 0), Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, 2, 0);
    test->createRectangularMesh(description1D.bottomLeft_, description1D.topRight_, description1D.numElementsInDIM_);
    testMesh(test);
    
    // dim 2
    
    delete test;
    description2D.numElementsInDIM_[0] = 2;
    description2D.numElementsInDIM_[1] = 3;
    
    Base::MeshManipulator<2> *test2 = new Base::MeshManipulator<2>(new Base::ConfigurationData(2, 1, 2, 0), Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, 2, 0);
    test2->createTriangularMesh(description2D.bottomLeft_, description2D.topRight_, description2D.numElementsInDIM_);
    
    testMesh(test2);
    
    delete test2;
    test2 = new Base::MeshManipulator<2>(new Base::ConfigurationData(2, 1, 2, 0), Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, 2, 0);
    test2->createRectangularMesh(description2D.bottomLeft_, description2D.topRight_, description2D.numElementsInDIM_);
    testMesh(test2);
    
    delete test2;
    description2D.numElementsInDIM_[0] = 3;
    description2D.numElementsInDIM_[1] = 2;
    
    test2 = new Base::MeshManipulator<2>(new Base::ConfigurationData(2, 1, 2, 0), Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, 2, 0);
    test2->createTriangularMesh(description2D.bottomLeft_, description2D.topRight_, description2D.numElementsInDIM_);
    
    testMesh(test2);
    
    delete test2;
    test2 = new Base::MeshManipulator<2>(new Base::ConfigurationData(2, 1, 2, 0), Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, 2, 0);
    test2->createRectangularMesh(description2D.bottomLeft_, description2D.topRight_, description2D.numElementsInDIM_);
    testMesh(test2);
    
    // dim 3
    
    delete test2;
    description3D.numElementsInDIM_[0] = 2;
    description3D.numElementsInDIM_[1] = 2;
    description3D.numElementsInDIM_[2] = 3;
    
    Base::MeshManipulator<3> *test3 = new Base::MeshManipulator<3>(new Base::ConfigurationData(3, 1, 2, 0), Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, 2, 0);
    test3->createTriangularMesh(description3D.bottomLeft_, description3D.topRight_, description3D.numElementsInDIM_);
    
    testMesh(test3);
    
    delete test3;
    test3 = new Base::MeshManipulator<3>(new Base::ConfigurationData(3, 1, 2, 0), Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, 2, 0);
    test3->createRectangularMesh(description3D.bottomLeft_, description3D.topRight_, description3D.numElementsInDIM_);
    testMesh(test3);
    
    delete test3;
    description3D.numElementsInDIM_[0] = 2;
    description3D.numElementsInDIM_[1] = 3;
    description3D.numElementsInDIM_[2] = 2;
    
    test3 = new Base::MeshManipulator<3>(new Base::ConfigurationData(3, 1, 2, 0), Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, 2, 0);
    test3->createTriangularMesh(description3D.bottomLeft_, description3D.topRight_, description3D.numElementsInDIM_);
    
    testMesh(test3);
    
    delete test3;
    test3 = new Base::MeshManipulator<3>(new Base::ConfigurationData(3, 1, 2, 0), Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, 2, 0);
    test3->createRectangularMesh(description3D.bottomLeft_, description3D.topRight_, description3D.numElementsInDIM_);
    testMesh(test3);
    
    delete test3;
    description3D.numElementsInDIM_[0] = 3;
    description3D.numElementsInDIM_[1] = 2;
    description3D.numElementsInDIM_[2] = 2;
    
    test3 = new Base::MeshManipulator<3>(new Base::ConfigurationData(3, 1, 2, 0), Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, 2, 0);
    test3->createTriangularMesh(description3D.bottomLeft_, description3D.topRight_, description3D.numElementsInDIM_);
    
    testMesh(test3);
    
    delete test3;
    test3 = new Base::MeshManipulator<3>(new Base::ConfigurationData(3, 1, 2, 0), Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, 2, 0);
    test3->createRectangularMesh(description3D.bottomLeft_, description3D.topRight_, description3D.numElementsInDIM_);
    testMesh(test3);
}

