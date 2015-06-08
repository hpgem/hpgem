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

#include "Base/MeshManipulator.h"
#include "Base/RectangularMeshDescriptor.h"
#include "Integration/ElementIntegral.h"
#include "Integration/ElementIntegrandBase.h"

#include "unordered_set"
#include "Logger.h"

#include "Utilities/BasisFunctions1DH1ConformingLine.h"
#include "Utilities/BasisFunctions2DH1ConformingSquare.h"
#include "Utilities/BasisFunctions2DH1ConformingTriangle.h"
#include "Utilities/BasisFunctions3DH1ConformingCube.h"
#include "Utilities/BasisFunctions3DH1ConformingTetrahedron.h"
#include "Base/ConfigurationData.h"
#include "Base/L2Norm.h"
#include "Base/ElementCacheData.h"
#include "Base/CommandLineOptions.h"
#include <cmath>

//computes \int_\Omega (\nabla f)^2 dx by interpolating f and then integrating using basisFunctionDerivatives
template<std::size_t DIM>
void testMesh(Base::MeshManipulator<DIM>* test)
{
    
    class : public Integration::ElementIntegrandBase<LinearAlgebra::MiddleSizeVector, DIM>
    {
        
        void elementIntegrand(Base::PhysicalElement<DIM>& element, LinearAlgebra::MiddleSizeVector& ret)
        {
            std::size_t numBasisFuns = element.getElement()->getNrOfBasisFunctions();
            ret.resize(numBasisFuns);
            const Geometry::PointPhysical<DIM>& pPhys = element.getPointPhysical();
            for (std::size_t i = 0; i < numBasisFuns; ++i)
            {
                ret[i] = element.basisFunction(i);
                for (std::size_t j = 0; j < pPhys.size(); ++j)
                {
                    ret[i] *= pPhys[j];
                }
            }
        }
    } interpolation;
    
    class : public Integration::ElementIntegrandBase<LinearAlgebra::MiddleSizeMatrix, DIM>
    {
        
        void elementIntegrand(Base::PhysicalElement<DIM>& element, LinearAlgebra::MiddleSizeMatrix& ret)
        {
            std::size_t numBasisFuns = element.getElement()->getNrOfBasisFunctions();
            ret.resize(numBasisFuns, numBasisFuns);
            for (std::size_t i = 0; i < numBasisFuns; ++i)
            {
                for (std::size_t j = 0; j < numBasisFuns; ++j)
                {
                    ret(i, j) = element.basisFunction(i) * element.basisFunction(j);
                }
            }
        }
    } massMatrix;
    
    class : public Integration::ElementIntegrandBase<double, DIM>
    {
        void elementIntegrand(Base::PhysicalElement<DIM>& element, double& ret)
        {
            LinearAlgebra::SmallVector<DIM> temp1 = element.getSolutionDeriv()[0];
            ret = temp1 * temp1;
        }
    } integrating;
    
    std::cout.precision(14);
    Integration::ElementIntegral<DIM> elIntegral(false);
    double total = 0;
    double result;
    LinearAlgebra::MiddleSizeVector  expansion;
    LinearAlgebra::MiddleSizeMatrix M;
    for (Base::Element* element : test->getElementsList())
    {
        expansion = elIntegral.integrate(element, &interpolation);
        M = elIntegral.integrate(element, &massMatrix);
        
        M.solve(expansion);
        element->setTimeLevelData(0, expansion);
        
        result = elIntegral.integrate(element, &integrating);
        
        //std::cout<<result[0]<<std::endl;
        total += result;
    }
    
    std::cout << total << " " << std::endl;
    logger.assert_always((std::abs(total - 4. / 3. + 1. / 3. * test->dimension()) < 1e-12), "derivatives");
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
    //1D triangular meshes dont exist
    Base::MeshManipulator<1> *test = new Base::MeshManipulator<1>(new Base::ConfigurationData(1, 1, 2, 1), Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, 2, 0);
    test->createRectangularMesh(description1D.bottomLeft_, description1D.topRight_, description1D.numElementsInDIM_);
    testMesh(test);
    //we dont have user defined basis function sets to test this, so fall back to a hpGEM provided one
    test->setDefaultBasisFunctionSet(Utilities::createDGBasisFunctionSet1DH1Line(2));
    testMesh(test);
    test->useDefaultDGBasisFunctions();
    testMesh(test);
    
    //somewhere a lookup-table does not get cleaned when a meshManipulator is deleted
    delete test;
    description1D.numElementsInDIM_[0] = 3;
    test = new Base::MeshManipulator<1>(new Base::ConfigurationData(1, 1, 2, 1), Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, 2, 0);
    test->createRectangularMesh(description1D.bottomLeft_, description1D.topRight_, description1D.numElementsInDIM_);
    testMesh(test);
    //we dont have user defined basis function sets to test this, so fall back to a hpGEM provided one
    test->setDefaultBasisFunctionSet(Utilities::createDGBasisFunctionSet1DH1Line(2));
    testMesh(test);
    test->useDefaultDGBasisFunctions();
    testMesh(test);
    
    // dim 2
    
    delete test;
    description2D.numElementsInDIM_[0] = 2;
    description2D.numElementsInDIM_[1] = 3;
    
    Base::MeshManipulator<2> *test2 = new Base::MeshManipulator<2>(new Base::ConfigurationData(2, 1, 2, 1), Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, 2, 0);
    test2->createTriangularMesh(description2D.bottomLeft_, description2D.topRight_, description2D.numElementsInDIM_);
    
    testMesh(test2);
    //we dont have user defined basis function sets to test this, so fall back to a hpGEM provided one
    test2->setDefaultBasisFunctionSet(Utilities::createDGBasisFunctionSet2DH1Triangle(2));
    testMesh(test2);
    test2->useDefaultDGBasisFunctions();
    testMesh(test2);
    
    delete test2;
    test2 = new Base::MeshManipulator<2>(new Base::ConfigurationData(2, 1, 2, 1), Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, 2, 0);
    test2->createRectangularMesh(description2D.bottomLeft_, description2D.topRight_, description2D.numElementsInDIM_);
    testMesh(test2);
    //we dont have user defined basis function sets to test this, so fall back to a hpGEM provided one
    test2->setDefaultBasisFunctionSet(Utilities::createDGBasisFunctionSet2DH1Square(2));
    testMesh(test2);
    test2->useDefaultDGBasisFunctions();
    testMesh(test2);
    
    delete test2;
    description2D.numElementsInDIM_[0] = 3;
    description2D.numElementsInDIM_[1] = 2;
    
    test2 = new Base::MeshManipulator<2>(new Base::ConfigurationData(2, 1, 2, 1), Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, 2, 0);
    test2->createTriangularMesh(description2D.bottomLeft_, description2D.topRight_, description2D.numElementsInDIM_);
    
    testMesh(test2);
    //we dont have user defined basis function sets to test this, so fall back to a hpGEM provided one
    test2->setDefaultBasisFunctionSet(Utilities::createDGBasisFunctionSet2DH1Triangle(2));
    testMesh(test2);
    test2->useDefaultDGBasisFunctions();
    testMesh(test2);
    
    delete test2;
    test2 = new Base::MeshManipulator<2>(new Base::ConfigurationData(2, 1, 2, 1), Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, 2, 0);
    test2->createRectangularMesh(description2D.bottomLeft_, description2D.topRight_, description2D.numElementsInDIM_);
    testMesh(test2);
    //we dont have user defined basis function sets to test this, so fall back to a hpGEM provided one
    test2->setDefaultBasisFunctionSet(Utilities::createDGBasisFunctionSet2DH1Square(2));
    testMesh(test2);
    test2->useDefaultDGBasisFunctions();
    testMesh(test2);
    // dim 3
    
    delete test2;
    description3D.numElementsInDIM_[0] = 2;
    description3D.numElementsInDIM_[1] = 2;
    description3D.numElementsInDIM_[2] = 3;
    
    Base::ConfigurationData* configData = new Base::ConfigurationData(3, 1, 3, 1);
    Base::MeshManipulator<3> *test3 = new Base::MeshManipulator<3>(configData, Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, 3, 0);
    test3->createTriangularMesh(description3D.bottomLeft_, description3D.topRight_, description3D.numElementsInDIM_);
    
    testMesh(test3);
    //we dont have user defined basis function sets to test this, so fall back to a hpGEM provided one
    test3->setDefaultBasisFunctionSet(Utilities::createDGBasisFunctionSet3DH1Tetrahedron(3));
    testMesh(test3);
    test3->useDefaultDGBasisFunctions();
    testMesh(test3);
    
    delete test3;
    delete configData;
    
    configData = new Base::ConfigurationData(3, 1, 2, 1);
    test3 = new Base::MeshManipulator<3>(configData, Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, 2, 0);
    test3->createRectangularMesh(description3D.bottomLeft_, description3D.topRight_, description3D.numElementsInDIM_);
    testMesh(test3);
    //we dont have user defined basis function sets to test this, so fall back to a hpGEM provided one
    test3->setDefaultBasisFunctionSet(Utilities::createDGBasisFunctionSet3DH1Cube(2));
    testMesh(test3);
    test3->useDefaultDGBasisFunctions();
    testMesh(test3);
    
    delete test3;
    delete configData;
    
    description3D.numElementsInDIM_[0] = 2;
    description3D.numElementsInDIM_[1] = 3;
    description3D.numElementsInDIM_[2] = 2;
    
    configData = new Base::ConfigurationData(3, 1, 3, 1);
    test3 = new Base::MeshManipulator<3>(configData, Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, 3, 0);
    test3->createTriangularMesh(description3D.bottomLeft_, description3D.topRight_, description3D.numElementsInDIM_);
    
    testMesh(test3);
    //we dont have user defined basis function sets to test this, so fall back to a hpGEM provided one
    test3->setDefaultBasisFunctionSet(Utilities::createDGBasisFunctionSet3DH1Tetrahedron(3));
    testMesh(test3);
    test3->useDefaultDGBasisFunctions();
    testMesh(test3);
    
    delete test3;
    delete configData;
    
    configData = new Base::ConfigurationData(3, 1, 2, 1);
    test3 = new Base::MeshManipulator<3>(configData, Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, 2, 0);
    test3->createRectangularMesh(description3D.bottomLeft_, description3D.topRight_, description3D.numElementsInDIM_);
    testMesh(test3);
    //we dont have user defined basis function sets to test this, so fall back to a hpGEM provided one
    test3->setDefaultBasisFunctionSet(Utilities::createDGBasisFunctionSet3DH1Cube(2));
    testMesh(test3);
    test3->useDefaultDGBasisFunctions();
    testMesh(test3);
    
    delete test3;
    delete configData;
    
    description3D.numElementsInDIM_[0] = 3;
    description3D.numElementsInDIM_[1] = 2;
    description3D.numElementsInDIM_[2] = 2;
    
    configData = new Base::ConfigurationData(3, 1, 3, 1);
    test3 = new Base::MeshManipulator<3>(configData, Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, 3, 0);
    test3->createTriangularMesh(description3D.bottomLeft_, description3D.topRight_, description3D.numElementsInDIM_);
    
    testMesh(test3);
    //we dont have user defined basis function sets to test this, so fall back to a hpGEM provided one
    test3->setDefaultBasisFunctionSet(Utilities::createDGBasisFunctionSet3DH1Tetrahedron(3));
    testMesh(test3);
    test3->useDefaultDGBasisFunctions();
    testMesh(test3);
    
    delete test3;
    delete configData;
    
    configData = new Base::ConfigurationData(3, 1, 2, 1);
    test3 = new Base::MeshManipulator<3>(configData, Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, 2, 0);
    test3->createRectangularMesh(description3D.bottomLeft_, description3D.topRight_, description3D.numElementsInDIM_);
    testMesh(test3);
    //we dont have user defined basis function sets to test this, so fall back to a hpGEM provided one
    test3->setDefaultBasisFunctionSet(Utilities::createDGBasisFunctionSet3DH1Cube(2));
    testMesh(test3);
    test3->useDefaultDGBasisFunctions();
    testMesh(test3);
    
    delete test3;
    delete configData;
    
    return 0;
}

