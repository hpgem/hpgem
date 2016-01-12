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

//This test checks that the basis functions are computed correctly.
//Then it tries to provoke the memory management system into assigning a different basis function to the same place in memory
//by repeatedly deleting and adding the basis functions
//Then it checks that the basis functions are still computed correctly
//(regression test for the recompute cache)

//interestingly, using memory debugging tools like libgmalloc.dylib or valgrind is more likely to fix the test than to point at the source of error

#include "Base/MeshManipulator.h"
#include "Base/Element.h"
#include "Base/ConfigurationData.h"
#include "Base/CommandLineOptions.h"

template<std::size_t DIM>
void testData(const Base::MeshManipulator<DIM>& mesh)
{
    for(Base::Element* element : mesh.getElementsList())
    {
        for(std::size_t i = 0; i < 4; ++i)
        {
            const Geometry::PointReference<DIM>& node = element->getReferenceGeometry()->getReferenceNodeCoordinate(i);
            logger.assert_always(std::abs(element->getSolution(0, node)[0] - double(i)) < 1e-12, "Solution changed");
        }
    }
}

int main(int argc, char** argv)
{
    Base::parse_options(argc, argv);
    //this test should also be effective in 1D , but 2D has 3x as much 'wrong' basis functions for only a little extra effort
    Base::ConfigurationData* config = new Base::ConfigurationData(2, 1, 1);
    Base::MeshManipulator<2> mesh(config);
    Geometry::PointPhysical<2> bottomLeft(LinearAlgebra::SmallVector<2>{{0., 0.}});
    Geometry::PointPhysical<2> topRight(LinearAlgebra::SmallVector<2>{{1., 1.}});
    mesh.createRectangularMesh(bottomLeft, topRight, {{1, 1}});
    mesh.useDefaultDGBasisFunctions();
    for(Base::Element* element : mesh.getElementsList())
    {
        element->setNumberOfTimeIntegrationVectors(1);
        element->setTimeIntegrationSubvector(0, 0, {{0., 1., 2., 3.}});
    }
    testData(mesh);
    //increase the number of iterations if this test is failing inconsistently
    for(std::size_t i = 0; i < 1000; ++i)
    {
        mesh.useDefaultDGBasisFunctions();
        for(Base::Element* element : mesh.getElementsList())
        {
            element->setTimeIntegrationSubvector(0, 0, {{0., 1., 2., 3.}});
        }
        testData(mesh);
    }
    //for a single element conforming basis functions happen to be the same
    for(std::size_t i = 0; i < 1000; ++i)
    {
        mesh.useDefaultConformingBasisFunctions();
        for(Base::Element* element : mesh.getElementsList())
        {
            element->setTimeIntegrationSubvector(0, 0, {{0., 1., 2., 3.}});
        }
        testData(mesh);
    }
    delete config;
    return 0;
}
