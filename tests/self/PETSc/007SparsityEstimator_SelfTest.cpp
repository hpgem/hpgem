/*
 This file forms part of hpGEM. This package has been developed over a number of years by various people at the University of Twente and a full list of contributors can be found at
 http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found below.


 Copyright (c) 2020, University of Twente
 All rights reserved.

 Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "CMakeDefinitions.h"

#include "Base/MeshManipulator.h"
#include "Base/CommandLineOptions.h"

#include "Utilities/GlobalIndexing.h"
#include "Utilities/SparsityEstimator.h"

using Layout = Utilities::GlobalIndexing::Layout;

struct TestConfiguration
{
    TestConfiguration(Layout layout, bool faceCoupling)
        : layout_ (layout), faceCoupling_ (faceCoupling)
    {}

    Layout layout_;
    bool faceCoupling_;
};

TestConfiguration TEST_CONFIGURATIONS[] = {
        TestConfiguration(Layout::SEQUENTIAL, false),
        TestConfiguration(Layout::SEQUENTIAL, true),
        TestConfiguration(Layout::BLOCKED_PROCESSOR, false),
        TestConfiguration(Layout::BLOCKED_PROCESSOR, true),
        TestConfiguration(Layout::BLOCKED_GLOBAL, false),
        TestConfiguration(Layout::BLOCKED_GLOBAL, true),
};

void testWithDGBasis(std::size_t unknowns, std::string meshFile)
{
    Base::ConfigurationData config (unknowns);
    Base::MeshManipulator<3> mesh (&config);

    using namespace std::string_literals;
    std::stringstream filename;
    filename << Base::getCMAKE_hpGEM_SOURCE_DIR() + "/tests/files/"s << meshFile << ".hpgem";
    mesh.readMesh(filename.str());

    mesh.useDefaultDGBasisFunctions(0);
    for (std::size_t unknown = 1; unknown < unknowns; ++unknown)
    {
        // Use more than zeroth order DG
        mesh.useDefaultDGBasisFunctions(unknown, unknown);
    }

    for (TestConfiguration configuration : TEST_CONFIGURATIONS)
    {
        Utilities::GlobalIndexing indexing(&mesh, configuration.layout_);
        Utilities::SparsityEstimator estimator(mesh, indexing);

        std::vector<int> owned, nonOwned;
        estimator.computeSparsityEstimate(owned, nonOwned, configuration.faceCoupling_);
        logger.assert_always(owned.size() == indexing.getNumberOfLocalBasisFunctions(), "Wrong size owned");
        logger.assert_always(nonOwned.size() == indexing.getNumberOfLocalBasisFunctions(), "Wrong size owned");

        for (const Base::Element *element : mesh.getElementsList())
        {
            std::size_t numberOfOwnDoFs = element->getTotalLocalNumberOfBasisFunctions();
            std::size_t numberOfNonOwnDoFs = 0;
            // Skip face loop if there is no face coupling
            if (configuration.faceCoupling_)
            {
                for (const Base::Face *face : element->getFacesList())
                {
                    if (face->getFaceType() == Geometry::FaceType::SUBDOMAIN_BOUNDARY
                        || face->getFaceType() == Geometry::FaceType::PERIODIC_SUBDOMAIN_BC)
                    {
                        numberOfNonOwnDoFs += face->getPtrOtherElement(element)->getTotalLocalNumberOfBasisFunctions();
                    }
                    else if (face->isInternal())
                    {
                        numberOfOwnDoFs += face->getPtrOtherElement(element)->getTotalLocalNumberOfBasisFunctions();
                    }
                }
            }
            for (std::size_t unknown = 0; unknown < unknowns; ++unknown)
            {
                // Check for each of the basis functions
                std::size_t dof = indexing.getProcessorLocalIndex(element, unknown);
                for (std::size_t dofOff = 0; dofOff < element->getLocalNumberOfBasisFunctions(unknown); ++dofOff)
                {
                    logger.assert_always(owned[dof + dofOff] == numberOfOwnDoFs,
                            "Expected % owned non zero entries but got %",
                            numberOfOwnDoFs, owned[dof + dofOff]);
                    logger.assert_always(nonOwned[dof + dofOff] == numberOfNonOwnDoFs,
                            "Expected % non owned non zero entries but got %",
                            numberOfNonOwnDoFs, nonOwned[dof + dofOff]);
                }
            }
        }
    }
}

template<typename GEOM, typename T>
T& selectByOwner(const GEOM* geom, T& owned, T& nonOwned)
{
    if(geom->isOwnedByCurrentProcessor())
        return owned;
    else
        return nonOwned;
}


void testConformingWith1DMesh()
{
    Base::ConfigurationData config (1);
    Base::MeshManipulator<1> mesh (&config);

    using namespace std::string_literals;
    std::stringstream filename;
    filename << Base::getCMAKE_hpGEM_SOURCE_DIR() + "/tests/files/"s << "1Drectangular2mesh"s << ".hpgem";
    mesh.readMesh(filename.str());

    mesh.useDefaultConformingBasisFunctions(3);

    for (TestConfiguration configuration : TEST_CONFIGURATIONS)
    {
        Utilities::GlobalIndexing indexing(&mesh, configuration.layout_);
        Utilities::SparsityEstimator estimator(mesh, indexing);

        std::vector<int> owned, nonOwned;
        estimator.computeSparsityEstimate(owned, nonOwned, configuration.faceCoupling_);
        logger.assert_always(owned.size() == indexing.getNumberOfLocalBasisFunctions(), "Wrong size owned");
        logger.assert_always(nonOwned.size() == indexing.getNumberOfLocalBasisFunctions(), "Wrong size owned");

        // Check for Element based DoF
        for (const Base::Element* element : mesh.getElementsList())
        {
            // The following DoFs are contributing
            // The DoFs from the element and its neighbours
            // The DoFs from the faces to element
            // The DoFs from the faces of the neighbouring elements
            std::size_t numberOfOwnDoFs = element->getTotalLocalNumberOfBasisFunctions();
            std::size_t numberOfNonOwnDoFs = 0;
            for (const Base::Face * face : element->getFacesList())
            {
                selectByOwner(face, numberOfOwnDoFs, numberOfNonOwnDoFs) += face->getTotalLocalNumberOfBasisFunctions();
                // If there is no face at the other side or we do not couple through the face matrix
                // skip the computation of the DoFs on the other side
                if (!face->isInternal() || !configuration.faceCoupling_)
                    continue;
                // The neighbouring element
                const Base::Element* otherElement = face->getPtrOtherElement(element);
                selectByOwner(otherElement, numberOfOwnDoFs, numberOfNonOwnDoFs)
                    += element->getTotalLocalNumberOfBasisFunctions();
                // Other face of the neighbour
                const Base::Face* otherFace = otherElement->getFace(otherElement->getFace(0) == face ? 1 : 0);
                selectByOwner(otherFace, numberOfOwnDoFs, numberOfNonOwnDoFs)
                    += otherFace->getTotalLocalNumberOfBasisFunctions();
            }
            // Check
            for (std::size_t unknown = 0; unknown < indexing.getNumberOfUnknowns(); ++unknown)
            {
                // Check for each of the basis functions
                std::size_t dof = indexing.getProcessorLocalIndex(element, unknown);
                for (std::size_t dofOff = 0; dofOff < element->getLocalNumberOfBasisFunctions(unknown); ++dofOff)
                {
                    logger.assert_always(owned[dof + dofOff] == numberOfOwnDoFs,
                            "Expected % owned non zero entries but got %",
                            numberOfOwnDoFs, owned[dof + dofOff]);
                    logger.assert_always(nonOwned[dof + dofOff] == numberOfNonOwnDoFs,
                            "Expected % non owned non zero entries but got %",
                            numberOfNonOwnDoFs, nonOwned[dof + dofOff]);
                }
            }
        }
        // Check for Face based DoFs
        for (const Base::Face* face : mesh.getFacesList())
        {
            // The contributions are from
            // The DoFs on the face
            // The DoFs from the element adjacent to the face it its faces
            // The DoFs from the element adjecent to those and their faces.
            std::vector<const Base::Element*> elements ({face->getPtrElementLeft()});
            if (face->isInternal())
                elements.emplace_back(face->getPtrElementRight());
            std::size_t numberOfOwnDoFs = face->getTotalLocalNumberOfBasisFunctions();
            std::size_t numberOfNonOwnDoFs = 0;
            for (const Base::Element* element : elements)
            {
                selectByOwner(element, numberOfOwnDoFs, numberOfNonOwnDoFs)
                    += element->getTotalLocalNumberOfBasisFunctions();
                const Base::Face* otherFace = element->getFace(element->getFace(0) == face ? 1 : 0);
                selectByOwner(otherFace, numberOfOwnDoFs, numberOfNonOwnDoFs)
                    += otherFace->getTotalLocalNumberOfBasisFunctions();
                // Next element
                if (!otherFace->isInternal() || !configuration.faceCoupling_)
                    continue;
                const Base::Element* nextElement = otherFace->getPtrOtherElement(element);
                selectByOwner(nextElement, numberOfOwnDoFs, numberOfNonOwnDoFs)
                    += nextElement->getTotalLocalNumberOfBasisFunctions();
                const Base::Face* nextFace = nextElement->getFace(nextElement->getFace(0) == otherFace ? 1 : 0);
                selectByOwner(nextFace, numberOfOwnDoFs, numberOfNonOwnDoFs)
                    += nextFace->getTotalLocalNumberOfBasisFunctions();
            }
            // Check
            for (std::size_t unknown = 0; unknown < indexing.getNumberOfUnknowns(); ++unknown)
            {
                // Check for each of the basis functions
                std::size_t dof = indexing.getProcessorLocalIndex(face, unknown);
                for (std::size_t dofOff = 0; dofOff < face->getLocalNumberOfBasisFunctions(unknown); ++dofOff)
                {
                    logger.assert_always(owned[dof + dofOff] == numberOfOwnDoFs,
                            "Expected % owned non zero entries but got %",
                            numberOfOwnDoFs, owned[dof + dofOff]);
                    logger.assert_always(nonOwned[dof + dofOff] == numberOfNonOwnDoFs,
                            "Expected % non owned non zero entries but got %",
                            numberOfNonOwnDoFs, nonOwned[dof + dofOff]);
                }
            }
        }
    }
}

int main(int argc, char** argv)
{
    using namespace std::string_literals;
    Base::parse_options(argc, argv);

    testWithDGBasis(1, "3Drectangular1mesh"s);
    testWithDGBasis(3, "3Drectangular1mesh"s);
    testWithDGBasis(1, "3Dtriangular1mesh"s);

    testConformingWith1DMesh();
}
