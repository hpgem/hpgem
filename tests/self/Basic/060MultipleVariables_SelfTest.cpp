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

/* Test the functions that are added to deal with multiple variables for a system of PDE's.
 
Functions to test:
-ElementData::convertToSingleIndex
-Element::getSolution
-Face::convertToSingleIndex
 */

#include <iostream>

#include "Base/AssembleBasisFunctionSet.hpp"
#include "Base/BasisFunctionSet.hpp"
#include "Base/Element.hpp"
#include "Base/Face.hpp"
#include "Base/FaceCacheData.hpp"
#include "Base/Node.hpp"
#include "Geometry/PointPhysical.hpp"
#include "Geometry/PointReference.hpp"

int main()
{
    std::cout << "Test functions for multiple variables.\n";
    
    // Build two square elements and one face in between them.
    const std::size_t dimension = 2;
    const std::size_t nrOfUnknowns = 5;
    const std::size_t nrOfTimeLevels = 2;
    const std::size_t nrOfBasisFunctions = 3;
    const std::size_t basisFunctionOrder = 1;
    const std::size_t elementIdLeft = 0;
    const std::size_t elementIdRight = 1;
    const std::size_t elementLeftFaceId = 3;
    const std::size_t elementRightFaceId = 0;
    const std::size_t faceId = 0;
    const std::size_t nrOfFaceMatrices = 1;
    
    std::vector<std::size_t> pointIndicesLeft;
    pointIndicesLeft.push_back(0);
    pointIndicesLeft.push_back(1);
    pointIndicesLeft.push_back(2);
    pointIndicesLeft.push_back(3);
    
    std::vector<std::size_t> pointIndicesRight;
    pointIndicesRight.push_back(2);
    pointIndicesRight.push_back(3);
    pointIndicesRight.push_back(4);
    pointIndicesRight.push_back(5);
    
    LinearAlgebra::NumericalVector coords0(dimension);
    coords0(0) = 0;
    coords0(1) = 0;
    LinearAlgebra::NumericalVector coords1(dimension);
    coords1(0) = 0;
    coords1(1) = 1;
    LinearAlgebra::NumericalVector coords2(dimension);
    coords2(0) = 1;
    coords2(1) = 0;
    LinearAlgebra::NumericalVector coords3(dimension);
    coords3(0) = 1;
    coords3(1) = 1;
    LinearAlgebra::NumericalVector coords4(dimension);
    coords4(0) = 2;
    coords4(1) = 0;
    LinearAlgebra::NumericalVector coords5(dimension);
    coords5(0) = 2;
    coords5(1) = 1;
    
    Geometry::PointPhysical point0(coords0);
    Geometry::PointPhysical point1(coords1);
    Geometry::PointPhysical point2(coords2);
    Geometry::PointPhysical point3(coords3);
    Geometry::PointPhysical point4(coords4);
    Geometry::PointPhysical point5(coords5);
    
    std::vector<Geometry::PointPhysical> pointsPhysical;
    pointsPhysical.push_back(point0);
    pointsPhysical.push_back(point1);
    pointsPhysical.push_back(point2);
    pointsPhysical.push_back(point3);
    pointsPhysical.push_back(point4);
    pointsPhysical.push_back(point5);
    
    Base::BasisFunctionSet basisFunctionSet(basisFunctionOrder);
    AssembleBasisFunctionSet_2D_Ord1_A0(basisFunctionSet);
    std::vector<const Base::BasisFunctionSet *> basisFunctionSetVector(1);
    basisFunctionSetVector[0] = &basisFunctionSet;
    std::vector<const Base::BasisFunctionSet *> *pBasisFunctionSetVector = &basisFunctionSetVector;
    
    std::cout << "Build elements.\n";
    Base::Element elementLeft(pointIndicesLeft, pBasisFunctionSetVector, pointsPhysical, nrOfUnknowns, nrOfTimeLevels, nrOfBasisFunctions, elementIdLeft);
    Base::Element elementRight(pointIndicesRight, pBasisFunctionSetVector, pointsPhysical, nrOfUnknowns, nrOfTimeLevels, nrOfBasisFunctions, elementIdRight);
    
    std::cout << "Build notes.\n";
    Base::Node node0(0);
    node0.addElement(&elementLeft, 0);
    
    Base::Node node1(1);
    node1.addElement(&elementLeft, 1);
    
    Base::Node node2(2);
    node2.addElement(&elementLeft, 2);
    node2.addElement(&elementRight,0);
    
    Base::Node node3(3);
    node3.addElement(&elementLeft, 3);
    node3.addElement(&elementRight,1);
    
    Base::Node node4(4);
    node4.addElement(&elementRight,2);
    
    Base::Node node5(5);
    node5.addElement(&elementRight,3);
    
    std::cout << "Build face.\n";
    Base::Face face(&elementLeft, elementLeftFaceId, &elementRight, elementRightFaceId, faceId, nrOfFaceMatrices);
    
    
    // Declare indices.
    std::size_t i;
    std::size_t j;
    std::size_t iVB;
    std::size_t iV;
    std::size_t iB;
    std::size_t iTimeLevel = 0;
    
    
    std::cout << "Do tests for the elements.\n";
    // Test ElementData::convertToSingleIndex.
    /* Every pair (iV,iB) should be mapped to a unique index iVB. Here iV = 0 .. nrOfUnknowns - 1 is the index corresponding to the variable and iB = 0 .. nrOfBasisFunctions - 1 is the index corresponding to the basis function.
     */
    std::vector<bool> checkMappingElementIndex(nrOfBasisFunctions * nrOfUnknowns, false);
    for(iV = 0; iV < nrOfUnknowns; iV++)
    {
        for(iB = 0; iB < nrOfBasisFunctions; iB++)
        {
            iVB = elementLeft.convertToSingleIndex(iB,iV);
            std::cout << "iV: " << iV << " ,iB: " << iB << " ,iVB: " << iVB << "\n";
            checkMappingElementIndex[iVB] = true;
        }
    }
    for(iVB = 0; iVB < nrOfBasisFunctions * nrOfUnknowns; iVB++)
    {
        logger.assert_always(checkMappingElementIndex[iVB], "Element mapping index check failed: %", iVB);
    }
    
    
    // Test Element::getSolution.
    LinearAlgebra::NumericalVector expansionCoefficients(nrOfBasisFunctions * nrOfUnknowns);
    for(iV = 0; iV < nrOfUnknowns; iV++)
    {
        for(iB = 0; iB < nrOfBasisFunctions; iB++)
        {
            iVB = elementLeft.convertToSingleIndex(iB,iV);
            expansionCoefficients(iVB) = iV + 1;   // Some arbitrary value
        }
    }
    elementRight.setTimeLevelDataVector(iTimeLevel, expansionCoefficients);
    LinearAlgebra::NumericalVector testVector = elementRight.getTimeLevelDataVector(iTimeLevel);
    logger.assert_always(testVector == expansionCoefficients, "Expansion coefficients incorrect: % != %", testVector, expansionCoefficients);
    
    Geometry::PointReference pointReference(coords0);
    LinearAlgebra::NumericalVector solutionVector(nrOfUnknowns);
    elementRight.getSolution(iTimeLevel, pointReference, solutionVector);
    for(iV = 0; iV < nrOfUnknowns; iV++)
    {
        logger.assert_always(solutionVector(iV) == (iV + 1) * solutionVector(0), "Solution vector test failed (%): % != (% + 1) * %",
                                     iV, solutionVector(iV), iV, solutionVector(0));
    }
    
    
    std::cout << "Do tests for the face.\n";
    // Test Face::convertToSingleIndex.
    /* Every triple (iS,iV,iB) should be mapped to a unique index i. Here iS = LEFT,RIGHT is the index corresponding to the side of the face, iV = 0 .. nrOfUnknowns - 1 is the index corresponding to the variable and iB = 0 .. nrOfBasisFunctions - 1 is the index corresponding to the basis function.
     */
    std::vector<bool> checkMappingFaceIndex(nrOfBasisFunctions * nrOfUnknowns * 2, false);
    std::vector<Base::Side> allSides;
    allSides.push_back(Base::Side::LEFT);
    allSides.push_back(Base::Side::RIGHT);
    for(Base::Side iS : allSides)
    {
        for(iV = 0; iV < nrOfUnknowns; iV++)
        {
            for(iB = 0; iB < nrOfBasisFunctions; iB++)
            {
                i = face.convertToSingleIndex(iS,iB,iV);
                if(iS == Base::Side::LEFT)
                    std::cout << "iS: LEFT" ;
                else
                    std::cout << "iS: RIGHT" ;
                std::cout << " ,iV: " << iV << " ,iB: " << iB << " ,i: " << i << "\n";
                checkMappingFaceIndex[i] = true;
            }
        }
    }
    for(i = 0; i < nrOfBasisFunctions * nrOfUnknowns * 2; i++)
    {
        logger.assert_always(checkMappingFaceIndex[i], "Face mapping index failed: %", i);
    }
    
    std::cout << "Finished tests.\n";
    return 0;
}