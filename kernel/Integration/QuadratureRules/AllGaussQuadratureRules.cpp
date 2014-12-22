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

#include "AllGaussQuadratureRules.hpp"
#include "GaussQuadratureRulesForCube.hpp"
#include "GaussQuadratureRulesForHypercube.hpp"
#include "GaussQuadratureRulesForLine.hpp"
#include "GaussQuadratureRulesForPoint.hpp"
#include "GaussQuadratureRulesForPyramid.hpp"
#include "GaussQuadratureRulesForSquare.hpp"
#include "GaussQuadratureRulesForTetrahedron.hpp"
#include "GaussQuadratureRulesForTriangle.hpp"
#include "GaussQuadratureRulesForTriangularPrism.hpp"

namespace QuadratureRules
{
    //add all the rules here
    AllGaussQuadratureRules::AllGaussQuadratureRules()
    {
        //*************************POINT QUADRATURES****************************
        addRule(&Cn0_inf_1::Instance());
        //*************************LINE QUADRATURES*****************************
        addRule(&Cn1_1_1::Instance());
        addRule(&Cn1_3_4::Instance());
        addRule(&Cn1_5_9::Instance());
        addRule(&C1_7_x::Instance());
        addRule(&C1_9_25::Instance());
        addRule(&C1_11_36::Instance());
        //*************************SQUARE QUADRATURES*****************************
        addRule(&Cn2_1_1::Instance());
        addRule(&Cn2_3_4::Instance());
        addRule(&Cn2_5_9::Instance());
        addRule(&C2_7_4::Instance());
        addRule(&C2_9_5::Instance());
        addRule(&C2_11_6::Instance());
        //*************************TRIANGLE QUADRATURES*****************************
        addRule(&Tn2_1_1::Instance());
        addRule(&Tn2_2_1::Instance());
        addRule(&Tn2_3_1::Instance());
        addRule(&Tn2_4_1::Instance());
        addRule(&T2_5_1::Instance());
        addRule(&T2_6_1::Instance());
        addRule(&T2_7_1::Instance());
        addRule(&T2_8_1::Instance());
        addRule(&T2_9_1::Instance());
        addRule(&T2_10_1::Instance());
        addRule(&T2_11_1::Instance());
        //*************************CUBE QUADRATURES*****************************
        addRule(&Cn3_1_1::Instance());
        addRule(&Cn3_3_4::Instance());
        addRule(&Cn3_5_9::Instance());
        addRule(&C3_7_2::Instance());
        addRule(&C3_9_2::Instance());
        addRule(&C3_11_2::Instance());
        //*************************PYRAMID QUADRATURES*****************************
        addRule(&Pyramid_1_1::Instance());
        addRule(&Pyramid_3_1::Instance());
        addRule(&Pyramid_5_1::Instance());
        addRule(&Pyramid_7_1::Instance());
        //*************************TETRAHEDRON QUADRATURES*****************************
        addRule(&Tn3_1_1::Instance());
        addRule(&Tn3_2_1::Instance());
        addRule(&Tn3_3_1::Instance());
        addRule(&Tn3_4_1::Instance());
        addRule(&T3_5_1::Instance());
        addRule(&T3_6_1::Instance());
        addRule(&T3_7_1::Instance());
        addRule(&T3_8_1::Instance());
        addRule(&T3_9_1::Instance());
        addRule(&T3_10_1::Instance());
        //*************************TRIANGULARPRISM QUADRATURES*****************************
        addRule(&TriPrism_1_1::Instance());
        addRule(&TriPrism_3_1::Instance());
        addRule(&TriPrism_5_1::Instance());
        addRule(&TriPrism_7_1::Instance());
        //*************************HYPERCUBE QUADRATURES*****************************
        addRule(&Cn4_1_1::Instance());
        addRule(&Cn4_3_4::Instance());
    }
    
    AllGaussQuadratureRules& AllGaussQuadratureRules::instance()
    {
        static AllGaussQuadratureRules theInstance;
        return theInstance;
    }
    
    void AllGaussQuadratureRules::addRule(const GaussQuadratureRule* rule)
    {
        std::vector<const GaussQuadratureRule*>& listForGeometry = listOfRules_[rule->forReferenceGeometry()];
        std::vector<const GaussQuadratureRule*>::iterator it = listForGeometry.begin();
        while (it != listForGeometry.end())
        {
            if ((*it)->order() < rule->order())
                ++it;
            else
                break;
        }
        listForGeometry.insert(it, rule);
    }
    
    const GaussQuadratureRule* AllGaussQuadratureRules::getRule(const Geometry::ReferenceGeometry* referenceGeometry, int order)
    {
        for (const GaussQuadratureRule* rule : listOfRules_[referenceGeometry])
        {
            if (rule->order() >= order)
            {
                return rule;
            }
        }
        throw "Tried to find a quadrature rule but didn't find one";
    }

}
