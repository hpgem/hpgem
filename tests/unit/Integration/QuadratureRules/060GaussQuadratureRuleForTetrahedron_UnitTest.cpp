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

//naming convention: <Digit><ClassName>_UnitTest.cpp where <Digit> is a number that will make sure
//the unit tests are ordered such that the first failing unit test indicate the culprit class and
//other 'unit' tests may assume correct execution of all prior unit tests
#include "Integration/QuadratureRules/GaussQuadratureRulesForTetrahedron.h"
#include "Logger.h"
#include <typeinfo>

#include "Utilities/BasisFunctions3DH1ConformingTetrahedron.h"
#include "Geometry/ReferenceTetrahedron.h"
#include "Base/BasisFunctionSet.h"
#include "Geometry/PointReference.h"
#include "LinearAlgebra/NumericalVector.h"
#include <cmath>

void testRule(QuadratureRules::GaussQuadratureRule& test, std::size_t expectedOrder)
{
    std::cout << test.getName() << std::endl;
    logger.assert_always((test.dimension() == 3), "dimension");
    logger.assert_always((test.order() >= expectedOrder), "order");
    logger.assert_always((typeid(*test.forReferenceGeometry()) == typeid(Geometry::ReferenceTetrahedron)), "forReferenceGeometry");
    Geometry::PointReference point(3);
    std::cout.precision(14);
    Base::BasisFunctionSet* functions = Utilities::createDGBasisFunctionSet3DH1Tetrahedron(expectedOrder);
    for (std::size_t i = 0; i < functions->size(); ++i)
    {
        double integrated = 0;
        for (std::size_t j = 0; j < test.nrOfPoints(); ++j)
        {
            point = test.getPoint(j);
            integrated += test.weight(j) * functions->eval(i, point);
        }
        if (i < 4)
        {
            logger.assert_always((std::abs(integrated - 1. / 24.) < 1e-12), "integration");
        }
        else if (i < 10)
        {
            logger.assert_always((std::abs(integrated + 0.102062072616 / 5.) < 1e-12), "integration");
        }
        else if (15 < i && i < 20)
        {
            logger.assert_always((std::abs(integrated - 1. / 120.) < 1e-12), "integration");
        }
        else if (19 < i && i < 26)
        {
            logger.assert_always((std::abs(integrated - 0.00408315786255) < 1e-12), "integration");
        }
        else if (i == 34)
        {
            logger.assert_always((std::abs(integrated + 0.0029160592176) < 1e-12), "integration");
        }
        else if (i == 41 || i == 43 || i == 44 || i == 46 || i == 47 || i == 49 || i == 50 || i == 52)
        {
            logger.assert_always((std::abs(integrated + 0.00204579272096) < 1e-12), "integration");
        }
        else if (i == 42 || i == 45 || i == 48 || i == 51)
        {
            logger.assert_always((std::abs(integrated + 0.00049603174603) < 1e-12), "integration");
        }
        else if (55 < i && i < 62)
        {
            logger.assert_always((std::abs(integrated + 0.00081430829164) < 1e-12), "integration");
        }
        else if (i == 63 || i == 67 || i == 71 || i == 75)
        {
            logger.assert_always((std::abs(integrated - 0.00016303129914) < 1e-12), "integration");
        }
        else if (i == 64 || i == 68 || i == 72 || i == 76)
        {
            logger.assert_always((std::abs(integrated + 0.00016303129914) < 1e-12), "integration");
        }
        else if (i == 78 || i == 80 || i == 83)
        {
            logger.assert_always((std::abs(integrated - 0.00080425836687) < 1e-12), "integration");
        }
        else if (i == 79 || i == 81)
        {
            logger.assert_always((std::abs(integrated - 0.000135002741555) < 1e-12), "integration");
        }
        else if (i == 90 || i == 94 || i == 95 || i == 99 || i == 100 || i == 104 || i == 105 || i == 109)
        {
            logger.assert_always((std::abs(integrated - 0.000598391942328) < 1e-12), "integration");
        }
        else if (i == 91 || i == 93 || i == 96 || i == 98 || i == 101 || i == 103 || i == 106 || i == 108)
        {
            logger.assert_always((std::abs(integrated - 0.00026619856875) < 1e-12), "integration");
        }
        else if (i == 92 || i == 97 || i == 102 || i == 107)
        {
            logger.assert_always((std::abs(integrated - 0.000564236111111) < 1e-12), "integration");
        }
        else if (i == 111 || i == 117)
        {
            logger.assert_always((std::abs(integrated + 0.000039934349489) < 1e-12), "integration");
        }
        else if (i == 112 || i == 114)
        {
            logger.assert_always((std::abs(integrated - 0.000039934349489) < 1e-12), "integration");
        }
        else if (119 < i && i < 126)
        {
            logger.assert_always((std::abs(integrated - 0.000152824374304) < 1e-12), "integration");
        }
        else if (i == 127 || i == 133 || i == 139 || i == 145)
        {
            logger.assert_always((std::abs(integrated + 0.000145660464112) < 1e-12), "integration");
        }
        else if (i == 128 || i == 134 || i == 140 || i == 146)
        {
            logger.assert_always((std::abs(integrated - 0.000038774764434) < 1e-12), "integration");
        }
        else if (i == 129 || i == 135 || i == 141 || i == 147)
        {
            logger.assert_always((std::abs(integrated + 0.000038774764434) < 1e-12), "integration");
        }
        else if (i == 130 || i == 136 || i == 142 || i == 148)
        {
            logger.assert_always((std::abs(integrated - 0.000145660464112) < 1e-12), "integration");
        }
        else if (i == 150 || i == 154 || i == 164)
        {
            logger.assert_always((std::abs(integrated + 0.000288356788986) < 1e-12), "integration");
        }
        else if (i == 151 || i == 153 || i == 155 || i == 162)
        {
            logger.assert_always((std::abs(integrated + 0.000083976221839) < 1e-12), "integration");
        }
        else if (i == 152 || i == 159)
        {
            logger.assert_always((std::abs(integrated + 0.000237328683166) < 1e-12), "integration");
        }
        else if (i == 156)
        {
            logger.assert_always((std::abs(integrated + 0.000009373640639) < 1e-12), "integration");
        }
        else if (i == 157 || i == 160)
        {
            logger.assert_always((std::abs(integrated + 0.000042181382877) < 1e-12), "integration");
        }
        else if (i == 161)
        {
            logger.assert_always((std::abs(integrated + 0.000210481547061) < 1e-12), "integration");
        }
        else if (i == 171 || i == 177 || i == 178 || i == 184 || i == 185 || i == 191 || i == 192 || i == 198)
        {
            logger.assert_always((std::abs(integrated + 0.000166374105469) < 1e-12), "integration");
        }
        else if (i == 172 || i == 176 || i == 179 || i == 183 || i == 186 || i == 190 || i == 193 || i == 197)
        {
            logger.assert_always((std::abs(integrated + 0.000559879010298) < 1e-12), "integration");
        }
        else if (i == 173 || i == 175 || i == 180 || i == 182 || i == 187 || i == 189 || i == 194 || i == 196)
        {
            logger.assert_always((std::abs(integrated + 0.000181525744621) < 1e-12), "integration");
        }
        else if (i == 174 || i == 181 || i == 188 || i == 195)
        {
            logger.assert_always((std::abs(integrated + 0.000174343885281) < 1e-12), "integration");
        }
        else if (i == 200 || i == 217)
        {
            logger.assert_always((std::abs(integrated - 0.0000414135675539) < 1e-12), "integration");
        }
        else if (i == 201 || i == 214)
        {
            logger.assert_always((std::abs(integrated + 0.0000146120596555) < 1e-12), "integration");
        }
        else if (i == 202 || i == 210)
        {
            logger.assert_always((std::abs(integrated - 0.0000146120596555) < 1e-12), "integration");
        }
        else if (i == 203 || i == 205)
        {
            logger.assert_always((std::abs(integrated + 0.0000414135675539) < 1e-12), "integration");
        }
        else if (i == 207)
        {
            logger.assert_always((std::abs(integrated + 0.0000152501816173) < 1e-12), "integration");
        }
        else if (i == 211)
        {
            logger.assert_always((std::abs(integrated - 0.0000152501816173) < 1e-12), "integration");
        }
        else if (219 < i && i < 226)
        {
            logger.assert_always((std::abs(integrated + 0.0000445921151835) < 1e-12), "integration");
        }
        else if (i == 227 || i == 235 || i == 243 || i == 251)
        {
            logger.assert_always((std::abs(integrated - 0.0000805455174651) < 1e-12), "integration");
        }
        else if (i == 228 || i == 236 || i == 244 || i == 252)
        {
            logger.assert_always((std::abs(integrated + 0.0001965625984962) < 1e-12), "integration");
        }
        else if (i == 229 || i == 237 || i == 245 || i == 253)
        {
            logger.assert_always((std::abs(integrated - 0.0000535536976881) < 1e-12), "integration");
        }
        else if (i == 230 || i == 238 || i == 246 || i == 254)
        {
            logger.assert_always((std::abs(integrated + 0.0000535536976881) < 1e-12), "integration");
        }
        else if (i == 231 || i == 239 || i == 247 || i == 255)
        {
            logger.assert_always((std::abs(integrated - 0.0001965625984962) < 1e-12), "integration");
        }
        else if (i == 232 || i == 240 || i == 248 || i == 256)
        {
            logger.assert_always((std::abs(integrated + 0.0000805455174651) < 1e-12), "integration");
        }
        else if (i == 258 || i == 264 || i == 285)
        {
            logger.assert_always((std::abs(integrated - 0.0001018829162026) < 1e-12), "integration");
        }
        else if (i == 259 || i == 263 || i == 265 || i == 283)
        {
            logger.assert_always((std::abs(integrated - 0.0001494113288829) < 1e-12), "integration");
        }
        else if (i == 260 || i == 262 || i == 271 || i == 280)
        {
            logger.assert_always((std::abs(integrated - 0.0000915303367133) < 1e-12), "integration");
        }
        else if (i == 261 || i == 276)
        {
            logger.assert_always((std::abs(integrated - 0.0000588501024355) < 1e-12), "integration");
        }
        else if (i == 266)
        {
            logger.assert_always((std::abs(integrated - 0.0000107562089804) < 1e-12), "integration");
        }
        else if (i == 267 || i == 272)
        {
            logger.assert_always((std::abs(integrated - 0.0000264460001519) < 1e-12), "integration");
        }
        else if (i == 268 || i == 277)
        {
            logger.assert_always((std::abs(integrated - 0.0000068926475018) < 1e-12), "integration");
        }
        else if (i == 269 || i == 281)
        {
            logger.assert_always((std::abs(integrated - 0.0000187103803273) < 1e-12), "integration");
        }
        else if (i == 273)
        {
            logger.assert_always((std::abs(integrated - 0.0000690224423240) < 1e-12), "integration");
        }
        else if (i == 274 || i == 278)
        {
            logger.assert_always((std::abs(integrated - 0.0000251763019278) < 1e-12), "integration");
        }
        else if (i == 275 || i == 282)
        {
            logger.assert_always((std::abs(integrated - 0.0000672472006580) < 1e-12), "integration");
        }
        else
        {
            logger.assert_always((std::abs(integrated) < 1e-12), "integration");
        }
        
    }
    
    delete functions;
}

int main()
{
    
    testRule(QuadratureRules::Tn3_1_1::Instance(), 1);
    testRule(QuadratureRules::Tn3_2_4::Instance(), 2);
    testRule(QuadratureRules::Tn3_3_5::Instance(), 3);
    testRule(QuadratureRules::Tn3_4_11::Instance(), 4);
    testRule(QuadratureRules::T3_5_14::Instance(), 5);
    testRule(QuadratureRules::T3_6_24::Instance(), 6);
    testRule(QuadratureRules::T3_7_31::Instance(), 7);
    testRule(QuadratureRules::T3_8_43::Instance(), 8);
    testRule(QuadratureRules::T3_9_53::Instance(), 9);
    testRule(QuadratureRules::T3_10_126::Instance(), 10); ///\TODO implement 11th order quadrature
            
    return 0;
}

