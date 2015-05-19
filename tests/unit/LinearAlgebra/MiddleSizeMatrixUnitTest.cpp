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

#include <iostream>
#include <cstdlib>
#include <cmath>
#include "LinearAlgebra/MiddleSizeMatrix.h"
#include "LinearAlgebra/MiddleSizeVector.h"
#include "Logger.h"

using LinearAlgebra::MiddleSizeMatrix;
using LinearAlgebra::MiddleSizeVector;

int main(int argc, char** argv)
{
    //constructors
    LinearAlgebra::MiddleSizeVector vec0{0., 1.}, vec1{2., 3.}, vec2{4., 5.}, vec3{6., 7.};
    MiddleSizeMatrix A0, A22(2, 2), A23(2, 3), A32(3, 2), destroy(3, 3, 1), count0({vec0, vec1}), count1({vec2, vec3}), copy(count0), bla({count0}), merge({count0, count1});
    logger.assert_always(destroy.getNRows() == 3, "Rows in a matrix");
    logger.assert_always(destroy.getNCols() == 3, "Columns in a matrix");
    logger.assert_always(destroy.size() == 9, "Size of a matrix");
    for(std::size_t i = 0; i < destroy.size(); ++i)
    {
        logger.assert_always(std::abs(destroy[i] - 1.) < 1e-12, "Entry of a matrix");
    }
    for(std::size_t i = 0; i < destroy.getNRows(); ++i)
    {
        for(std::size_t j = 0; j < destroy.getNCols(); ++j)
        {
            logger.assert_always(std::abs(destroy(i, j) - 1.) < 1e-12, "Entry of a matrix");
        }
    }
    MiddleSizeMatrix moved(std::move(destroy));
    logger.assert_always(moved.getNRows() == 3, "Rows in a matrix");
    logger.assert_always(moved.getNCols() == 3, "Columns in a matrix");
    logger.assert_always(moved.size() == 9, "Size of a matrix");
    for(std::size_t i = 0; i < moved.size(); ++i)
    {
        logger.assert_always(std::abs(moved[i] - 1.) < 1e-12, "Entry of a matrix");
    }
    for(std::size_t i = 0; i < moved.getNRows(); ++i)
    {
        for(std::size_t j = 0; j < moved.getNCols(); ++j)
        {
            logger.assert_always(std::abs(moved(i, j) - 1.) < 1e-12, "Entry of a matrix");
        }
    }
    logger.assert_always(A0.getNRows() == 0, "Rows in a matrix");
    logger.assert_always(A0.getNCols() == 0, "Columns in a matrix");
    logger.assert_always(A0.size() == 0, "Size of a matrix");
    logger.assert_always(A22.getNRows() == 2, "Rows in a matrix");
    logger.assert_always(A22.getNCols() == 2, "Columns in a matrix");
    logger.assert_always(A22.size() == 4, "Size of a matrix");
    logger.assert_always(A23.getNRows() == 2, "Rows in a matrix");
    logger.assert_always(A23.getNCols() == 3, "Columns in a matrix");
    logger.assert_always(A23.size() == 6, "Size of a matrix");
    logger.assert_always(A32.getNRows() == 3, "Rows in a matrix");
    logger.assert_always(A32.getNCols() == 2, "Columns in a matrix");
    logger.assert_always(A32.size() == 6, "Size of a matrix");
    logger.assert_always(count0.getNRows() == 2, "Rows in a matrix");
    logger.assert_always(count0.getNCols() == 2, "Columns in a matrix");
    logger.assert_always(count0.size() == 4, "Size of a matrix");
    for(std::size_t i = 0; i < count0.size(); ++i)
    {
        logger.assert_always(std::abs(count0[i] - i) < 1e-12, "Entry of a matrix");
    }
    logger.assert_always(std::abs(count0(0, 0) - 0.) < 1e-12, "Entry of a matrix");
    logger.assert_always(std::abs(count0(1, 0) - 1.) < 1e-12, "Entry of a matrix");
    logger.assert_always(std::abs(count0(0, 1) - 2.) < 1e-12, "Entry of a matrix");
    logger.assert_always(std::abs(count0(1, 1) - 3.) < 1e-12, "Entry of a matrix");
    logger.assert_always(count1.getNRows() == 2, "Rows in a matrix");
    logger.assert_always(count1.getNCols() == 2, "Columns in a matrix");
    logger.assert_always(count1.size() == 4, "Size of a matrix");
    for(std::size_t i = 0; i < count1.size(); ++i)
    {
        logger.assert_always(std::abs(count1[i] - 4. - i) < 1e-12, "Entry of a matrix");
    }
    logger.assert_always(copy.getNRows() == 2, "Rows in a matrix");
    logger.assert_always(copy.getNCols() == 2, "Columns in a matrix");
    logger.assert_always(copy.size() == 4, "Size of a matrix");
    for(std::size_t i = 0; i < copy.size(); ++i)
    {
        logger.assert_always(std::abs(copy[i] - i) < 1e-12, "Entry of a matrix");
    }
    logger.assert_always(bla.getNRows() == 2, "Rows in a matrix");
    logger.assert_always(bla.getNCols() == 2, "Columns in a matrix");
    logger.assert_always(bla.size() == 4, "Size of a matrix");
    for(std::size_t i = 0; i < bla.size(); ++i)
    {
        logger.assert_always(std::abs(bla[i] - i) < 1e-12, "Entry of a matrix");
    }
    logger.assert_always(merge.getNRows() == 2, "Rows in a matrix");
    logger.assert_always(merge.getNCols() == 4, "Columns in a matrix");
    logger.assert_always(merge.size() == 8, "Size of a matrix");
    for(std::size_t i = 0; i < merge.size(); ++i)
    {
        logger.assert_always(std::abs(merge[i] - i) < 1e-12, "Entry of a matrix");
    }
    A23(0, 0) = 0.0;
    A23(1, 0) = 0.1;
    A23(0, 1) = 0.2;
    A23(1, 1) = 0.3;
    A23(0, 2) = 0.4;
    A23(1, 2) = 0.5;
    A32[0] = 0.0;
    A32[1] = 0.1;
    A32[2] = 0.2;
    A32[3] = 0.3;
    A32[4] = 0.4;
    A32[5] = 0.5;

    //out-of-place operators
    logger.assert_always(std::abs((count0*count1)(0, 0) - 10.) < 1e-12, "multiply");
    logger.assert_always(std::abs((count0*count1)(1, 0) - 19.) < 1e-12, "multiply");
    logger.assert_always(std::abs((count0*count1)(0, 1) - 14.) < 1e-12, "multiply");
    logger.assert_always(std::abs((count0*count1)(1, 1) - 27.) < 1e-12, "multiply");
    logger.assert_always(std::abs((count1*count0)(0, 0) - 6.) < 1e-12, "multiply");
    logger.assert_always(std::abs((count1*count0)(1, 0) - 7.) < 1e-12, "multiply");
    logger.assert_always(std::abs((count1*count0)(0, 1) - 26.) < 1e-12, "multiply");
    logger.assert_always(std::abs((count1*count0)(1, 1) - 31.) < 1e-12, "multiply");
    logger.assert_always(std::abs((count0*A23)(0, 0) - .2) < 1e-12, "multiply");
    logger.assert_always(std::abs((count0*A23)(1, 0) - .3) < 1e-12, "multiply");
    logger.assert_always(std::abs((count0*A23)(0, 1) - .6) < 1e-12, "multiply");
    logger.assert_always(std::abs((count0*A23)(1, 1) - 1.1) < 1e-12, "multiply");
    logger.assert_always(std::abs((count0*A23)(0, 2) - 1.) < 1e-12, "multiply");
    logger.assert_always(std::abs((count0*A23)(1, 2) - 1.9) < 1e-12, "multiply");
    logger.assert_always(std::abs((A32*count0)(0, 0) - .3) < 1e-12, "multiply");
    logger.assert_always(std::abs((A32*count0)(1, 0) - .4) < 1e-12, "multiply");
    logger.assert_always(std::abs((A32*count0)(2, 0) - .5) < 1e-12, "multiply");
    logger.assert_always(std::abs((A32*count0)(0, 1) - .9) < 1e-12, "multiply");
    logger.assert_always(std::abs((A32*count0)(1, 1) - 1.4) < 1e-12, "multiply");
    logger.assert_always(std::abs((A32*count0)(2, 1) - 1.9) < 1e-12, "multiply");
    logger.assert_always(std::abs((A32*A23)(0, 0) - .03) < 1e-12, "multiply");
    logger.assert_always(std::abs((A32*A23)(1, 0) - .04) < 1e-12, "multiply");
    logger.assert_always(std::abs((A32*A23)(2, 0) - .05) < 1e-12, "multiply");
    logger.assert_always(std::abs((A32*A23)(0, 1) - .09) < 1e-12, "multiply");
    logger.assert_always(std::abs((A32*A23)(1, 1) - .14) < 1e-12, "multiply");
    logger.assert_always(std::abs((A32*A23)(2, 1) - .19) < 1e-12, "multiply");
    logger.assert_always(std::abs((A32*A23)(0, 2) - .15) < 1e-12, "multiply");
    logger.assert_always(std::abs((A32*A23)(1, 2) - .24) < 1e-12, "multiply");
    logger.assert_always(std::abs((A32*A23)(2, 2) - .33) < 1e-12, "multiply");
    logger.assert_always(std::abs((A23*A32)(0, 0) - .1) < 1e-12, "multiply");
    logger.assert_always(std::abs((A23*A32)(1, 0) - .13) < 1e-12, "multiply");
    logger.assert_always(std::abs((A23*A32)(0, 1) - .28) < 1e-12, "multiply");
    logger.assert_always(std::abs((A23*A32)(1, 1) - .40) < 1e-12, "multiply");
    logger.assert_always(std::abs((vec1*count0)*vec1 - 45) < 1e-12, "multiply");
    logger.assert_always(std::abs(vec1*(count0*vec1) - 45) < 1e-12, "multiply");
    MiddleSizeVector size3 = {3., 4., 5.};
    logger.assert_always(std::abs(size3*(A32*vec0) - 5) < 1e-12, "multiply");
    logger.assert_always(std::abs((size3*A32)*vec0 - 5) < 1e-12, "multiply");
    logger.assert_always(std::abs((2*count0*2)(0, 0) - 0.) < 1e-12, "multiply");
    logger.assert_always(std::abs((2*count0*2)(1, 0) - 4.) < 1e-12, "multiply");
    logger.assert_always(std::abs((2*count0*2)(0, 1) - 8.) < 1e-12, "multiply");
    logger.assert_always(std::abs((2*count0*2)(1, 1) - 12.) < 1e-12, "multiply");
    logger.assert_always(std::abs((count0/2.)(0, 0) - 0.) < 1e-12, "divide");
    logger.assert_always(std::abs((count0/2.)(1, 0) - .5) < 1e-12, "divide");
    logger.assert_always(std::abs((count0/2.)(0, 1) - 1.) < 1e-12, "divide");
    logger.assert_always(std::abs((count0/2.)(1, 1) - 1.5) < 1e-12, "divide");
    logger.assert_always(std::abs((count0+count1)(0, 0) - 4.) < 1e-12, "add");
    logger.assert_always(std::abs((count0+count1)(1, 0) - 6.) < 1e-12, "add");
    logger.assert_always(std::abs((count0+count1)(0, 1) - 8.) < 1e-12, "add");
    logger.assert_always(std::abs((count0+count1)(1, 1) - 10.) < 1e-12, "add");
    logger.assert_always(std::abs((count1+count0)(0, 0) - 4.) < 1e-12, "add");
    logger.assert_always(std::abs((count1+count0)(1, 0) - 6.) < 1e-12, "add");
    logger.assert_always(std::abs((count1+count0)(0, 1) - 8.) < 1e-12, "add");
    logger.assert_always(std::abs((count1+count0)(1, 1) - 10.) < 1e-12, "add");
    logger.assert_always(std::abs((count0-count1)(0, 0) + 4.) < 1e-12, "subtract");
    logger.assert_always(std::abs((count0-count1)(1, 0) + 4.) < 1e-12, "subtract");
    logger.assert_always(std::abs((count0-count1)(0, 1) + 4.) < 1e-12, "subtract");
    logger.assert_always(std::abs((count0-count1)(1, 1) + 4.) < 1e-12, "subtract");
    logger.assert_always(std::abs((count1-count0)(0, 0) - 4.) < 1e-12, "subtract");
    logger.assert_always(std::abs((count1-count0)(1, 0) - 4.) < 1e-12, "subtract");
    logger.assert_always(std::abs((count1-count0)(0, 1) - 4.) < 1e-12, "subtract");
    logger.assert_always(std::abs((count1-count0)(1, 1) - 4.) < 1e-12, "subtract");

    //assignent operators
    destroy = 4;
    logger.assert_always(destroy.getNRows() == 1, "Rows in a matrix");
    logger.assert_always(destroy.getNCols() == 1, "Columns in a matrix");
    logger.assert_always(destroy.size() == 1, "Size of a matrix");
    logger.assert_always(std::abs(destroy[0] - 4.) < 1e-12, "Entry of a matrix");
    MiddleSizeMatrix extra = copy;
    logger.assert_always(extra.getNRows() == 2, "Rows in a matrix");
    logger.assert_always(extra.getNCols() == 2, "Columns in a matrix");
    logger.assert_always(extra.size() == 4, "Size of a matrix");
    for(std::size_t i = 0; i < extra.size(); ++i)
    {
        logger.assert_always(std::abs(extra[i] - i) < 1e-12, "Entry of a matrix");
    }
    copy = count1;
    logger.assert_always(copy.getNRows() == 2, "Rows in a matrix");
    logger.assert_always(copy.getNCols() == 2, "Columns in a matrix");
    logger.assert_always(copy.size() == 4, "Size of a matrix");
    for(std::size_t i = 0; i < copy.size(); ++i)
    {
        logger.assert_always(std::abs(copy[i] - i - 4.) < 1e-12, "Entry of a matrix");
    }
    A0 = std::move(destroy);
    logger.assert_always(A0.getNRows() == 1, "Rows in a matrix");
    logger.assert_always(A0.getNCols() == 1, "Columns in a matrix");
    logger.assert_always(A0.size() == 1, "Size of a matrix");
    logger.assert_always(std::abs(A0[0] - 4.) < 1e-12, "Entry of a matrix");
    count0*=count1;
    logger.assert_always(count0.getNRows() == 2, "Rows in a matrix");
    logger.assert_always(count0.getNCols() == 2, "Columns in a matrix");
    logger.assert_always(count0.size() == 4, "Size of a matrix");
    logger.assert_always(std::abs((count0)(0, 0) - 10.) < 1e-12, "multiply");
    logger.assert_always(std::abs((count0)(1, 0) - 19.) < 1e-12, "multiply");
    logger.assert_always(std::abs((count0)(0, 1) - 14.) < 1e-12, "multiply");
    logger.assert_always(std::abs((count0)(1, 1) - 27.) < 1e-12, "multiply");
    count1*=A23;
    logger.assert_always(count1.getNRows() == 2, "Rows in a matrix");
    logger.assert_always(count1.getNCols() == 3, "Columns in a matrix");
    logger.assert_always(count1.size() == 6, "Size of a matrix");
    logger.assert_always(std::abs((count1)(0, 0) - .6) < 1e-12, "multiply");
    logger.assert_always(std::abs((count1)(1, 0) - .7) < 1e-12, "multiply");
    logger.assert_always(std::abs((count1)(0, 1) - 2.6) < 1e-12, "multiply");
    logger.assert_always(std::abs((count1)(1, 1) - 3.1) < 1e-12, "multiply");
    logger.assert_always(std::abs((count1)(0, 2) - 4.6) < 1e-12, "multiply");
    logger.assert_always(std::abs((count1)(1, 2) - 5.5) < 1e-12, "multiply");
    count0*=4;
    logger.assert_always(count0.getNRows() == 2, "Rows in a matrix");
    logger.assert_always(count0.getNCols() == 2, "Columns in a matrix");
    logger.assert_always(count0.size() == 4, "Size of a matrix");
    logger.assert_always(std::abs((count0)(0, 0) - 40.) < 1e-12, "multiply");
    logger.assert_always(std::abs((count0)(1, 0) - 76.) < 1e-12, "multiply");
    logger.assert_always(std::abs((count0)(0, 1) - 56.) < 1e-12, "multiply");
    logger.assert_always(std::abs((count0)(1, 1) - 108.) < 1e-12, "multiply");
    count0/=2;
    logger.assert_always(count0.getNRows() == 2, "Rows in a matrix");
    logger.assert_always(count0.getNCols() == 2, "Columns in a matrix");
    logger.assert_always(count0.size() == 4, "Size of a matrix");
    logger.assert_always(std::abs((count0)(0, 0) - 20.) < 1e-12, "divide");
    logger.assert_always(std::abs((count0)(1, 0) - 38.) < 1e-12, "divide");
    logger.assert_always(std::abs((count0)(0, 1) - 28.) < 1e-12, "divide");
    logger.assert_always(std::abs((count0)(1, 1) - 54.) < 1e-12, "divide");
    copy+=extra;
    logger.assert_always(copy.getNRows() == 2, "Rows in a matrix");
    logger.assert_always(copy.getNCols() == 2, "Columns in a matrix");
    logger.assert_always(copy.size() == 4, "Size of a matrix");
    logger.assert_always(std::abs((copy)(0, 0) - 4.) < 1e-12, "add");
    logger.assert_always(std::abs((copy)(1, 0) - 6.) < 1e-12, "add");
    logger.assert_always(std::abs((copy)(0, 1) - 8.) < 1e-12, "add");
    logger.assert_always(std::abs((copy)(1, 1) - 10.) < 1e-12, "add");
    extra-=count0;
    logger.assert_always(extra.getNRows() == 2, "Rows in a matrix");
    logger.assert_always(extra.getNCols() == 2, "Columns in a matrix");
    logger.assert_always(extra.size() == 4, "Size of a matrix");
    logger.assert_always(std::abs((extra)(0, 0) + 20.) < 1e-12, "subtract");
    logger.assert_always(std::abs((extra)(1, 0) + 37.) < 1e-12, "subtract");
    logger.assert_always(std::abs((extra)(0, 1) + 26.) < 1e-12, "subtract");
    logger.assert_always(std::abs((extra)(1, 1) + 51.) < 1e-12, "subtract");
    extra.axpy(3., copy);
    logger.assert_always(extra.getNRows() == 2, "Rows in a matrix");
    logger.assert_always(extra.getNCols() == 2, "Columns in a matrix");
    logger.assert_always(extra.size() == 4, "Size of a matrix");
    logger.assert_always(std::abs((extra)(0, 0) + 8.) < 1e-12, "ax+y");
    logger.assert_always(std::abs((extra)(1, 0) + 19.) < 1e-12, "ax+y");
    logger.assert_always(std::abs((extra)(0, 1) + 2.) < 1e-12, "ax+y");
    logger.assert_always(std::abs((extra)(1, 1) + 21.) < 1e-12, "ax+y");
    A0.resize(3, 7);
    logger.assert_always(A0.getNRows() == 3, "Rows in a matrix");
    logger.assert_always(A0.getNCols() == 7, "Columns in a matrix");
    logger.assert_always(A0.size() == 21, "Size of a matrix");
    logger.assert_always(std::abs(A0[0] - 4.) < 1e-12, "Entry of a matrix");

    //wedge stuff
    logger.assert_always(std::abs((MiddleSizeMatrix(vec0).computeWedgeStuffVector())*(MiddleSizeMatrix(vec0).computeWedgeStuffVector()) - vec0 * vec0) < 1e-12, "norm of wedge stuff vector");
    logger.assert_always(std::abs((MiddleSizeMatrix(vec0).computeWedgeStuffVector()) * vec0) < 1e-12, "direction of wedge stuff vector");
    logger.assert_always(std::abs((MiddleSizeMatrix(vec1).computeWedgeStuffVector())*(MiddleSizeMatrix(vec1).computeWedgeStuffVector()) - vec1 * vec1) < 1e-12, "norm of wedge stuff vector");
    logger.assert_always(std::abs((MiddleSizeMatrix(vec1).computeWedgeStuffVector()) * vec1) < 1e-12, "direction of wedge stuff vector");
    logger.assert_always(std::abs((MiddleSizeMatrix(vec2).computeWedgeStuffVector())*(MiddleSizeMatrix(vec2).computeWedgeStuffVector()) - vec2 * vec2) < 1e-12, "norm of wedge stuff vector");
    logger.assert_always(std::abs((MiddleSizeMatrix(vec2).computeWedgeStuffVector()) * vec2) < 1e-12, "direction of wedge stuff vector");
    logger.assert_always(std::abs((MiddleSizeMatrix(vec3).computeWedgeStuffVector())*(MiddleSizeMatrix(vec3).computeWedgeStuffVector()) - vec3 * vec3) < 1e-12, "norm of wedge stuff vector");
    logger.assert_always(std::abs((MiddleSizeMatrix(vec3).computeWedgeStuffVector()) * vec3) < 1e-12, "direction of wedge stuff vector");
    MiddleSizeVector vec3D0{0., 1., 2.}, vec3D1{3., 4., 5.}, vec3D2{0., -1., 2.};
    ///\todo test that the norm of the 3D wedge stuff vector equals the area of the triangle formed by nodes {0, 0, 0}, v1 and v2
    logger.assert_always(std::abs((MiddleSizeMatrix({vec3D0, vec3D1}).computeWedgeStuffVector()) * vec3D0) < 1e-12, "direction of wedge stuff vector");
    logger.assert_always(std::abs((MiddleSizeMatrix({vec3D0, vec3D1}).computeWedgeStuffVector()) * vec3D1) < 1e-12, "direction of wedge stuff vector");
    logger.assert_always(std::abs((MiddleSizeMatrix({vec3D1, vec3D2}).computeWedgeStuffVector()) * vec3D1) < 1e-12, "direction of wedge stuff vector");
    logger.assert_always(std::abs((MiddleSizeMatrix({vec3D1, vec3D2}).computeWedgeStuffVector()) * vec3D2) < 1e-12, "direction of wedge stuff vector");
    logger.assert_always(std::abs((MiddleSizeMatrix({vec3D0, vec3D2}).computeWedgeStuffVector()) * vec3D0) < 1e-12, "direction of wedge stuff vector");
    logger.assert_always(std::abs((MiddleSizeMatrix({vec3D0, vec3D2}).computeWedgeStuffVector()) * vec3D2) < 1e-12, "direction of wedge stuff vector");
    MiddleSizeVector vec4D0{0., 1., 2., 3.}, vec4D1{4., 5., 6., 7.}, vec4D2{0., -1., 2., -3.}, vec4D3{0., -1., -2., 3.};
    ///\todo test that the norm of the 4D wedge stuff vector equals the area of the tetrahedron formed by nodes {0, 0, 0}, v1 and v2 and v3
    logger.assert_always(std::abs((MiddleSizeMatrix({vec4D0, vec4D1, vec4D2}).computeWedgeStuffVector()) * vec4D0) < 1e-12, "direction of wedge stuff vector");
    logger.assert_always(std::abs((MiddleSizeMatrix({vec4D0, vec4D1, vec4D2}).computeWedgeStuffVector()) * vec4D1) < 1e-12, "direction of wedge stuff vector");
    logger.assert_always(std::abs((MiddleSizeMatrix({vec4D0, vec4D1, vec4D2}).computeWedgeStuffVector()) * vec4D2) < 1e-12, "direction of wedge stuff vector");
    logger.assert_always(std::abs((MiddleSizeMatrix({vec4D0, vec4D1, vec4D3}).computeWedgeStuffVector()) * vec4D0) < 1e-12, "direction of wedge stuff vector");
    logger.assert_always(std::abs((MiddleSizeMatrix({vec4D0, vec4D1, vec4D3}).computeWedgeStuffVector()) * vec4D1) < 1e-12, "direction of wedge stuff vector");
    logger.assert_always(std::abs((MiddleSizeMatrix({vec4D0, vec4D1, vec4D3}).computeWedgeStuffVector()) * vec4D3) < 1e-12, "direction of wedge stuff vector");
    logger.assert_always(std::abs((MiddleSizeMatrix({vec4D0, vec4D2, vec4D3}).computeWedgeStuffVector()) * vec4D0) < 1e-12, "direction of wedge stuff vector");
    logger.assert_always(std::abs((MiddleSizeMatrix({vec4D0, vec4D2, vec4D3}).computeWedgeStuffVector()) * vec4D2) < 1e-12, "direction of wedge stuff vector");
    logger.assert_always(std::abs((MiddleSizeMatrix({vec4D0, vec4D2, vec4D3}).computeWedgeStuffVector()) * vec4D3) < 1e-12, "direction of wedge stuff vector");
    logger.assert_always(std::abs((MiddleSizeMatrix({vec4D1, vec4D2, vec4D3}).computeWedgeStuffVector()) * vec4D1) < 1e-12, "direction of wedge stuff vector");
    logger.assert_always(std::abs((MiddleSizeMatrix({vec4D1, vec4D2, vec4D3}).computeWedgeStuffVector()) * vec4D2) < 1e-12, "direction of wedge stuff vector");
    logger.assert_always(std::abs((MiddleSizeMatrix({vec4D1, vec4D2, vec4D3}).computeWedgeStuffVector()) * vec4D3) < 1e-12, "direction of wedge stuff vector");

    copy.concatenate(extra);
    logger.assert_always(copy.getNRows() == 4, "Rows in a matrix");
    logger.assert_always(copy.getNCols() == 2, "Columns in a matrix");
    logger.assert_always(copy.size() == 8, "Size of a matrix");
    logger.assert_always(std::abs((copy)(0, 0) - 4.) < 1e-12, "concatenate");
    logger.assert_always(std::abs((copy)(1, 0) - 6.) < 1e-12, "concatenate");
    logger.assert_always(std::abs((copy)(0, 1) - 8.) < 1e-12, "concatenate");
    logger.assert_always(std::abs((copy)(1, 1) - 10.) < 1e-12, "concatenate");
    logger.assert_always(std::abs((copy)(2, 0) + 8.) < 1e-12, "concatenate");
    logger.assert_always(std::abs((copy)(3, 0) + 19.) < 1e-12, "concatenate");
    logger.assert_always(std::abs((copy)(2, 1) + 2.) < 1e-12, "concatenate");
    logger.assert_always(std::abs((copy)(3, 1) + 21.) < 1e-12, "concatenate");

    logger.assert_always(MiddleSizeMatrix({vec3D0, vec3D1}).getColumn(0).size() == 3, "getColumn");
    logger.assert_always(std::abs(((MiddleSizeMatrix({vec3D0, vec3D1}).getColumn(0)) - vec3D0)[0]) < 1e-12, "getColumn");
    logger.assert_always(std::abs(((MiddleSizeMatrix({vec3D0, vec3D1}).getColumn(0)) - vec3D0)[1]) < 1e-12, "getColumn");
    logger.assert_always(std::abs(((MiddleSizeMatrix({vec3D0, vec3D1}).getColumn(0)) - vec3D0)[2]) < 1e-12, "getColumn");
    logger.assert_always(MiddleSizeMatrix({vec3D0, vec3D1}).getColumn(1).size() == 3, "getColumn");
    logger.assert_always(std::abs(((MiddleSizeMatrix({vec3D0, vec3D1}).getColumn(1)) - vec3D1)[0]) < 1e-12, "getColumn");
    logger.assert_always(std::abs(((MiddleSizeMatrix({vec3D0, vec3D1}).getColumn(1)) - vec3D1)[1]) < 1e-12, "getColumn");
    logger.assert_always(std::abs(((MiddleSizeMatrix({vec3D0, vec3D1}).getColumn(1)) - vec3D1)[2]) < 1e-12, "getColumn");
    logger.assert_always(MiddleSizeMatrix({vec3D0, vec3D1}).getRow(0).size() == 2, "getRow");
    logger.assert_always(std::abs(((MiddleSizeMatrix({vec3D0, vec3D1}).getRow(0)))[0] - 0.) < 1e-12, "getRow");
    logger.assert_always(std::abs(((MiddleSizeMatrix({vec3D0, vec3D1}).getRow(0)))[1] - 3.) < 1e-12, "getRow");
    logger.assert_always(MiddleSizeMatrix({vec3D0, vec3D1}).getRow(1).size() == 2, "getRow");
    logger.assert_always(std::abs(((MiddleSizeMatrix({vec3D0, vec3D1}).getRow(1)))[0] - 1.) < 1e-12, "getRow");
    logger.assert_always(std::abs(((MiddleSizeMatrix({vec3D0, vec3D1}).getRow(1)))[1] - 4.) < 1e-12, "getRow");
    logger.assert_always(MiddleSizeMatrix({vec3D0, vec3D1}).getRow(2).size() == 2, "getRow");
    logger.assert_always(std::abs(((MiddleSizeMatrix({vec3D0, vec3D1}).getRow(2)))[0] - 2.) < 1e-12, "getRow");
    logger.assert_always(std::abs(((MiddleSizeMatrix({vec3D0, vec3D1}).getRow(2)))[1] - 5.) < 1e-12, "getRow");

    ///\todo figure out a way to test a LU factorisation

    MiddleSizeVector duplicate = vec2;
    count0.solve(vec2);
    logger.assert_always(std::abs((vec2 - count0.inverse() * duplicate)[0]) < 1e-12, "inverse and solve");
    logger.assert_always(std::abs((vec2 - count0.inverse() * duplicate)[1]) < 1e-12, "inverse and solve");
    logger.assert_always(std::abs((count0 - count0.inverse().inverse())[0]) < 1e-12, "inverse");
    logger.assert_always(std::abs((count0 - count0.inverse().inverse())[1]) < 1e-12, "inverse");
    logger.assert_always(std::abs((count0 - count0.inverse().inverse())[2]) < 1e-12, "inverse");
    logger.assert_always(std::abs((count0 - count0.inverse().inverse())[3]) < 1e-12, "inverse");
    logger.assert_always(std::abs((count0 - count0.transpose().transpose())[0]) < 1e-12, "transpose");
    logger.assert_always(std::abs((count0 - count0.transpose().transpose())[1]) < 1e-12, "transpose");
    logger.assert_always(std::abs((count0 - count0.transpose().transpose())[2]) < 1e-12, "transpose");
    logger.assert_always(std::abs((count0 - count0.transpose().transpose())[3]) < 1e-12, "transpose");
    logger.assert_always(std::abs(A23(0, 0) - A23.transpose()(0, 0)) < 1e-12, "transpose");
    logger.assert_always(std::abs(A23(1, 0) - A23.transpose()(0, 1)) < 1e-12, "transpose");
    logger.assert_always(std::abs(A23(0, 1) - A23.transpose()(1, 0)) < 1e-12, "transpose");
    logger.assert_always(std::abs(A23(1, 1) - A23.transpose()(1, 1)) < 1e-12, "transpose");
    logger.assert_always(std::abs(A23(0, 2) - A23.transpose()(2, 0)) < 1e-12, "transpose");
    logger.assert_always(std::abs(A23(1, 2) - A23.transpose()(2, 1)) < 1e-12, "transpose");

    auto data = A23.data();
    logger.assert_always(std::abs(data[0] - A23[0]) < 1e-12, "data");
    logger.assert_always(std::abs(data[1] - A23[1]) < 1e-12, "data");
    logger.assert_always(std::abs(data[2] - A23[2]) < 1e-12, "data");
    logger.assert_always(std::abs(data[3] - A23[3]) < 1e-12, "data");
    logger.assert_always(std::abs(data[4] - A23[4]) < 1e-12, "data");
    logger.assert_always(std::abs(data[5] - A23[5]) < 1e-12, "data");
    data[2] = 17.3;
    logger.assert_always(std::abs(A23[2] - 17.3) < 1e-12, "data");

    std::cout << A32 << std::endl;


}



//there is only comment below this line
//---------------------------------------------------------------------------------------


/*int main(int argc, char* argv[])
{
    
    //First test the empty constructor
    LinearAlgebra::Matrix AA1;
    
    //Second test the constuctor which defines the size
    LinearAlgebra::Matrix AA2(2, 3);
    
    //Third test the constuctor to initiate to a set value e.g. 2 by 3 matrix with all values initised to 1.5
    LinearAlgebra::Matrix AA3(2, 3, 1.5);
    
    //Check the copy constructor
    LinearAlgebra::Matrix AA4(AA3);
    
    //Check resizing the matrix
    AA4.resize(2, 5);
    
    //Check the Matrix assigment
    AA2 = AA3;
    
    //Now test Matrix times vector
    LinearAlgebra::Matrix BB1(2, 2);
    BB1(0, 0) = 1.0;
    BB1(0, 1) = 2.0;
    BB1(1, 0) = 3.0;
    BB1(1, 1) = 4.0;
    
    logger.assert_always(2.0 == BB1[2], "Test the [] operator, Matrix should column-major.");
    logger.assert_always(3.0 == BB1[1], "Test the [] operator, Matrix should column-major.");
    logger.assert_always(BB1(0,1) == 2.0, "Test the () operator. Expected 2.0, but got %", BB1(0,1));
    
    cout << "This is BB1 \n" << BB1 << "\n";
    
    //now test divide
    
    BB1 /= 2.0;
    cout << "This is BB1 divided by 2 \n" << BB1 << "\n";
    
    cout << "and by 2.0 again but using inline divide \n" << BB1 / 2.0 << "\n";
    
    //Now test output
    
    cout << BB1 << "\n";
    
    cout << AA3 << endl;
    
    LinearAlgebra::NumericalVector B2(2);
    LinearAlgebra::NumericalVector B3;
    B2(0) = 1.0;
    B2(1) = 2.0;
    
    cout << "axpy" << endl;
    B3 = BB1 * B2;
    
    cout << "axpy" << endl;
    //Test wedge product both quick form if vector exist and create vector form.
    //Not supposed to be used for square matrices, turned of for now
    //B2=BB1.computeWedgeStuffVector();
    //B3=BB1.computeWedgeStuffVector();
    
    //Matrix assigment test
    AA4 = 2.0;
    
    //Matrix times matrix test
    LinearAlgebra::Matrix CC1(2, 3, 2);
    LinearAlgebra::Matrix CC2(3, 2, 2);
    LinearAlgebra::Matrix CC3(2, 2);
    
    CC3 = CC1 * CC2;
    
    cout << CC1 << std::endl;
    
    LinearAlgebra::Matrix CC2_fix(3, 2, 4);
    std::cout << "CC2_fix PRE" << std::endl;
    std::cout << CC2_fix << std::endl;
    
    CC2.axpy(2.0, CC2_fix);
    
    std::cout << "CC2_fix POST" << std::endl;
    std::cout << CC2_fix << std::endl;
    
    std::cout << CC2 << std::endl;
    
    CC3 = BB1;
    CC3 = BB1.LUfactorisation();
    
    cout << "\n Now the LU factorisation \n";
    
    cout << "Before : \n";
    
    cout << BB1;
    
    cout << "After : \n";
    cout << CC3;
    
    cout << "\n Now the inverse" << endl;
    
    CC3(0, 0) = 0.8147;
    CC3(0, 1) = 0.1270;
    CC3(1, 0) = 0.9058;
    CC3(1, 1) = 0.9134;
    
    LinearAlgebra::Matrix CC4(2, 2);
    
    CC4 = CC3.inverse();
    
    cout << CC4;
    
    CC3 = CC3 * CC4;
    
    cout << CC3;
    
    cout << "\n Now test the solution of Ax=B \n ";
    
    CC3.solve(CC4);
    
    cout << CC4;
    
    cout << "\n Now test the inverse of a 3 by3 \n";
    
    LinearAlgebra::Matrix DD(3, 3);
    
    DD(0, 0) = 0.6324;
    DD(0, 1) = 0.5469;
    DD(0, 2) = 0.1576;
    DD(1, 0) = 0.0975;
    DD(1, 1) = 0.9575;
    DD(1, 2) = 0.9706;
    DD(2, 0) = 0.2785;
    DD(2, 1) = 0.9649;
    DD(2, 2) = 0.9572;
    
    LinearAlgebra::Matrix ans(3, 3);
    
    ans = DD.inverse();
    
    cout << ans << std::endl;
    
    return 0;
    
}

//	LinearAlgebra::Matrix<double> AA(3,2);
//	
//	LinearAlgebra::Matrix<double> BB(2,3);
//	
//	LinearAlgebra::Matrix<double> CC(2,4);
//	
//	LinearAlgebra::Matrix<double> DD(3,2);
//	
//	LinearAlgebra::Matrix<double> ANS;
//    
//    LinearAlgebra::NumericalVector x(2);
//    LinearAlgebra::NumericalVector y(2);
//    
//    LinearAlgebra::Matrix<double> A(2,2);
//    
//    
//	x(0)=1.0;
//    x(1)=2.0;
//    
//    A(0,0)=1.0;
//    A(0,1)=2.0;
//    A(1,0)=3.0;
//    A(1,1)=4.0;
//	
//	AA(0,0)=1;
//	AA(0,1)=2;
//	AA(1,0)=3;
//	AA(1,1)=4;
//	AA(2,0)=5;
//	AA(2,1)=6;
//    
//    cout << x;
//    
//    y=AA*x;
//    
//    cout << "Help" << endl;
//    cout << y;
//    cout << "and now " <<endl;
//	
//	BB(0,0)=7;
//	BB(0,1)=8;
//	BB(0,2)=9;
//	BB(1,0)=10;
//	BB(1,1)=11;
//	BB(1,2)=12;
//	
//	CC(0,0)=7;
//	CC(0,1)=8;
//	CC(0,2)=9;
//	CC(0,3)=10;
//	CC(1,0)=11;
//	CC(1,1)=12;
//	CC(1,2)=13;
//	CC(1,3)=14;
//	
//	DD(0,0)=6;
//	DD(0,1)=5;
//	DD(1,0)=4;
//	DD(1,1)=3;
//	DD(2,0)=2;
//	DD(2,1)=1;
//	
//	//AA*=BB;
//	cout << AA;
//	//AA*=CC;
//	cout << AA;
//	
//	ANS=AA+DD;
//	cout << ANS;
//	//cout << DD;
//	//cout << AA;
//	
//	
// 	;

*/
