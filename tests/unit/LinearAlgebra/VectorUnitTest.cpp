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
#include <stdlib.h>
#include <math.h>
#include "LinearAlgebra/NumericalVector.h"

using std::cout;
using std::endl;

int main(int argc, char* argv[])
{
    
    LinearAlgebra::NumericalVector a(3);
    LinearAlgebra::NumericalVector b(3);
    
    cout << a.size();
    
    a[0] = 2.0;
    a[1] = 3.0;
    a[2] = 1.0;
    
    cout << a << endl;
    
    b = 3.0 * a;
    
    cout << b << endl;
    
    a.axpy(3.0, b);
    cout << a << endl;
    
    LinearAlgebra::NumericalVector c;
    
    LinearAlgebra::NumericalVector AA(3);
    
    LinearAlgebra::NumericalVector BB(3);
    
    LinearAlgebra::NumericalVector CC(3);
    
    LinearAlgebra::NumericalVector DD(3);
    
    AA(0) = 1;
    AA(1) = 2;
    AA(2) = 3;
    
    BB(0) = 7;
    BB(1) = 8;
    BB(2) = 9;
    
    cout << AA << endl;
    
    DD(0) = 1;
    DD(1) = 2;
    DD(2) = 3;
    
    DD = AA;
    
    AA += BB;
    
    cout << AA << endl;
    
    CC = AA + BB;
    
    cout << CC << endl;
    
}

