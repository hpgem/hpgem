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
#include "Geometry/PointPhysical.h"
#include "Base/Element.h"
#include "Geometry/PhysicalGeometry.h"

#include <vector>
#include "CMakeDefinitions.h"
#include <sstream>
#include "Output/TecplotDiscontinuousSolutionWriter.h"
#include "Base/ElementCacheData.h"
#include "Base/FaceCacheData.h"

using namespace std;

class Dummy
{
public:
    Dummy()
    {
    }
    void operator()(const Base::Element& el, const Geometry::PointReference& p, ostream& os)
    {
    }
};

//Also test the copy constructor of Mesh/MeshManipulator here after this test has been fixed.
int main()
{
    /*    Base::ConfigurationData config(3,1,1);
     
     config.numberOfUnknowns_       = 1;
     config.numberOfTimeLevels_     = 1;
     config.numberOfBasisFunctions_ = 1;
     
     Base::MeshManipulator<2> myTwoDDemoMesh(&config, 1,1);
     
     std::stringstream filename;
     
     filename << getCMAKE_hpGEM_SOURCE_DIR() << "/tests/files/centaurQuadMinimum.hyb";

     myTwoDDemoMesh.readCentaurMesh(filename.str());
     
     std::cout << myTwoDDemoMesh << std::endl;
     
     
     std::ofstream file2D;
     file2D.open ("SavedCentaurQuadMinimum.dat");
     
     
     Output::TecplotDiscontinuousSolutionWriter<2> out(file2D,"QuadMinimum Test Mesh","01","xy");
     Dummy d;
     out.write(&myTwoDDemoMesh,"holi",false, d);
     
     file2D.close();
     
     //Now do it again with a more complicated mesh
     
     Base::MeshManipulator<2> triQuadTwoDDemoMesh(&config,1,1);
     
     std::stringstream filename2;
     
     filename2 << getCMAKE_hpGEM_SOURCE_DIR() << "/tests/files/centaurMinimum.hyb";
     
     cout << "filename " << filename2.str() <<endl;
     
     triQuadTwoDDemoMesh.readCentaurMesh(filename2.str());
     
     std::cout << triQuadTwoDDemoMesh << std::endl;
     
     file2D.open ("SavedCentaurMinimum.dat");
     
     Output::TecplotDiscontinuousSolutionWriter<2> out2(file2D,"Minimum Test Mesh","01","x,y");
     out2.write(&triQuadTwoDDemoMesh,"holi",false,d);
     
     file2D.close();
     
     */
}
