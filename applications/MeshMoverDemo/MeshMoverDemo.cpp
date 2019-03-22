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

#include "Base/HpgemAPIBase.h"
#include "MeshMover.h"
#include "Base/GlobalData.h"
#include "Base/ConfigurationData.h"
#include "Output/TecplotSingleElementWriter.h"
#include "Output/TecplotDiscontinuousSolutionWriter.h"
#include "Output/TecplotPhysicalGeometryIterator.h"
#include "Output/VTKSpecificTimeWriter.h"
using Base::ConfigurationData;
using Base::GlobalData;

//Note: the intended use of the prototype classes is to merge Dummy with MeshMoverExampleProblem
class Dummy : public Output::TecplotSingleElementWriter<2>
{
public:
    Dummy()
    {
    }
    void writeToTecplotFile(const Base::Element* el, const Geometry::PointReference<2>& p, std::ostream& os)
    {
    }
};

///Demo of how a mesh mover can be used to create a non-rectangular domain. 
class MeshMoverExampleProblem : public Base::HpgemAPIBase<2>
{
    
public:
    MeshMoverExampleProblem(GlobalData* const global, const ConfigurationData* config)
            : Base::HpgemAPIBase<2>(global, config)
    {
    }
    
    ///\brief Make the initial mesh.
    /// Make the initial rectangular mesh and add it to the application. Also 
    /// initialise the mesh mover of MeshMover.h and couple it to the mesh.
    bool initialise()
    {
        
        std::size_t id = addMesh("mesh.hpgem");
        
        //Set up the move of the mesh; note that the mesh mover gets deleted in the mesh manipulator
        const MeshMover* meshMover = new MeshMover;
        initialiseMeshMover(meshMover, id);
        
        return true;
    }
    
    ///\brief Write the mesh to output files. 
    ///The file with extension .dat is meant to be used for Tecplot, the files 
    ///with extensions .pvtu and .vtu can be read in for example Paraview or 
    ///other applications that can read VTK files.
    ///Note that we are only interested in the mesh, so the VTK file just plots
    ///the value 0 in each point.
    void output()
    {
        std::ofstream file2D;
        file2D.open("movedMesh.dat");
        Output::TecplotDiscontinuousSolutionWriter<2> out(file2D, "RectangularMesh", "01", "xy");
        
        Dummy d;
        out.write(meshes_[0], "holi", false, &d);
        
        Output::VTKSpecificTimeWriter<2> vtkOut("movedMesh", meshes_[0]);
        std::function<double(Base::Element*, const Geometry::PointReference<2>&, std::size_t)> fun = 
        [](Base::Element* elt, const Geometry::PointReference<2>& pRef, std::size_t i){return 0.;};
        vtkOut.write(fun, "zero");
    }
    
    ///\brief Move the mesh as described by the mesh mover.
    void solve()
    {
        
        meshes_[0]->move();
        
    }
    
};

int main(int argc, char **argv)
{
    Base::parse_options(argc, argv);
    
    GlobalData* globalData = new GlobalData();
    
    ConfigurationData* config = new ConfigurationData(1, 1);
    
    config->numberOfUnknowns_ = 1;
    config->numberOfTimeLevels_ = 1;
    config->numberOfBasisFunctions_ = 1;
    
    globalData->numberOfUnknowns_ = 10;
    globalData->numberOfTimeLevels_ = 1;
    
    MeshMoverExampleProblem problem(globalData, config);
    
    problem.initialise();
    
    problem.solve();
    
    problem.output();
}
