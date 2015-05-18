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

#ifndef VTKTIMEDEPENDENTWRITER_HPP
#define	VTKTIMEDEPENDENTWRITER_HPP

#include <functional>
#include <fstream>
#include "Base/MeshManipulator.h"
#include "VTKSpecificTimeWriter.h"

namespace Output
{
    
    ///\class VTKSpecificTimeWriter
    ///\brief writes steady state data or a single time level
    ///
    ///this produces multiple files, you do not have to append anything, just load the .pvd into paraview
    ///once you have written at time x, with timelevel n, you must remain using timelevel n for time x
    ///once you write at time y, where y!=x, you can no longer write to time x
    ///also I'm not sure if VTK allows negative time steps, so be careful with those
    ///If you thing this class produces too much output files, use an existing folder as part of your filename
    //class is final because the destructor would be the only virtual function
    class VTKTimeDependentWriter final
    {
    public:
        ///\brief write front matter and open file stream
        ///\param baseFileName name of the file WITHOUT extentions
        ///\param mesh the mesh containing the data you want to output
        ///if you want to write from multiple meshes, simply have paraview load both output files
        VTKTimeDependentWriter(std::string baseFileName, Base::MeshManipulator* mesh);

        ///\brief write end matter and close the file stream
        ~VTKTimeDependentWriter();

        ///\brief write a scalar field
        void write(std::function<double(Base::Element*, const Geometry::PointReference&, std::size_t)>, const std::string& name, double time, std::size_t timelevel = 0);

        ///\brief write a vector field
        void write(std::function<LinearAlgebra::MiddleSizeVector(Base::Element*, const Geometry::PointReference&, std::size_t)>, const std::string& name, double time, std::size_t timelevel = 0);

        ///\brief write an order 2 tensor field
        void write(std::function<LinearAlgebra::MiddleSizeMatrix(Base::Element*, const Geometry::PointReference&, std::size_t)>, const std::string& name, double time, std::size_t timelevel = 0);

        ///\brief do not copy the writer to prevent havoc when destructing all the copies
        VTKTimeDependentWriter(const VTKTimeDependentWriter& orig) = delete;
        VTKTimeDependentWriter operator=(const VTKTimeDependentWriter& orig) = delete;
    private:
        std::string baseName_;
        std::ofstream masterFile_;
        const Base::MeshManipulator* mesh_;
        VTKSpecificTimeWriter* currentFile_;
        std::size_t timelevel_;
        double time_;
        std::size_t numberOfFilesWritten_;
    };

}
#endif	/* VTKTIMEDEPENDENTWRITER_HPP */

