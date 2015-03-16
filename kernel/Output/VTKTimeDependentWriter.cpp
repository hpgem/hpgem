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

#include "VTKTimeDependentWriter.h"
#include "Base/MpiContainer.h"
#include "Logger.h"
#include "base64.h"
#include "Base/CommandLineOptions.h"

Output::VTKTimeDependentWriter::VTKTimeDependentWriter(std::string baseFileName, Base::MeshManipulator* mesh)
        : baseName_(baseFileName), mesh_(mesh), currentFile_(nullptr), time_(0), numberOfFilesWritten_(0)
{
    logger.assert(mesh!=nullptr,"Invalid mesh passed");
    std::size_t id = Base::MPIContainer::Instance().getProcessorID();
    if (id == 0)
    {
        masterFile_.open(baseFileName + ".pvd");
        if (!masterFile_.good())
        {
            if (baseFileName.find('/') != std::string::npos)
            {
                logger(FATAL, "failed to open main paraview output file %.pvd, does the directory % exist?", baseFileName, baseFileName.substr(0, baseFileName.find_last_of('/') + 1));
            }
            else
            {
                logger(FATAL, "failed to open main paraview output file %.pvd", baseFileName);
            }
            exit(1);
        }
        masterFile_ << "<?xml version=\"1.0\"?>" << std::endl;
        masterFile_ << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"" << (Detail::isBigEndian() ? "BigEndian" : "LittleEndian") << "\">" << std::endl;
        masterFile_ << "  <Collection>" << std::endl;
    }
}

Output::VTKTimeDependentWriter::~VTKTimeDependentWriter()
{
    masterFile_ << "  </Collection>" << std::endl;
    masterFile_ << "</VTKFile>" << std::endl;
    masterFile_.flush();
    masterFile_.close();
    if (currentFile_ != nullptr)
    {
        delete currentFile_;
    }
    else
    {
        logger(ERROR, "no time levels written!");
    }
}

//3x the same function, but I dont like this mess in the header, so cant template
void Output::VTKTimeDependentWriter::write(std::function<double(Base::Element*, const Geometry::PointReference&, std::size_t)> f, const std::string& name, double time, std::size_t timelevel)
{
    if (time != time_ || currentFile_ == nullptr)
    {
        if (currentFile_ != nullptr)
        {
            delete currentFile_;
        }
        //convert the current time to some unique number as part of the base name for the file for the new timelevel
        std::string fileName = baseName_ + std::to_string(numberOfFilesWritten_);
        numberOfFilesWritten_++;
        currentFile_ = new VTKSpecificTimeWriter {fileName, mesh_, timelevel};
        if (fileName.find('/') != std::string::npos)
        {
            fileName = fileName.substr(fileName.find_last_of('/') + 1);
        }
        masterFile_ << "    <DataSet timestep=\"" << time << "\" group=\"\" part=\"0\" file=\"" << fileName << ".pvtu\"/>" << std::endl;
        time_ = time;
        timelevel_ = timelevel;
    }
    logger.assert(timelevel == timelevel_, "Timelevel isn't as expected. % != %", timelevel, timelevel_);
    currentFile_->write(f, name);
}

//3x the same function, but I dont like this mess in the header, so cant template
void Output::VTKTimeDependentWriter::write(std::function<LinearAlgebra::NumericalVector(Base::Element*, const Geometry::PointReference&, std::size_t)> f, const std::string& name, double time, std::size_t timelevel)
{
    if (time != time_ || currentFile_ == nullptr)
    {
        if (currentFile_ != nullptr)
        {
            delete currentFile_;
        }
        //convert the current time to some unique number as part of the base name for the file for the new timelevel
        std::string fileName = baseName_ + Detail::toBase64(&time, sizeof(double));
        currentFile_ = new VTKSpecificTimeWriter {fileName, mesh_, timelevel};
        masterFile_ << "<DataSet timestep=\"" << time - time_ << "\" group=\"\" part=\"0\" file=\"" << fileName << ".pvtu\"/>" << std::endl;
        time_ = time;
        timelevel_ = timelevel;
    }
    logger.assert(timelevel == timelevel_, "Timelevel isn't as expected. % != %", timelevel, timelevel_);
    currentFile_->write(f, name);
}

//3x the same function, but I dont like this mess in the header, so cant template
void Output::VTKTimeDependentWriter::write(std::function<LinearAlgebra::Matrix(Base::Element*, const Geometry::PointReference&, std::size_t)> f, const std::string& name, double time, std::size_t timelevel)
{
    if (time != time_ || currentFile_ == nullptr)
    {
        if (currentFile_ != nullptr)
        {
            delete currentFile_;
        }
        //convert the current time to some unique number as part of the base name for the file for the new timelevel
        std::string fileName = baseName_ + Detail::toBase64(&time, sizeof(double));
        currentFile_ = new VTKSpecificTimeWriter {fileName, mesh_, timelevel};
        masterFile_ << "<DataSet timestep=\"" << time - time_ << "\" group=\"\" part=\"0\" file=\"" << fileName << ".pvtu\"/>" << std::endl;
        time_ = time;
        timelevel_ = timelevel;
    }
    logger.assert(timelevel == timelevel_, "Timelevel isn't as expected. % != %", timelevel, timelevel_);
    currentFile_->write(f, name);
}

