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

#ifndef TECPLOTDISCONTINUOUSSOLUTIONWRITER_HH
#define TECPLOTDISCONTINUOUSSOLUTIONWRITER_HH

#include <sstream>
#include <ostream>
#include <vector>
#include <utility>
#include <string>
#include <functional>

#include "Logger.h"

namespace Base
{
    template<std::size_t DIM>
    class MeshManipulator;
    class Element;
}

namespace Geometry
{
    template<std::size_t DIM>
    class PointReference;
}

namespace Output
{
    template<std::size_t DIM>
    class TecplotSingleElementWriter;
    
    //! \brief This class prints the nodes and the solution in every element in Tecplot format.
    template<std::size_t DIM>
    class TecplotDiscontinuousSolutionWriter
    {
    public:

        TecplotDiscontinuousSolutionWriter(std::ostream& output, const std::string& fileTitle, const std::string& dimensionsToWrite, const std::string& variableString);
        
        ///No copy constructor, since we don't want to assign the same ostream to different writers at the same time.
        TecplotDiscontinuousSolutionWriter(const TecplotDiscontinuousSolutionWriter &other) = delete;

        /// Write a zone with data from the current mesh to the stream held by the object.
        ///\deprecated please use the other write instead
        void write(const Base::MeshManipulator<DIM>* mesh, const std::string& zoneTitle, const bool sameGeometry, TecplotSingleElementWriter<DIM>* writeDataClass, const double time = 0);
        
        /// Write a zone with data from the current mesh to the stream held by the object.
        void write(const Base::MeshManipulator<DIM>* mesh, const std::string& zoneTitle, const bool sameGeometry, std::function<void(const Base::Element*, const Geometry::PointReference<DIM>&, std::ostream&)>, const double time = 0);

        ~TecplotDiscontinuousSolutionWriter()
        {
            //assertion is technically checking internal state, but the writability of the filesystem may change outside the influence of this class
            logger.assert_always(output_.good(), "Something is not so good about the output stream");
            output_.flush();
            delete[] dimensionNumbers;
        }
        
    private:
        
        std::string makeTecplotVariableString(const std::string& s) const;

        std::ostream& output_;

        std::size_t previousNumberOfElements_;

        std::size_t previousNumberOfNodes_;

        std::string elementType_[5];

        const std::size_t numberOfDimensionsToWrite_;

        std::size_t* dimensionNumbers;
    };
}

#include "TecplotDiscontinuousSolutionWriter_Impl.h"
#endif
