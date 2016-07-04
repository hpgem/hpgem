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

#include "CommandLineOptions.h"
#include <limits>

namespace Base
{

    ///\bug Workaround for Bug 60352 in (at least) gcc 4.8.2 (should read auto& numberOfSnapshots = ...)
    CommandLineOption<std::size_t>& numberOfSnapshots = Base::register_argument<std::size_t>(0, "nOutputFrames", "Number of frames to output", false, 1);
    CommandLineOption<double>& endTime = Base::register_argument<double>(0, "endTime", "end time of the simulation", false, 1);
    CommandLineOption<double>& startTime = Base::register_argument<double>(0, "startTime", "start time of the simulation", false, 0);
    CommandLineOption<double>& dt = Base::register_argument<double>(0, "dt", "maximum time step of the simulation", false, 0.01);
    CommandLineOption<double>& error = Base::register_argument<double>(0, "error", "maximum acceptable relative error per time step", false, std::numeric_limits<double>::infinity());
    CommandLineOption<std::string>& outputName = Base::register_argument<std::string>(0, "outFile", "Name of the output file (without extentions)", false, "output");
}
