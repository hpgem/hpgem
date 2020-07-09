/*
This file forms part of hpGEM. This package has been developed over a number of
years by various people at the University of Twente and a full list of
contributors can be found at http://hpgem.org/about-the-code/team

This code is distributed using BSD 3-Clause License. A copy of which can found
below.


Copyright (c) 2020, Univesity of Twenete
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#ifndef HPGEM_APP_ABSTRACTEIGENVALUERESULT_H
#define HPGEM_APP_ABSTRACTEIGENVALUERESULT_H

#include <fstream>
#include <iostream>

#include "EigenvalueProblem.h"

using namespace hpgem;


/// Result of solving the EigenvalueProblem
template <std::size_t DIM>
class AbstractEigenvalueResult {
   public:
    virtual ~AbstractEigenvalueResult() = default;

    /// The problem that was solved
    virtual const EigenvalueProblem<DIM>& originalProblem() const = 0;
    /// A list of frequencies (in increasing order) obtained at a certain
    /// k-point
    virtual const std::vector<double> frequencies(std::size_t point) const = 0;

    void writeFrequencies(const std::string& fileName) const {
        std::ofstream file;
        file.open(fileName);
        writeFrequencies(file, ',');
        file.close();
    }

    /// Write the frequencies to the standard out.
    void printFrequencies() const { writeFrequencies(std::cout, '\t'); }

    /// Write the frequencies to a stream
    ///
    /// \param stream The output stream to write to
    /// \param separator The separator between frequencies for the same k-point
    void writeFrequencies(std::ostream& stream, char separator) const {
        for (std::size_t i = 0;
             i < originalProblem().getPath().totalNumberOfSteps(); ++i) {
            std::vector<double> freqs = frequencies(i);
            stream << i;
            for (double& freq : freqs) {
                stream << separator << freq;
            }
            stream << std::endl;
        }
    }
};

#endif  // HPGEM_APP_ABSTRACTEIGENVALUERESULT_H
