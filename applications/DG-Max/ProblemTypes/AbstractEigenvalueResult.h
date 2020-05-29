
#ifndef HPGEM_BASEEIGENVALUEOUTPUT_H
#define HPGEM_BASEEIGENVALUEOUTPUT_H

#include "EigenValueProblem.h"

#include <iostream>
#include <fstream>

template <std::size_t DIM>
class AbstractEigenvalueResult {
   public:
    virtual ~AbstractEigenvalueResult() = default;

    virtual const EigenValueProblem<DIM>& originalProblem() const = 0;
    virtual const std::vector<double> frequencies(std::size_t point) const = 0;

    void writeFrequencies(std::string fileName) const {
        std::ofstream file;
        file.open(fileName);
        writeFrequencies(file, ',');
        file.close();
    }

    void printFrequencies() const { writeFrequencies(std::cout, '\t'); }

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

#endif  // HPGEM_BASEEIGENVALUEOUTPUT_H
