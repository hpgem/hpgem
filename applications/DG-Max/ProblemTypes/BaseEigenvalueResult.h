
#ifndef HPGEM_BASEEIGENVALUEOUTPUT_H
#define HPGEM_BASEEIGENVALUEOUTPUT_H

#include "EigenValueProblem.h"

template<std::size_t DIM>
class BaseEigenvalueResult
{
public:
    virtual ~BaseEigenvalueResult() = default;

    virtual const EigenValueProblem<DIM>& originalProblem()  const = 0;
    virtual const std::vector<double> frequencies(std::size_t point)  const = 0;

    void printFrequencies()
    {
        for(std::size_t i = 0; i < originalProblem().getPath().totalNumberOfSteps(); ++i)
        {
            std::vector<double> freqs = frequencies(i);
            std::cout << i;
            for(double& freq : freqs)
            {
                std::cout << "\t" << freq;
            }
            std::cout << std::endl;
        }
    }
};

#endif //HPGEM_BASEEIGENVALUEOUTPUT_H
