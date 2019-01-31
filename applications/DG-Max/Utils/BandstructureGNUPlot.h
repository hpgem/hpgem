
#ifndef HPGEM_GNUOUTPUT_H
#define HPGEM_GNUOUTPUT_H

#include "KSpacePath.h"
#include "HomogeneousBandStructure.h"
#include "ProblemTypes/BaseEigenvalueResult.h"

#include <string>
#include <vector>

/// Utility class to generate plots of a band structure using gnuplot.
/// \tparam DIM The dimension to work in.
template<std::size_t DIM>
class BandstructureGNUPlot
{
public:
    BandstructureGNUPlot(const KSpacePath<DIM>& path,
            const std::vector<std::string>& pointNames,
            const HomogeneousBandStructure<DIM>& structure,
            const BaseEigenvalueResult<DIM>* computedSpectrum = nullptr);
    void plot(std::string fileName);
private:
    std::string band(const typename HomogeneousBandStructure<DIM>::Line& line,
        double x1, double x2, bool titled);

    std::map<int, std::vector<std::tuple<double, double>>> groupSpectrum();

    const HomogeneousBandStructure<DIM>& structure_;
    const BaseEigenvalueResult<DIM>* computedSpectrum_;
    const KSpacePath<DIM>& path_;
    const std::vector<std::string>& pointNames_;
};


#endif //HPGEM_GNUOUTPUT_H
