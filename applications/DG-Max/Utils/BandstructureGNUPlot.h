
#ifndef HPGEM_APP_BANDSTRUCTUREGNUPLOT_H
#define HPGEM_APP_BANDSTRUCTUREGNUPLOT_H

#include "KSpacePath.h"
#include "BandStructure.h"
#include "ProblemTypes/AbstractEigenvalueResult.h"

#include <string>
#include <vector>

/// Utility class to generate plots of a band structure using gnuplot.
/// \tparam DIM The dimension to work in.
template <std::size_t DIM>
class BandstructureGNUPlot {
   public:
    BandstructureGNUPlot(
        const KSpacePath<DIM>& path, const std::vector<std::string>& pointNames,
        const BandStructure<DIM>& structure,
        const AbstractEigenvalueResult<DIM>* computedSpectrum = nullptr);
    void plot(std::string fileName);

   private:
    struct Line {
        Line(std::string styling, std::string title, std::string data)
            : styling_(styling),
              data_(data),
              title_(title),
              deduplicate_(false),
              priority_(0) {}

        std::string styling_;
        std::string data_;
        std::string title_;
        bool deduplicate_;
        int priority_;

        bool operator==(const Line& other) { return this == &other; }
    };

    Line band(const typename BandStructure<DIM>::LineSet& line,
              std::size_t lineIndex, double x1, double x2, bool titled);

    std::map<int, std::vector<std::tuple<double, double>>> groupSpectrum();

    const BandStructure<DIM>& structure_;
    const AbstractEigenvalueResult<DIM>* computedSpectrum_;
    const KSpacePath<DIM>& path_;
    const std::vector<std::string>& pointNames_;
};

#endif // HPGEM_APP_BANDSTRUCTUREGNUPLOT_H
