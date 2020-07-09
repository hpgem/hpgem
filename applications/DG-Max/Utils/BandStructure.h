
#ifndef HPGEM_APP_BANDSTRUCTURE_H
#define HPGEM_APP_BANDSTRUCTURE_H

#include "LinearAlgebra/SmallVector.h"
#include <map>
#include <memory>
#include <sstream>

template <std::size_t DIM>
class BandStructure {
   public:
    class LineSet {
       public:
        virtual ~LineSet() = default;
        virtual std::size_t multiplicity(std::size_t line) const = 0;
        virtual double frequency(std::size_t line,
                                 double interpolation) const = 0;
        virtual std::size_t numberOfLines() const = 0;

        /// \brief Numeric type of the line, used to determine the line type in
        /// plotting. Use -1 if unused.
        virtual int lineType(std::size_t line) const { return -1; }

        /// \brief Title of the line used for plotting
        virtual std::string lineTitle(std::size_t line) const {
            std::stringstream out;
            out << multiplicity(line);
            return out.str();
        }
    };

    virtual ~BandStructure() = default;
    virtual std::map<double, std::size_t> computeSpectrum(
        LinearAlgebra::SmallVector<DIM> kpoint, double omegaMax) const = 0;
    virtual std::unique_ptr<LineSet> computeLines(
        LinearAlgebra::SmallVector<DIM> point1,
        LinearAlgebra::SmallVector<DIM> point2, double maxFrequency) const = 0;
};

#endif  // HPGEM_APP_BANDSTRUCTURE_H
