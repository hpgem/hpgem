
#include "BandstructureGNUPlot.h"
#include "Logger.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <set>

template <std::size_t DIM>
BandstructureGNUPlot<DIM>::BandstructureGNUPlot(
    const KSpacePath<DIM>& path, const std::vector<std::string>& pointNames,
    const BandStructure<DIM>& structure,
    const BaseEigenvalueResult<DIM>* computedSpectrum)
    : path_(path),
      pointNames_(pointNames),
      structure_(structure),
      computedSpectrum_(computedSpectrum) {
    logger.assert_always(path_.numberOfCornerPoints() == pointNames_.size() ||
                             pointNames_.empty(),
                         "Mismatch in path length and number of point names.");
    // It might even be better to check if the two paths match, or take the path
    // from the computed spectrum.
    logger.assert_always(
        computedSpectrum_ == nullptr ||
            computedSpectrum_->originalProblem()
                    .getPath()
                    .totalNumberOfSteps() == path.totalNumberOfSteps(),
        "Computed and path step length mismatch.");
}

template <std::size_t DIM>
void BandstructureGNUPlot<DIM>::plot(std::string fileName) {
    const double maxFreq = M_PI;

    std::ofstream out(fileName);
    // Header
    out << "set term postscript eps color enhanced" << std::endl;
    out << "set output 'bands.eps'" << std::endl;

    std::vector<double> xpoints(path_.numberOfCornerPoints());
    double totalDistance = 0;
    xpoints[0] = 0;
    for (std::size_t i = 1; i < path_.numberOfCornerPoints(); ++i) {
        totalDistance += (path_.kcorner(i - 1) - path_.kcorner(i)).l2Norm();
        xpoints[i] = totalDistance;
    }
    // Set x,y range
    out << "xmin=" << 0.0 << ";xmax=" << totalDistance << std::endl;
    out << "ymin=" << 0.0 << ";ymax=" << maxFreq << std::endl;
    out << "set xrange [xmin:xmax]" << std::endl;
    out << "set yrange [ymin:ymax]" << std::endl;
    // Move the legend/key to the bottom, reducing overlap
    out << "set key bottom" << std::endl;

    // Set the x-ticks
    if (!pointNames_.empty()) {
        out << "set xtics (";
        for (std::size_t i = 0; i < path_.numberOfCornerPoints(); ++i) {
            if (i != 0) {
                out << ", ";
            }
            out << "'" << pointNames_[i] << "'" << xpoints[i];
        }
        out << ")" << std::endl;
    }
    // Set y-label
    out << "set ylabel \"Reduced frequency\"" << std::endl;

    std::vector<Line> gnulines;

    for (std::size_t i = 1; i < path_.numberOfCornerPoints(); ++i) {
        std::unique_ptr<typename BandStructure<DIM>::LineSet> lines =
            structure_.computeLines(path_.kcorner(i - 1), path_.kcorner(i),
                                    maxFreq);
        for (std::size_t j = 0; j < lines->numberOfLines(); ++j) {
            Line pline = band(*lines, j, xpoints[i - 1], xpoints[i], "");
            gnulines.emplace_back(pline);
        }
    }

    // Frequency data
    std::map<int, std::vector<std::tuple<double, double>>> points;
    if (computedSpectrum_ != nullptr) {
        points = groupSpectrum();

        for (auto& ps : points) {
            std::ostringstream data;
            std::ostringstream title;
            std::ostringstream style;

            for (auto p : ps.second) {
                data << std::get<0>(p) << " " << std::get<1>(p) << std::endl;
            }
            data << std::endl;

            int difference = ps.first;
            int pointType;
            if (difference < 0)
                pointType = 2;  //'x'
            else if (difference > 0)
                pointType = 1;  // '+'
            else
                pointType = 7;  // 'o'

            style << "points pointtype " << pointType;
            if (difference == 0) {
                // the points are rather big for something that is correct so
                // reduce their size a bit.
                style << " pointsize 0.3";
            }
            style << " linecolor " << std::abs(difference) + 1;
            title << '"' << difference << '"';
            gnulines.emplace_back(style.str(), title.str(), data.str());
            gnulines.back().priority_ = difference;
        }
    }

    out << "$data << EOD" << std::endl;
    for (auto& line : gnulines) {
        // data should be ended with a newline internally.
        out << line.data_ << std::endl << std::endl;
    }
    out << "EOD" << std::endl;

    // Plot the actual structure
    out << "plot ";
    bool first = true;
    int index = 0;

    // Group everything by priority

    std::map<int, std::vector<Line*>> plotOrder;
    for (Line& line : gnulines) {
        std::vector<Line*>& plines = plotOrder[line.priority_];
        plines.emplace_back(&line);
    }

    // For deduplication of the titles, note that we do not consider style.
    std::set<std::string> titles;

    for (auto& plines : plotOrder) {
        for (Line* line : plines.second) {
            bool title = true;
            if (line->deduplicate_) {
                if (titles.find(line->title_) == titles.end()) {
                    titles.emplace(line->title_);
                } else {
                    title = false;
                }
            }

            if (!first) {
                out << ",\\" << std::endl << "  ";
            } else {
                first = false;
            }
            int index = std::find(gnulines.begin(), gnulines.end(), *line) -
                        gnulines.begin();

            out << "$data index " << index << " with " << line->styling_;
            if (title) {
                out << " title " << line->title_;
            } else {
                out << " notitle";
            }
        }
    }
    out << std::endl;
    out.close();
}

template <std::size_t DIM>
typename BandstructureGNUPlot<DIM>::Line BandstructureGNUPlot<DIM>::band(
    const typename BandStructure<DIM>::LineSet& line, std::size_t lineIndex,
    double x1, double x2, bool titled) {
    std::ostringstream data;
    std::ostringstream title;
    std::ostringstream style;

    // TODO: Magic number of points
    const std::size_t POINTS = 81;
    const double dx = (x2 - x1) / (POINTS - 1);
    for (std::size_t i = 0; i < POINTS; ++i) {
        double x = i * dx + x1;
        double f = line.frequency(lineIndex, i / (POINTS - 1.0));
        data << x << " " << f << std::endl;
    }
    std::size_t multiplicity = line.multiplicity(lineIndex);
    title << '"' << line.lineTitle(lineIndex) << '"';
    int dashType = line.lineType(lineIndex);
    if (dashType == -1) {
        dashType = (1 + (multiplicity - 0) / 8);
    }

    style << "lines linecolor " << (multiplicity + 1) << " dashtype "
          << dashType;
    Line result(style.str(), title.str(), data.str());
    result.deduplicate_ = true;
    result.priority_ = multiplicity - 100;
    return result;
}

// Find the element whose neighbour is nearest to the given value
template <typename T>
typename std::map<double, T>::iterator findNearest(std::map<double, T>& map,
                                                   double value) {
    if (map.empty()) {
        return map.end();
    }

    // First iterator greater than value
    auto iter = map.upper_bound(value);
    if (iter == map.end()) {
        // No entry greater than value, so the last element is the nearest
        iter--;
        return iter;
    } else if (iter == map.begin()) {
        // value is smaller than the smallest entry in the map
        return iter;
    } else {
        // Actual comparison needed
        auto prev = iter;
        prev--;
        if (std::abs(iter->first - value) < std::abs(prev->first - value)) {
            return iter;
        } else {
            return prev;
        }
    }
}

// Compute the distance in between the key and the key of its nearest neighbour.
template <typename T>
double separation(const std::map<double, T>& map,
                  const typename std::map<double, T>::iterator iter) {
    if (map.size() <= 1) {
        return std::numeric_limits<double>::infinity();
    } else {
        double sep = std::numeric_limits<double>::infinity();
        if (iter != map.begin()) {
            auto prev = iter;
            prev--;
            sep = std::min(sep, std::abs(iter->first - prev->first));
        }
        auto next = iter;
        next++;
        if (next != map.end()) {
            sep = std::min(sep, std::abs(next->first - iter->first));
        }
        return sep;
    }
}

template <std::size_t DIM>
std::map<int, std::vector<std::tuple<double, double>>>
    BandstructureGNUPlot<DIM>::groupSpectrum() {
    const double TRESHOLD = 0.15;
    std::map<int, std::vector<std::tuple<double, double>>> points;
    std::size_t steps =
        computedSpectrum_->originalProblem().getPath().totalNumberOfSteps();
    double x = 0;
    for (std::size_t i = 0; i < steps; ++i) {
        if (i > 0) {
            x += computedSpectrum_->originalProblem().getPath().dk(i).l2Norm();
        }
        std::vector<double> freqs = computedSpectrum_->frequencies(i);
        // Theoretical spectrum for comparison
        std::map<double, std::size_t> theoretical = structure_.computeSpectrum(
            computedSpectrum_->originalProblem().getPath().k(i), freqs.size());
        // Deduplicate computed spectrum.
        std::map<double, std::size_t> computed;
        double last = 0;
        for (double f : freqs) {
            double sep = separation(theoretical, findNearest(theoretical, f));
            if (std::abs(f - last) < std::min(TRESHOLD, sep / 2) &&
                !computed.empty()) {
                (--computed.end())->second++;
            } else {
                computed[f] = 1;
                last = f;
            }
        }
        // Compute difference with theoretical multiplicity
        for (auto& cfreq : computed) {
            auto liter = findNearest(theoretical, cfreq.first);
            std::size_t theoreticalMultiplicity = 0;
            double sep = separation(theoretical, liter);
            if (liter != theoretical.end() &&
                std::abs(liter->first - cfreq.first) <
                    std::min(TRESHOLD, sep / 2)) {
                theoreticalMultiplicity = liter->second;
            }
            // Put it in the correct bin.
            points[cfreq.second - theoreticalMultiplicity].emplace_back(
                std::make_tuple(x, cfreq.first));
        }
    }
    return points;
}

template class BandstructureGNUPlot<2>;
template class BandstructureGNUPlot<3>;