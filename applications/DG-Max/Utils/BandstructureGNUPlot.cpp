
#include "BandstructureGNUPlot.h"
#include "Logger.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <set>

template <std::size_t DIM>
BandstructureGNUPlot<DIM>::BandstructureGNUPlot(
        const KSpacePath<DIM> &path,
        const std::vector<std::string>& pointNames,
        const HomogeneousBandStructure<DIM>& structure,
        const BaseEigenvalueResult<DIM>* computedSpectrum)
    : path_ (path)
    , pointNames_ (pointNames)
    , structure_ (structure)
    , computedSpectrum_ (computedSpectrum)
{
    logger.assert_always(path_.numberOfCornerPoints() == pointNames_.size()
        || pointNames_.empty(),
        "Mismatch in path length and number of point names.");
    // It might even be better to check if the two paths match, or take the path
    // from the computed spectrum.
    logger.assert_always(computedSpectrum_ == nullptr
        || computedSpectrum_->originalProblem().getPath().totalNumberOfSteps() == path.totalNumberOfSteps(),
        "Computed and path step length mismatch.");
}

template <std::size_t DIM>
void BandstructureGNUPlot<DIM>::plot(std::string fileName)
{
    const double maxFreq = 16;

    std::ofstream out (fileName);
    // Header
    out << "set term postscript eps color enhanced" << std::endl;
    out << "set output 'bands.eps'" << std::endl;

    std::vector<double> xpoints (path_.numberOfCornerPoints());
    double totalDistance = 0;
    xpoints[0] = 0;
    for(std::size_t i = 1; i < path_.numberOfCornerPoints(); ++i)
    {
        totalDistance += (path_.kcorner(i-1) - path_.kcorner(i)).l2Norm();
        xpoints[i] = totalDistance;
    }
    // Set x,y range
    out << "xmin=" << 0.0 << ";xmax=" << totalDistance << std::endl;
    out << "ymin=" << 0.0  << ";ymax=" << maxFreq << std::endl;
    out << "set xrange [xmin:xmax]" << std::endl;
    out << "set yrange [ymin:ymax]" << std::endl;
    // Move the legend/key to the bottom, reducing overlap
    out << "set key bottom" <<  std::endl;

    // Set the x-ticks
    if(!pointNames_.empty())
    {
        out << "set xtics (";
        for(std::size_t i = 0; i < path_.numberOfCornerPoints(); ++i)
        {
            if(i != 0)
            {
                out << ", ";
            }
            out << "'" << pointNames_[i] << "'" << xpoints[i];
        }
        out << ")" << std::endl;
    }

    // Frequency data
    std::map<int, std::vector<std::tuple<double, double>>> points;
    if(computedSpectrum_ != nullptr)
    {
        points = groupSpectrum();

        out << "$bandpoints << EOD" << std::endl;
        for(auto& ps : points)
        {
            for(auto p : ps.second)
            {
                out << std::get<0>(p) << " " << std::get<1>(p) << std::endl;
            }
            // Dataset separation
            out << std::endl;
            out << std::endl;
        }

        out << "EOD" << std::endl;
    }


    // Plot the actual structure
    std::map<std::size_t, std::vector<std::string>> results;
    std::set<std::size_t> hasKey;
    out << "plot sample ";
    for(std::size_t i = 1; i < path_.numberOfCornerPoints(); ++i)
    {
        std::vector<typename HomogeneousBandStructure<DIM>::Line> lines =
                structure_.computeLines(path_.kcorner(i-1), path_.kcorner(i), maxFreq);
        for(auto& line : lines)
        {
            std::size_t multiplicity = line.multiplicity();
            std::vector<std::string>& resultBin = results[multiplicity];
            bool title = resultBin.empty();
            resultBin.emplace_back(band(line, xpoints[i-1], xpoints[i], title));
        }
    }
    bool first = true;
    for(auto& bin : results)
    {
        for(auto& line : bin.second)
        {
            if(first)
            {
                first = false;
            }
            else
            {
                out << ",\\" << std::endl << "  ";
            }
            out << line;
        }
    }
    if(computedSpectrum_ != nullptr)
    {
        int index = 0;
        for(auto& ps : points)
        {
            int difference = ps.first;
            int pointType;
            if (difference < 0)
                pointType = 2; //'x'
            else if (difference > 0)
                pointType = 1; // '+'
            else
                pointType = 7; // 'o'

            out << ",\\" << std::endl
                << "  $bandpoints index " << index
                << " with points pointtype " << pointType;
            if(difference == 0)
            {
                // the points are rather big for something that is correct so
                // reduce their size a bit.
                out << " pointsize 0.3";
            }
            out << " linecolor " << std::abs(difference) + 1
                << " title \"" << difference << "\"";
            index++;
        }
    }
    out << std::endl;
    out.close();
}

template<std::size_t DIM>
std::string BandstructureGNUPlot<DIM>::band(
        const typename HomogeneousBandStructure<DIM>::Line &line,
        double x1, double x2, bool titled)
{
    double x = line.getX();
    double y = line.getY();
    double l = line.getL();
    std::ostringstream out;

    // print the range
    out << "[" << x1 << ":" << x2 << "]";
    // actual function
    // - xpoints to correct for the offset on the x-axis,
    // +x from the actual shape
    out << " sqrt((x + " << (x - x1) << ")**2 + " << y*y << ")";
    std::size_t multiplicity = line.multiplicity();

    if(!titled)
    {
        out << " notitle";
    }
    else
    {
        out << " title \"" << multiplicity << '"';
    }

    // mult + 1 so 0 can also have its unique color, which is used in the
    // computed spectrum.
    out << " linecolor " << (multiplicity + 1) << " dashtype " << (1 + (multiplicity - 1) / 8);
    return out.str();
}

// Find the element whose neighbour is nearest to the given value
template<typename T>
typename std::map<double, T>::iterator findNearest (
        std::map<double, T>& map, double value)
{
    if(map.empty())
    {
        return map.end();
    }

    // First iterator greater than value
    auto iter = map.upper_bound(value);
    if (iter == map.end())
    {
        // No entry greater than value, so the last element is the nearest
        iter--;
        return iter;
    }
    else if (iter == map.begin())
    {
        // value is smaller than the smallest entry in the map
        return iter;
    }
    else
    {
        // Actual comparison needed
        auto prev = iter;
        prev--;
        if (std::abs(iter->first - value) < std::abs(prev->first - value))
        {
            return iter;
        }
        else
        {
            return prev;
        }
    }
}

// Compute the distance in between the key and the key of its nearest neighbour.
template<typename T>
double separation (const std::map<double, T>& map, const typename std::map<double, T>::iterator iter)
{
    if(map.size() <= 1)
    {
        return std::numeric_limits<double>::infinity();
    }
    else
    {
        double sep = std::numeric_limits<double>::infinity();
        if(iter != map.begin())
        {
            auto prev = iter;
            prev--;
            sep = std::min(sep, std::abs(iter->first - prev->first));
        }
        auto next = iter;
        next++;
        if(next != map.end())
        {
            sep = std::min(sep, std::abs(next->first - iter->first));
        }
        return sep;
    }

}


template<std::size_t DIM>
std::map<int, std::vector<std::tuple<double, double>>> BandstructureGNUPlot<DIM>::groupSpectrum()
{
    const double TRESHOLD = 0.15;
    std::map<int, std::vector<std::tuple<double, double>>> points;
    std::size_t steps = computedSpectrum_->originalProblem().getPath().totalNumberOfSteps();
    double x = 0;
    for(std::size_t i = 0; i < steps; ++i)
    {
        if (i > 0)
        {
            x += computedSpectrum_->originalProblem().getPath().dk(i).l2Norm();
        }
        std::vector<double> freqs = computedSpectrum_->frequencies(i);
        // Theoretical spectrum for comparison
        std::map<double, std::size_t> theoretical = structure_.computeSpectrum(
                computedSpectrum_->originalProblem().getPath().k(i),
                freqs.size());
        // Deduplicate computed spectrum.
        std::map<double, std::size_t> computed;
        double last = 0;
        for(double f : freqs)
        {
            double sep = separation(theoretical, findNearest(theoretical, f));
            if(std::abs(f - last) < std::min(TRESHOLD, sep/2) && !computed.empty())
            {
                (--computed.end())->second++;
            }
            else
            {
                computed[f] = 1;
                last = f;
            }
        }
        // Compute difference with theoretical multiplicity
        for(auto& cfreq : computed)
        {
            auto liter = findNearest(theoretical, cfreq.first);
            std::size_t theoreticalMultiplicity = 0;
            double sep = separation(theoretical, liter);
            if(liter != theoretical.end()
                && std::abs(liter->first - cfreq.first) < std::min(TRESHOLD, sep/2))
            {
                theoreticalMultiplicity = liter->second;
            }
            // Put it in the correct bin.
            points[cfreq.second - theoreticalMultiplicity]
                    .emplace_back(std::make_tuple(x, cfreq.first));
        }
    }
    return points;
}

template class BandstructureGNUPlot<2>;
template class BandstructureGNUPlot<3>;