/* 
 * File:   discontinuousSolutionWriter.h
 * Author: irana
 *
 * Created on November 10, 2014, 2:58 PM
 */

#ifndef DISCONTINUOUSSOLUTIONWRITER_HH
#define	DISCONTINUOUSSOLUTIONWRITER_HH

#include <ostream>
#include <vector>

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
    ///class to write the data of a single element.
    template<std::size_t DIM>
    class SingleElementWriter
    {
    public:
        ///function that actually writes the data for one node on one element.
        /// it is purely virtual, since there is no default for what needs to be
        /// written.
        virtual void writeOutput(const Base::Element*, const Geometry::PointReference<DIM>&, std::ostream&) = 0;
    };
    
    /// \brief This class prints the solution in every element.
    ///
    /// Class to print the coordinates and solution in every node in every element.
    /// Example usage for Laplace problem: 
    /// DiscontinuousSolutionWriter testWriter(outStream, "title", "01", "u");
    /// testWriter.write(mesh);
    /// Then open gnuplot and plot with the command: plot "output.dat" every ::3 (1D) or splot "output.dat" (for 2D)
    template<std::size_t DIM>
    class GNUPlotDiscontinuousSolutionWriter
    {
    public:
        
        ///Constructor: Initialise the output stream and write the header.
        GNUPlotDiscontinuousSolutionWriter(std::ostream& output, const std::string& fileTitle, const std::string& dimensionsToWrite, const std::string& resultVariableName);
        
        ///No copy constructor, since we don't want to assign the same ostream to different writers at the same time.
        GNUPlotDiscontinuousSolutionWriter(const GNUPlotDiscontinuousSolutionWriter &other) = delete;

        /// Write the data to the stream ouput_.
        void write(const Base::MeshManipulator<DIM>* mesh, SingleElementWriter<DIM>* writeDataClass);

        ///Destructor: just flush the output stream, the rest will be destructed 
        /// automatically.
        ~GNUPlotDiscontinuousSolutionWriter()
        {
            output_.flush();
        }
        
    private:
        
        ///stream where all output will be written to.
        std::ostream& output_;

        ///Number of physical dimensions of the domain of the problem.
        const std::size_t nDimensionsToWrite_;
    };
}

#include "GNUPlotDiscontinuousSolutionWriter_Impl.h"

#endif	/* DISCONTINUOUSSOLUTIONWRITER_HH */

