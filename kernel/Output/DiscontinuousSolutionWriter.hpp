/* 
 * File:   discontinuousSolutionWriter.hpp
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
  class MeshManipulator;
  class Element;
}

namespace Geometry
{
  class PointReference;
}

namespace Output
{

  ///class to write the data of a single element.
  class SingleElementWriter
  {
  public:
    ///function that acutally writes the data for one node on one element.
    virtual void writeOutput(const Base::Element*, const Geometry::PointReference&, std::ostream&)=0;
  } ;

  /// \brief This class prints the solution in every element.
  ///
  /// Class to print the coordinates and solution in every node in every element.
  /// Example usage for Laplace problem: 
  /// DiscontinuousSolutionWriter testWriter(outStream, "title", "01", "u");
  /// testWriter.write(mesh);
  /// Then open gnuplot and plot with the command: splot "output.dat"
  class DiscontinuousSolutionWriter
  {
  public:

    typedef std::string StringT;

    ///Constructor: Initialise the output stream and write the header.
    DiscontinuousSolutionWriter(
                                std::ostream& output,
                                const std::string& fileTitle,
                                const std::string& dimensionsToWrite,
                                const std::string& resultVariableName);

    /// Write the data to the stream ouput_.
    void write(const Base::MeshManipulator* mesh, SingleElementWriter* writeDataClass);

    /// TODO: Perfect this deconstructor. Irana: there's no pointers, do we need better dtor?
    ~DiscontinuousSolutionWriter()
    {
      output_.flush();
    }

  private:

    ///stream where all output will be written to.
    std::ostream& output_;

    ///Number of physical dimensions of the domain of the problem.
    const unsigned int nDimensionsToWrite_;
  } ;
}

#endif	/* DISCONTINUOUSSOLUTIONWRITER_HH */

