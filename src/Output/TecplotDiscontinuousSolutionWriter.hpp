/*
 * TecplotDiscontinuousSolutionWriter.hpp
 *
 *  Created on: Feb 16, 2013
 *      Author: nicorivas
 */
#ifndef TECPLOTDISCONTINUOUSSOLUTIONWRITER_HH
#define TECPLOTDISCONTINUOUSSOLUTIONWRITER_HH

#include <sstream>
#include <ostream>
using std::ostream;
#include <vector>
using std::vector;
#include <utility>
using std::pair;
#include <string>

// Package includes:
#include "Geometry/PointPhysical.hpp"
//using Geometry::PhysSpacePoint;
#include "Geometry/PointReference.hpp"
//using Geometry::RefSpacePoint;
#include "Base/MeshManipulator.hpp"
//using Geometry::Mesh;

#include "TecplotPhysicalGeometryIterator.hpp"

namespace Output
{
    template <unsigned int DIM>
    class TecplotDiscontinuousSolutionWriter
    {

        //! \brief This class prints the nodes and the solution in every element in Tecplot format.

    public:

        typedef std::string StringT;

    public:

        TecplotDiscontinuousSolutionWriter(
                ostream& output,
                const StringT& fileTitle,
                int* dimensionsToWrite,
                const StringT& variableString);

        /// Write a zone with data from the current mesh to the stream held by the object.
        //template <class WriteFunctor>
        void write(const Base::MeshManipulator<DIM>& mesh,
                   const StringT& zoneTitle,
                   const bool sameGeometry//,
                   //WriteFunctor& writeDataFunc
                   );

        /// TODO: Perfect this deconstructor.
        ~TecplotDiscontinuousSolutionWriter()
        {
            output_.flush();
        }

    private:

        StringT makeTecplotVariableString(const StringT& s) const;

        ostream& output_;

        unsigned int previousNrOfElements_;

        unsigned int previousNrOfNodes_;

        StringT elementType_[5];

        int* dimensionsToWrite_;

        const unsigned int nDimensionsToWrite_;
    };
}
#include "TecplotDiscontinuousSolutionWriter_Impl.hpp"
#endif
