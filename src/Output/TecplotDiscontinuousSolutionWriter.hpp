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
using std::string;
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
                const std::string& fileTitle,
                const std::string& dimensionsToWrite,
                const std::string& variableString);

        /// Write a zone with data from the current mesh to the stream held by the object.
        template <typename WriteFunctor>
        void write(const Base::MeshManipulator<DIM>* mesh,
                   const std::string& zoneTitle,
                   const bool sameGeometry,
                   WriteFunctor& writeDataFunc
                   );

	/// Write a zone with data from the current mesh to the stream held by the object. (class member write function)
	/// Has the exact same behaviour as the other write function except it uses OBJ::writeDataFunc() as the write function
        template <typename OBJ, typename WriteFunction>
        void write(const Base::MeshManipulator<DIM>* mesh,
                   const std::string& zoneTitle,
                   const bool sameGeometry,
                   WriteFunction& writeDataFunc,
                   OBJ* objPtr                   
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

        std::string elementType_[5];

        int* dimensionsToWrite_;

        const unsigned int nDimensionsToWrite_;
        
        unsigned int* dimNrs;
    };
}
#include "TecplotDiscontinuousSolutionWriter_Impl.hpp"
#endif
