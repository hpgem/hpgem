#include "TecplotDiscontinuousSolutionWriter.hpp"

namespace Output
{
    template<unsigned int DIM>
    TecplotDiscontinuousSolutionWriter<DIM>::TecplotDiscontinuousSolutionWriter(
            ostream& output,
            const StringT& fileTitle,
            int* dimensionsToWrite,
            const StringT& variableString):
                output_(output),
                previousNrOfElements_(0),
                previousNrOfNodes_(0),
                dimensionsToWrite_(dimensionsToWrite),
                nDimensionsToWrite_(sizeof(dimensionsToWrite)/sizeof(dimensionsToWrite[0]))
    {
        // element type for a given dimension as index
        elementType_[0] = "UNKNOWN";
        elementType_[1] = "FELINESEG";
        elementType_[2] = "FEQUADRILATERAL";
        elementType_[3] = "FEBRICK";
        elementType_[4] = "UNKNOWN";

        /// See if the dimensionsToWrite array makes sense.
        if (nDimensionsToWrite_ > DIM)
        {
            throw "TecplotDiscontinuousSolutionWriter: Printing more dimensions that there are";
        }
        for (int i = 0; i < nDimensionsToWrite_; ++i)
        {
            if (dimensionsToWrite_[0] > DIM)
                throw "TecplotDiscontinuousSolutionWriter: Requested an invalid dimension to print";
        }

        output_ << "TITLE = \"" << fileTitle << "\"\n";
        output_ << "VARIABLES = ";
        for (unsigned int i = 0; i < nDimensionsToWrite_; ++i)
        {
            output_ << "\"x" << dimensionsToWrite_[i] << "\", ";
        }
        output_ << makeTecplotVariableString(variableString) << "\n";
    }

    /*!
     * The variables should be the same as in the specification for the file given to the ctor.
     *
     * The function will iterate over the elements of the mesh m.
     *
     * The zoneTitle will go to the corresponding tecplot variable, e.g. this can be the time
     * in case several timesteps are put into one file (as different zones).
     *
     * If sameGeometry is true then the node coordinates will not be written but reused from
     * the first zone. Use this if the mesh is not changing.
     *
     * WriteFunctor must have an operator()(EType&, const Point<dim>&, ostream&); the point
     * given to it is in the coordinates of the reference element.
     *
     */
    //template <class WriteFunctor>
    template<unsigned int DIM>
    void TecplotDiscontinuousSolutionWriter<DIM>::write(
            const Base::MeshManipulator<DIM>& mesh,
            const StringT& zoneTitle,
            const bool sameGeometry//,
            //WriteFunctor& writeDataFunc
            )
    {
        long int posNumberOfNodes(0);
        long int posNumberOfElements(0);

        // Zone header.
        output_ << "ZONE T = \""
                << zoneTitle
                << "\""
                << ", ZONETYPE = "
                << elementType_[nDimensionsToWrite_]
                << ", DATAPACKING = POINT";
        output_ << ", N = ";

        if (sameGeometry)
        {
            // If the same geometry is used then we write the old counts for elements and nodes
            // into the header
            output_ << previousNrOfNodes_    << ", E = "
                    << previousNrOfElements_ << ", D = (";
            for (unsigned int ii = 1; ii <= nDimensionsToWrite_; ++ii)
            {
                output_ << ii << ",";
            }
            output_ << "FECONNECT)\n";
        }
        else
        {
            // If the geometry is new then we have to first leave the number of elements and
            // nodes open; however we keep the file positions for filling them in later in this
            // function
            posNumberOfNodes = output_.tellp();
            output_ << "          , E = ";
            posNumberOfElements = output_.tellp();
            output_ << "          \n";
        }

        unsigned int totalNrOfVertices = 0;
        unsigned int totalNrOfElements = 0;

        // Iterate over elements and write the solution at the vertices.
        // We do this by getting the element list from the mesh, and then iterating over the
        // elements.

        typedef Base::MeshManipulator<DIM> MeshType;
        typedef Base::Element<DIM>  ElementT;
            //typename typedef Base::MeshManipulator<DIM>::ListOfElementsT ListOfElementsT;
        typedef std::list<ElementT*> ListOfElementsT;

        Geometry::PointPhysical<DIM> pPhys;
        Geometry::PointReference<DIM> pRef;

        unsigned int nrOfNodes; // i.e. on one element
        TecplotPhysicalGeometryIterator& nodeIt = TecplotPhysicalGeometryIterator::Instance();

        const ListOfElementsT& elements = mesh.getElementsList();

        // 1. Element cycle, print physical coordinates.
        for (typename ListOfElementsT::const_iterator iterator = elements.begin(), end = elements.end();
                iterator != end;
                ++iterator)
        {
            totalNrOfElements++;
            nrOfNodes = 0;

            // Tell the TecplotPhysicalGeometryIterator which shape is to be iterated next

            nodeIt.acceptG((*iterator)->getPhysicalGeometry());

            // Cycle through nodes
            while (nodeIt.more())
            {
                const unsigned int localNode = nodeIt.getNodeNr();

                nrOfNodes++;

                if (!sameGeometry)
                {
                    // First write the (possibly reduced) coordinates of the point;
                    // note: PHYSICAL coordinates here!
                    (*iterator)->getPhysicalGeometry()->getNodeCoordinates(localNode, pPhys);

                    for (unsigned int i = 0; i < nDimensionsToWrite_; ++i)
                    {
                        output_.precision(6);
                        output_.width(12);
                        output_ << pPhys[dimensionsToWrite_[i]] << ' ';
                    }
                }

                // For the solution data, write function of the user, however we pass a local
                // coordinate of the current reference element
                (*iterator)->getReferenceGeometry()->getNode(localNode, pRef);

                // For safety we give more precision here, user may change it within the write
                // function
                output_.precision(8);
                output_.width(16);
                output_ << "0.0";
                //writeDataFunc(*it, pRef, output_); // TODO: ?
                output_ << "\n";

            } // 'nodes of element' loop

            /*
            if (((unsigned int) 1 << nrOfDimensionsWritten) != nrOfNodes)
            {
                throw "TecplotDiscontinuousSolutionWriter: wrong number of nodes written";
            }
            */
            totalNrOfVertices += nrOfNodes;
        }

        // 2. Print global node numbers per element.
        if (!sameGeometry)
        {
            unsigned int countNumberOfVertices = 1;

            // Node count PER ELEMENT (one global node can count several times locally).
            // (Basically we just write an ascending series of numbers).
            for (unsigned int elementCounter = 0; elementCounter < totalNrOfElements; elementCounter++)
            {
                for (unsigned int j = 0; j < ((unsigned int) 1 << nDimensionsToWrite_); ++j) // number of vertices
                {
                    output_ << countNumberOfVertices++ << " ";
                }
                output_ << "\n";
            }

            previousNrOfElements_ = totalNrOfElements;
            previousNrOfNodes_ = totalNrOfVertices;

            // Write the number of elements and nodes to the blank space at the beginning of the zone
            const long int currentFilePos = output_.tellp();
            output_.seekp(posNumberOfNodes);
            output_ << previousNrOfNodes_;
            output_.seekp(posNumberOfElements);
            output_ << previousNrOfElements_;
            output_.seekp(currentFilePos);
        }
        output_.flush();
    } // Write function


    template<unsigned int DIM>
    std::string TecplotDiscontinuousSolutionWriter<DIM>::makeTecplotVariableString(const std::string& s) const
    {
        std::string res(s);
        std::string::size_type pos = 0;
        const std::string::size_type one(1);
        const std::string replacement("\", \"");

        while ((pos = res.find(',', pos)) != std::string::size_type(-1))
        {
            res.replace(pos, one, replacement);
            pos += replacement.length();
        }

        res = "\"" + res + "\"";

        return res;
    }
}
