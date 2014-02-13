/*
 * TecplotDiscontinuousSolutionWriter_Impl.hpp
 *
 *  Created on: Feb 3, 2014
 *      Author: brinkf
 */

#ifndef TECPLOTDISCONTINUOUSSOLUTIONWRITER_IMPL_HPP_
#define TECPLOTDISCONTINUOUSSOLUTIONWRITER_IMPL_HPP_

namespace Output{




   /* template <typename OBJ, typename WriteFunction>
    void TecplotDiscontinuousSolutionWriter::write(
            const Base::MeshManipulator* mesh,
            const std::string& zoneTitle,
            const bool sameGeometry,
            WriteFunction& writeDataFunc,
	    OBJ* objPtr)
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

        typedef Base::MeshManipulator MeshType;
        typedef Base::Element  ElementT;
            //typename typedef Base::MeshManipulator<DIM>::ListOfElementsT ListOfElementsT;
        typedef std::list<ElementT*> ListOfElementsT;

        Geometry::PointPhysical pPhys(nDimensionsToWrite_);
        Geometry::PointReference pRef(nDimensionsToWrite_);

        unsigned int nrOfNodes; // i.e. on one element
        TecplotPhysicalGeometryIterator& nodeIt = TecplotPhysicalGeometryIterator::Instance();

        const ListOfElementsT& elements = mesh->getElementsList();

        // 1. Element cycle, print physical coordinates.

       // mesh->outputMesh(std::cout);
        for (typename ListOfElementsT::const_iterator iterator = elements.begin(), end = elements.end();
                iterator != end;
                ++iterator)
        {
            totalNrOfElements++;
            nrOfNodes = 0;
            // Tell the TecplotPhysicalGeometryIterator which shape is to be iterated next
            const Base::Element&   el = *(*iterator);
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
                        output_ << pPhys[dimNrs[i]] << ' ';
                    }
                }

                // For the solution data, write function of the user, however we pass a local
                // coordinate of the current reference element
                (*iterator)->getReferenceGeometry()->getNode(localNode, pRef);

                // For safety we give more precision here, user may change it within the write
                // function
                output_.precision(8);
                output_.width(16);
                    //output_ << "0.0";
                (objPtr->*writeDataFunc)(**iterator, pRef, output_);
                output_ << "\n";

            } // 'nodes of element' loop

            /*
            if (((unsigned int) 1 << nrOfDimensionsWritten) != nrOfNodes)
            {
                throw "TecplotDiscontinuousSolutionWriter: wrong number of nodes written";
            }
            *//*
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
    } // Write function*/
}



#endif /* TECPLOTDISCONTINUOUSSOLUTIONWRITER_IMPL_HPP_ */
