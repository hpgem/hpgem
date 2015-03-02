/*
 This file forms part of hpGEM. This package has been developed over a number of years by various people at the University of Twente and a full list of contributors can be found at
 http://hpgem.org/about-the-code/team
 
 This code is distributed using BSD 3-Clause License. A copy of which can found below.
 
 
 Copyright (c) 2014, University of Twente
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 
 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "TecplotDiscontinuousSolutionWriter.hpp"
#include "Base/Element.hpp"
#include "TecplotPhysicalGeometryIterator.hpp"
#include "Base/MeshManipulator.hpp"
#include "Geometry/PointPhysical.hpp"
#include "Geometry/PointReference.hpp"
#include "Geometry/PhysicalGeometry.hpp"
#include "Geometry/ReferenceGeometry.hpp"
#include "TecplotSingleElementWriter.hpp"
#include "Base/ElementCacheData.hpp"
#include "Base/FaceCacheData.hpp"

#include <list>

namespace Output
{
    
    TecplotDiscontinuousSolutionWriter::TecplotDiscontinuousSolutionWriter(
    		std::ostream& output, const std::string& fileTitle,
            const std::string& dimensionsToWrite,
            const std::string& variableString):
                output_(output),
                previousNrOfElements_(0),
                previousNrOfNodes_(0),
                nDimensionsToWrite_(dimensionsToWrite.length())
    {
        // element type for a given dimension as index
        elementType_[0] = "UNKNOWN";
        elementType_[1] = "FELINESEG";
        elementType_[2] = "FEQUADRILATERAL";
        elementType_[3] = "FEBRICK";
        elementType_[4] = "UNKNOWN";

        /// See if the dimensionsToWrite array makes sense.\todo think of a way to check this
        /*if (nDimensionsToWrite_ > DIM)
        {
            throw "TecplotDiscontinuousSolutionWriter: Printing more dimensions that there are";
        }*/
        dimNrs = new std::size_t[nDimensionsToWrite_];
        for (std::size_t i = 0; i < nDimensionsToWrite_; ++i)
        {
            std::istringstream istr(dimensionsToWrite.substr(i, 1));
		    istr >> dimNrs[i];
                //            if (dimensionsToWrite_[0] > DIM)
                //throw "TecplotDiscontinuousSolutionWriter: Requested an invalid dimension to print";
        }

        output_ << "TITLE = \"" << fileTitle << "\"\n";
        output_ << "VARIABLES = ";
        for (std::size_t i = 0; i < nDimensionsToWrite_; ++i)
        {
            output_ << "\"x" << dimNrs[i] << "\", ";
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
     * WriteFunctor(EType&, const Point<dim>&, ostream&) can also be a function
     *
     * Setting the variable time enables you to create animations in TecPlot.
     */
    void TecplotDiscontinuousSolutionWriter::write(
                const Base::MeshManipulator* mesh,
                const std::string& zoneTitle,
                const bool sameGeometry,
                TecplotSingleElementWriter* writeDataClass,
                const double time)
        {

            std::size_t posNumberOfNodes(0);
            std::size_t posNumberOfElements(0);

            // Zone header.
            output_ << "ZONE T = \""
                    << zoneTitle
                    << "\""
                    << ", STRANDID = 1"
                    << ", SOLUTIONTIME = "
                    << time
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
                for (std::size_t ii = 1; ii <= nDimensionsToWrite_; ++ii)
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

            std::size_t totalNrOfVertices = 0;
            std::size_t totalNrOfElements = 0;

            // Iterate over elements and write the solution at the vertices.
            // We do this by getting the element list from the mesh, and then iterating over the
            // elements.

            typedef Base::MeshManipulator MeshType;
            typedef Base::Element  ElementT;
                //typename typedef Base::MeshManipulator<DIM>::ListOfElementsT ListOfElementsT;
            typedef std::vector<ElementT*> ListOfElementsT;

            std::size_t nrOfNodes; // i.e. on one element
            TecplotPhysicalGeometryIterator& nodeIt = TecplotPhysicalGeometryIterator::Instance();

            const ListOfElementsT& elements = mesh->getElementsList();

            Geometry::PointPhysical pPhys((*elements.begin())->getPhysicalGeometry()->getNodePtr(0)->size());
            Geometry::PointReference pRef((*elements.begin())->getPhysicalGeometry()->getNodePtr(0)->size());

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
                    const std::size_t localNode = nodeIt.getNodeNr();

                    nrOfNodes++;

                    // For the solution data, write function of the user, however we pass a local
                    // coordinate of the current reference element
                    (*iterator)->getReferenceGeometry()->getNode(localNode, pRef);

                    if (!sameGeometry)
                    {
                        // First write the (possibly reduced) coordinates of the point;
                        // note: PHYSICAL coordinates here!

                    	(*iterator)->referenceToPhysical(pRef,pPhys);
                        //(*iterator)->getPhysicalGeometry()->getNodeCoordinates(localNode, pPhys);

                        for (std::size_t i = 0; i < nDimensionsToWrite_; ++i)
                        {
                            output_.precision(6);
                            output_.width(12);
                            output_ << pPhys[dimNrs[i]] << ' ';
                        }
                    }

                    // For safety we give more precision here, user may change it within the write
                    // function
                    output_.precision(8);
                    output_.width(16);
                        //output_ << "0.0";
                    writeDataClass->writeToTecplotFile(*iterator, pRef, output_);
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
                std::size_t countNumberOfVertices = 1;

                // Node count PER ELEMENT (one global node can count several times locally).
                // (Basically we just write an ascending series of numbers).
                for (std::size_t elementCounter = 0; elementCounter < totalNrOfElements; elementCounter++)
                {
                    for (std::size_t j = 0; j < ((std::size_t) 1 << nDimensionsToWrite_); ++j) // number of vertices
                    {
                        output_ << countNumberOfVertices++ << " ";
                    }
                    output_ << "\n";
                }

                previousNrOfElements_ = totalNrOfElements;
                previousNrOfNodes_ = totalNrOfVertices;

                // Write the number of elements and nodes to the blank space at the beginning of the zone
                const std::size_t currentFilePos = output_.tellp();
                output_.seekp(posNumberOfNodes);
                output_ << previousNrOfNodes_;
                output_.seekp(posNumberOfElements);
                output_ << previousNrOfElements_;
                output_.seekp(currentFilePos);
            }
            output_.flush();
        } // Write function

    std::string TecplotDiscontinuousSolutionWriter::makeTecplotVariableString(const std::string& s) const
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
