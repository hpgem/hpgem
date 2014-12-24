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


#include <list>
#include <cassert>

#include "GNUPlotDiscontinuousSolutionWriter.hpp"
#include "Base/Element.hpp"
#include "TecplotPhysicalGeometryIterator.hpp"
#include "Base/MeshManipulator.hpp"
#include "Geometry/PointPhysical.hpp"
#include "Geometry/PointReference.hpp"
#include "Geometry/PhysicalGeometry.hpp"
#include "Geometry/ReferenceGeometry.hpp"
#include "Base/ElementCacheData.hpp"


namespace Output
{
    ///Constructor: initialise the output stream and make the header. 
    ///@param[in] output             The output stream you want to write to
    ///@param[in] fileTitle          The title of your file/plot
    ///@param[in] dimensionsToWrite  Names of the variables, for example "01" gives x0 and x1
    ///@param[in] resultVariableName Names of the dependent variable, for example "position" or "u"
    GNUPlotDiscontinuousSolutionWriter::GNUPlotDiscontinuousSolutionWriter(
                                                             std::ostream& output,
                                                             const std::string& fileTitle,
                                                             const std::string& dimensionsToWrite,
                                                             const std::string& resultVariableName) :
    output_(output),
    nDimensionsToWrite_(dimensionsToWrite.length())
    {
        output_ << "TITLE = \"" << fileTitle << "\"\n";
        output_ << "VARIABLES = ";
        for (unsigned int i = 0; i < nDimensionsToWrite_; ++i)
        {
            output_ << "\"x" << dimensionsToWrite[i] << "\", ";
        }
        output_ << "\"" << resultVariableName << "\"" << "\n";
    }
    
    
    /**
     * Function that writes the values of the physical coordinates and the values
     * of the results to the stream output_. It does this by iterating over the elements
     * and the nodes in each element, find out its physical coordinates and then 
     * call the function writeToFile which has to be defined in the writeDataClass.
     * @param[in] mesh           Mesh on which the values are already computed
     * @param[in] writeDataClass Class which is a child of the class SingleElementWriter
     *                           which has defined the function writeToFile.
     */
    void GNUPlotDiscontinuousSolutionWriter::write(const Base::MeshManipulator* mesh,
                                            SingleElementWriter* writeDataClass)
    {
        
        typedef Base::Element ElementT;
        typedef std::vector<ElementT*> ListOfElementsT;

        //First assert that we have defined the correct number of dimensions in 
        //the constructor.
        assert(mesh->dimension() == nDimensionsToWrite_);
        
        //Write how many elements there are in this simulation.
        //It would be nice to also display the polynomial order here.
        output_ << "Number of elements = " << mesh->getNumberOfElements() << '\n';


        // Iterate over elements and write the solution at the vertices.
        // We do this by getting the element list from the mesh, and then iterating over the
        // elements.
        
        
        //(@Irana) Does GNUPlot se the same node ordering as Tecplot, for all element types? -FB
        //make an iterator that can iterate over all nodes.
        TecplotPhysicalGeometryIterator& nodeIt = TecplotPhysicalGeometryIterator::Instance();

        //construct the list of all elements.
        const ListOfElementsT& elements = mesh->getElementsList();

        //get the physical and reference coordinates of the first node of the first element.
        Geometry::PointPhysical pPhys((*elements.begin())->getPhysicalGeometry()->getNodePtr(0)->size());
        Geometry::PointReference pRef((*elements.begin())->getPhysicalGeometry()->getNodePtr(0)->size());

        //Element cycle, print physical coordinates:
        for (typename ListOfElementsT::const_iterator eltIterator = elements.begin();
            eltIterator != elements.end(); ++eltIterator)
        {
            // Tell the TecplotPhysicalGeometryIterator which shape is to be iterated next
            nodeIt.acceptG((*eltIterator)->getPhysicalGeometry());

            // Cycle through nodes
            while (nodeIt.more())
            {
                const unsigned int localNode = nodeIt.getNodeNr();

                // For the solution data, write function of the user, however we pass a local
                // coordinate of the current reference element
                (*eltIterator)->getReferenceGeometry()->getNode(localNode, pRef);

                // First write the (possibly reduced) coordinates of the point;
                // note: PHYSICAL coordinates here!
                (*eltIterator)->referenceToPhysical(pRef, pPhys);

                //write the physical coordinates of the point
                for (unsigned int i = 0; i < nDimensionsToWrite_; ++i)
                {
                    output_.precision(6);
                    output_.width(12);
                    output_ << pPhys[i] << ' ';
                }

                // For safety we give more precision here, user may change it within the write
                // function
                output_.precision(8);
                output_.width(16);
                writeDataClass->writeOutput(*eltIterator, pRef, output_);
                output_ << "\n";

            } // end of 'nodes of element' loop
        } // end of 'elements' loop
        output_.flush();
    } // end of write function 


}
