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

#include "TecplotDiscontinuousSolutionWriter.h"
#include "Base/Element.h"
#include "TecplotPhysicalGeometryIterator.h"
#include "Base/MeshManipulator.h"
#include "Geometry/PointPhysical.h"
#include "Geometry/PointReference.h"
#include "Geometry/PhysicalGeometry.h"
#include "Geometry/ReferenceGeometry.h"
#include "TecplotSingleElementWriter.h"
#include "Base/ElementCacheData.h"
#include "Base/FaceCacheData.h"

namespace Output
{
    
    TecplotDiscontinuousSolutionWriter::TecplotDiscontinuousSolutionWriter(std::ostream& output, const std::string& fileTitle, const std::string& dimensionsToWrite, const std::string& variableString)
            : output_(output), previousNrOfElements_(0), previousNrOfNodes_(0), nDimensionsToWrite_(dimensionsToWrite.length())
    {
        logger.assert_always(output.good(), "Something is not so good about the given output stream");
        // element type for a given dimension as index
        elementType_[0] = "UNKNOWN";
        elementType_[1] = "FELINESEG";
        elementType_[2] = "FEQUADRILATERAL";
        elementType_[3] = "FEBRICK";
        elementType_[4] = "UNKNOWN";
        
        dimNrs = new std::size_t[nDimensionsToWrite_];
        for (std::size_t i = 0; i < nDimensionsToWrite_; ++i)
        {
            std::istringstream istr(dimensionsToWrite.substr(i, 1));
            istr >> dimNrs[i];
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
    void TecplotDiscontinuousSolutionWriter::write(const Base::MeshManipulator* mesh, const std::string& zoneTitle, const bool sameGeometry, TecplotSingleElementWriter* writeDataClass, const double time)
    {
        logger.assert(mesh!=nullptr, "Invalid mesh passed to this writer");
        logger.assert(writeDataClass!=nullptr, "Invalid write class passed");
        std::function<void(const Base::Element*, const Geometry::PointReference&, std::ostream&)> function = 
        [=](const Base::Element* el, const Geometry::PointReference& pR, std::ostream& os){
            writeDataClass->writeToTecplotFile(el,pR,os);
        };
        write(mesh,zoneTitle,sameGeometry,function,time);
        
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
    void TecplotDiscontinuousSolutionWriter::write(const Base::MeshManipulator* mesh, const std::string& zoneTitle, const bool sameGeometry, std::function<void(const Base::Element*, const Geometry::PointReference&, std::ostream&)>writeDataFun, const double time)
    {
        logger.assert(mesh!=nullptr, "Invalid mesh passed to this writer");
        //assertion is technically checking internal state, but the writability of the filesystem may change outside the influence of this class
        logger.assert_always(output_.good(), "Something is not so good about the output stream");
        
        std::streamoff posNumberOfNodes(0);
        std::streamoff posNumberOfElements(0);
        
        // Zone header.
        output_ << "ZONE T = \"" << zoneTitle << "\"" << ", STRANDID = 1" << ", SOLUTIONTIME = " << time << ", ZONETYPE = " << elementType_[nDimensionsToWrite_] << ", DATAPACKING = POINT";
        output_ << ", N = ";
        
        if (sameGeometry)
        {
            // If the same geometry is used then we write the old counts for elements and nodes
            // into the header
            output_ << previousNrOfNodes_ << ", E = " << previousNrOfElements_ << ", D = (";
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
        
        using MeshType = Base::MeshManipulator;
        using ElementT = Base::Element;
        using ListOfElementsT = std::vector<ElementT*>;
        
        std::size_t nrOfNodes; // i.e. on one element
        TecplotPhysicalGeometryIterator& nodeIt = TecplotPhysicalGeometryIterator::Instance();
        
        const ListOfElementsT& elements = mesh->getElementsList();
        
        Geometry::PointPhysical pPhys((*elements.begin())->getPhysicalGeometry()->getNodePtr(0)->size());
        Geometry::PointReference pRef((*elements.begin())->getPhysicalGeometry()->getNodePtr(0)->size());
        
        // 1. Element cycle, print physical coordinates.
        
        for (typename ListOfElementsT::const_iterator iterator = elements.begin(), end = elements.end(); iterator != end; ++iterator)
        {
            totalNrOfElements++;
            nrOfNodes = 0;
            // Tell the TecplotPhysicalGeometryIterator which shape is to be iterated next
            nodeIt.acceptG((*iterator)->getPhysicalGeometry());
            
            // Cycle through nodes
            while (nodeIt.more())
            {
                const std::size_t localNode = nodeIt.getNodeNr();
                
                nrOfNodes++;
                
                // For the solution data, write function of the user, however we pass a local
                // coordinate of the current reference element
                pRef = (*iterator)->getReferenceGeometry()->getNode(localNode);
                
                if (!sameGeometry)
                {
                    // First write the (possibly reduced) coordinates of the point;
                    // note: PHYSICAL coordinates here!                    
                    pPhys = (*iterator)->referenceToPhysical(pRef);
                    
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
                writeDataFun(*iterator, pRef, output_);
                output_ << "\n";
                
            } // 'nodes of element' loop
            
            
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
                for (std::size_t j = 0; j < ( 1UL << nDimensionsToWrite_); ++j) // number of vertices
                {
                    output_ << countNumberOfVertices++ << " ";
                }
                output_ << "\n";
            }
            
            previousNrOfElements_ = totalNrOfElements;
            previousNrOfNodes_ = totalNrOfVertices;
            
            // Write the number of elements and nodes to the blank space at the beginning of the zone
            const std::streamoff currentFilePos = output_.tellp();
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
