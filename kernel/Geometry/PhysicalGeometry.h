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

#ifndef PHYSICALGEOMETRY_H_
#define PHYSICALGEOMETRY_H_

#include <vector>
#include <iostream>

#include "PointPhysical.h"
#include "PointReference.h"
#include "ReferenceGeometry.h"

namespace Geometry
{
    class ReferenceGeometry;
    class PointPhysical;
    
    class PhysicalGeometry
    {
        /*!\class PhysicalGeometry
         * \brief PhysicalGeometry describes an actual physical shape in real space.
         * \details
         * You shouldn't create a particular Physical<Shape> (although it is possible), but a
         * PhysicalGeometry. The reference geometry passed to the constructor selects what shape
         * this PhysicalGeometry has
         *
         * It contains only the global indexes of its points in globalNodeIndexes_. These global
         * indexes refer to the global node container, of which every PhysicalGeometry has a
         * reference: nodes_.
         *
         * It also contains a reference to the corresponding referenceGeometry.
         *
         * ~ Point is the name of a class.               ~
         * ~ Node is a point that belongs to a geometry. ~
         */
    public:
        
        using VectorOfPhysicalPointsT = std::vector<PointPhysical>;

    public:
        
        /// \brief Constructor gets indexes of the nodes, a reference to the node container, and a pointer to the corresponding reference geometry.
        
        PhysicalGeometry(const std::vector<std::size_t>& globalNodeIndexes, VectorOfPhysicalPointsT& nodes, const ReferenceGeometry * const refG)
                : nodes_(nodes), globalNodeIndexes_(globalNodeIndexes), refGeometry_(refG)
        {
            logger.assert(refG!=nullptr, "Invalid reference geometry passed");
        }
                
        PhysicalGeometry(const PhysicalGeometry& other) = delete;
        
        /// \brief Returns a constant pointer to the container of the global node indexes.
        const std::vector<std::size_t>& getNodeIndexes() const
        {
            return globalNodeIndexes_;
        }
        
        /// \brief Returns a pointer to the global container of nodes.
        VectorOfPhysicalPointsT& getNodes()
        {
            return nodes_;
        }
        
        /// \brief Returns a constant pointer of the global container of nodes
        const VectorOfPhysicalPointsT& getNodes() const
        {
            return nodes_;
        }
        
        /// \brief Returns the name of the particular geometry.
        std::string getName()
        {
            //skip "Reference" and put "Physical" instead
            return std::string("Physical").append(refGeometry_->getName().substr(9));
        }

        /// \brief Given a local index relative to globalNodeIndexes_, return the global node index.
        std::size_t getNodeIndex(std::size_t localIndex) const
        {
            logger.assert(localIndex < getNumberOfNodes(), "Asked for local index %, but this geometry only has % nodes",localIndex,getNumberOfNodes());
            return globalNodeIndexes_[localIndex];
        }
        
        /// \brief Given a global index, returns a pointer to the corresponding point.
        const PointPhysical* getNodePtr(const std::size_t globalIndex) const
        {
            logger.assert(globalIndex < nodes_.size(),"This mesh does not contain a node with index %",globalIndex);
            return &(nodes_[globalIndex]);
        }
        
        /// \brief Given a global index, returns a pointer to the corresponding point.
        PointPhysical* getNodePtr(const std::size_t globalIndex)
        {
            logger.assert(globalIndex < nodes_.size(),"This mesh does not contain a node with index %",globalIndex);
            return &(nodes_[globalIndex]);
        }
        
        /// \brief Returns the number of nodes of this geometry.
        std::size_t getNumberOfNodes() const
        {
            return globalNodeIndexes_.size();
        }
        
        /// \brief Given a local index, return the physical coordinates of the corresponding point.
        PointPhysical getLocalNodeCoordinates(const std::size_t localIndex) const;

        /// \brief Given a global index, return the physical coordinates of the corresponding point.
        PointPhysical getGlobalNodeCoordinates(const std::size_t globalIndex) const;

        /// \brief Given a local face index, return the global indices of the entities contained on that face.
        std::vector<std::size_t> getGlobalFaceNodeIndices(const std::size_t i) const
        {
            logger.assert(i < getNrOfFaces(), "Asked for face %, but there are only % faces", i, getNrOfFaces());
            std::vector<std::size_t> result = getLocalFaceNodeIndices(i);
            for(std::size_t i=0; i < result.size(); ++i)
            {
                result[i] = getNodeIndex(result[i]);
            }
            return result;
        }

        /// \brief Given a local face index, return the local indices of the entities contained on that face.
        std::vector<std::size_t> getLocalFaceNodeIndices(const std::size_t i) const
        {
            logger.assert(i < getNrOfFaces(), "Asked for face %, but there are only % faces", i, getNrOfFaces());
            return refGeometry_->getCodim1EntityLocalIndices(i);
        }

        /// \brief Returns the number of faces via a call to ReferenceGeometry->getNrOfCodim1Entities();
        std::size_t getNrOfFaces() const
        {
            return refGeometry_->getNrOfCodim1Entities();
        }

        /// \brief Returns a reference to the corresponding reference geometry.
        const ReferenceGeometry * const getRefGeometry() const
        {
            return refGeometry_;
        }
        
        /// \brief Output operator
        friend std::ostream& operator<<(std::ostream& os, const PhysicalGeometry& physicalGeometry)
        {
            os << "PhysicalGeometry=( ";
            
            for (std::size_t i = 0; i < physicalGeometry.getNumberOfNodes(); i++)
            {
                os << physicalGeometry.getNodeIndex(i) << " ";
            }
            os << ')' << std::endl;
            
            return os;
        }
        
    protected:
        /// Reference to the global node container.
        VectorOfPhysicalPointsT& nodes_;

        /// Reference to the container of global indexes of the nodes, relative to nodes_.
        std::vector<std::size_t> globalNodeIndexes_;

        const ReferenceGeometry * const refGeometry_;
    };

}
#endif /* PHYSICALGEOMETRY_H_ */
