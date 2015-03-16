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

#ifndef REFERENCETOPHYSICALM_H_
#define REFERENCETOPHYSICALM_H_

#include "MappingInterface.h"
#include "Logger.h"
#include <vector>

namespace Geometry
{
    class PhysicalGeometry;
    class PointPhysical;
    
    /*! ~OC~
     Second layer abstract base class (derived from Mapping) for mappings
     that go from a reference element to physical space (hence different point
     types in the two spaces, cf. member function transform).

     In the current design of mappings, these do not store the point
     coordinates, but this can neither be enforced, nor am I sure that this is
     essential. This should be reconsidered at a later stage. The access to
     physical space information is via a reference to the NodeContainer;
     mappings can index into that one with the global node numbers of the
     elements.

     The reinit-function is meant to alert an object of a change of the layout
     of the Element in physical space. Since the reference geometry of the
     Element does not change, but rather only (some of) the vertex positions,
     the Mapping can be adjusted to the new layout.*/

    class MappingReferenceToPhysical : public MappingInterface
    {
    public:
        /// \bug This is a work around for g++ bug 14258 which is fixed in modern compliers so at some point change back
        
        using VectorOfPointsT = const std::vector<PointPhysical>*;

    public:
        MappingReferenceToPhysical()
                : MappingInterface()
        {
        }
        
        // Sets.
        void setNodesPtr(VectorOfPointsT nodes)
        {
            logger.assert(nodes!=nullptr, "Invalid cooridantes passed");
            nodes_ = nodes;
        }
        
        // Methods.
        //! ~OC~ Transform a point from reference space to physical space.
        virtual PointPhysical transform(const PointReference&) const = 0;
        //! ~OC~ Recompute mapping after physical nodes have moved.
        ///\BUG will horribly break everything unless you happen to pass the same  physicalGeometry that you used to construct this mapping
        virtual void reinit(const PhysicalGeometry* const) = 0;
        const PointPhysical& getNodeCoordinates(const std::size_t index) const;

    private:
        ///\TODO fix this properly (for now just made it working)
        const std::vector<PointPhysical>* nodes_; /// Pointer to the global node container.
    };

}
;

#endif /* REFERENCETOPHYSICALM_H_ */
