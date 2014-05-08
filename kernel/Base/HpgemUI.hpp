/*
 This file forms part of hpGEM. This package has been developed over a number of years by various people at the University of Twente and a full list of contributors can be found at
 http://hpgem.org/about-the-code/team
 
 This code is distributed using BSD 3-Clause License. A copy of which can found below.
 
 
 Copyright (c) 2014, Univesity of Twenete
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 
 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef BASE_HPP
#define BASE_HPP

#include "Base/MeshMoverBase.hpp"
#include "Base/MeshManipulator.hpp"
#include "Base/GlobalNamespaceBase.hpp"
#include "Base/RectangularMeshDescriptor.hpp"
#include "Output/TecplotDiscontinuousSolutionWriter.hpp"
#include "Integration/ElementIntegral.hpp"
#include "Integration/FaceIntegral.hpp"

#include "vector"

namespace Base
{
    class HpgemUI
    {
    public:
        typedef typename MeshManipulator::ConstElementIterator     ConstElementIterator;
        typedef typename MeshManipulator::ElementIterator          ElementIterator;
        typedef typename MeshManipulator::ConstFaceIterator        ConstFaceIterator;
        typedef typename MeshManipulator::FaceIterator             FaceIterator;
        
    public:
        typedef Base::Element                                      ElementT;
        typedef Base::Face                                         FaceT;
        typedef RectangularMeshDescriptor                          RectangularMeshDescriptorT;
        typedef MeshManipulator                                    MeshManipulatorT;
        typedef MeshMoverBase                                      MeshMoverBaseT;
        typedef Geometry::PointPhysical                            PointPhysicalT;
        typedef Geometry::PointReference                           PointReferenceT;

        typedef unsigned int                                            MeshId;
        typedef std::vector<unsigned int>                               VectorOfUIntegers;
       
        typedef std::vector<MeshManipulatorT* >                         VectorOfMeshManipulatorT;
        typedef std::string                                             String;
                
        //typedef BasisFunctions<DIM>     BasisFunctionT;
        

    public:
             
        HpgemUI(GlobalData* const global, const ConfigurationData* config);

        virtual ~HpgemUI() ;
        
        

        /// \brief Gives the pointer of meshMoverBase class to mesh.
        virtual bool initialiseMeshMover(const MeshMoverBaseT* meshMoverBase, unsigned int meshID);

            /// Creating a mesh with in-house remesher.
        MeshId addMesh(const RectangularMeshDescriptorT& meshDescriptor, const MeshType& meshType = RECTANGULAR, int nrOfElementMatrixes=0, int nrOfElementVectors=0,int nrOfFaceMatrixes=0, int nrOfFaceVectors=0);
            /// Reading a mesh from a file, currently only Centaur is supported.
        MeshId addMesh(const String& fileName){}
        
        unsigned int getNumberOfElements(MeshId id)const {return meshes_[id]->getNumberOfElements();}
        
        ConstElementIterator    elementColBegin(MeshId mId=0)const;
        ConstElementIterator    elementColEnd(MeshId mId=0)const;
        
        ElementIterator         elementColBegin(MeshId mId=0);
        ElementIterator         elementColEnd(MeshId mId=0);
        
        ConstFaceIterator       faceColBegin(MeshId mId=0)const;
        ConstFaceIterator       faceColEnd(MeshId mId=0)const;
        
        FaceIterator            faceColBegin(MeshId mId=0);
        FaceIterator            faceColEnd(MeshId mId=0);
        

        /// \brief Virtual function that should be overwritten by specific problem, specifies initial conditions.
        //virtual void initialCondition() const;

    protected:
        VectorOfMeshManipulatorT                meshes_;
       
        GlobalData* const                       globalData_;
        const ConfigurationData* const          configData_;
    };
};
#endif
