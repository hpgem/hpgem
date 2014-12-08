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

#ifndef BaseSimplifiedHPP
#define BaseSimplifiedHPP

#include "Base/HpgemUI.hpp"
#include "Integration/ElementIntegrandBase.hpp"
#include "Integration/FaceIntegrandBase.hpp"
#include "Output/TecplotSingleElementWriter.hpp"
#include "Output/VTKTimeDependentWriter.hpp"
#include <functional>

namespace Integration 
{
	class FaceIntegral;
}

namespace Base
{
    class HpgemUISimplified : public HpgemUI,Integration::ElementIntegrandBase<LinearAlgebra::Matrix>,
                                     Integration::FaceIntegrandBase<LinearAlgebra::Matrix>,
                                     Integration::FaceIntegrandBase<LinearAlgebra::NumericalVector>,
                                     Integration::ElementIntegrandBase<LinearAlgebra::NumericalVector>,
                              public Output::TecplotSingleElementWriter
    {
    
    public:
        
        typedef typename MeshManipulator::ConstElementIterator     ConstElementIterator;
        typedef typename MeshManipulator::ElementIterator          ElementIterator;
        typedef typename MeshManipulator::ConstFaceIterator        ConstFaceIterator;
        
        typedef Base::Element                                     ElementT;
        typedef Base::Face                                         FaceT;
        typedef Geometry::PointPhysical                            PointPhysicalT;
        typedef Geometry::PointReference                           PointReferenceT;
        typedef Geometry::PointReference                         PointReferenceOnTheFaceT;
        typedef Integration::FaceIntegral                          FaceIntegralT;
        typedef RectangularMeshDescriptor                          RectangularMeshDescriptorT;
        

        /// You need the basis functions before creating the mesh, because the mesh manipulator
        /// needs to create the elements, and they need some basis functions variables.
    public:

        HpgemUISimplified(unsigned int DIMension, int polynomialOrder=2);

        /// \brief Where the user creates a mesh
        bool virtual initialise()=0;
        
        /// \brief User-defined element integrand
        virtual void elementIntegrand(const ElementT* element, const PointReferenceT& p, LinearAlgebra::Matrix& ret)=0;

        virtual void elementIntegrand(const ElementT* element, const PointReferenceT& p, LinearAlgebra::NumericalVector& ret)=0;
    
        /// \brief User-defined face integrand
        virtual void faceIntegrand(const FaceT* face, const LinearAlgebra::NumericalVector& normal,
                                   const PointReferenceT& p,  LinearAlgebra::Matrix& ret)=0;
        
        virtual void faceIntegrand(const FaceT* face, const LinearAlgebra::NumericalVector& normal,
                				   const PointReferenceT& p,  LinearAlgebra::NumericalVector& ret)=0;

        /// \brief User-defined initial conditions
        virtual double initialConditions(const PointPhysicalT& p)=0;

        /// \brief Does time integration.
        bool solve();
        
        virtual void beforeTimeIntegration(){}
        
        virtual void computeRhsLocal()
        {
#ifdef HPGEM_USE_MPI
            throw "If you want to call the function \'computeRhsLocal\', please implement it";
#endif
        }
        
        virtual void computeRhsFaces()
        {
#ifdef HPGEM_USE_MPI
            throw "If you want to call the function \'computeRhsFaces\', please implement it";
#endif
        }
        
        virtual void interpolate()
        {
#ifdef HPGEM_USE_MPI
            throw "If you want to call the function \'interpolate\', please implement it";
#endif
        }
        
        /**
         * Performs (MPI) synchronisation between meshes
         */
        virtual void synchronize(std::size_t meshID = 0);
        
        ///Preforms all the element integrations
        void doAllElementIntegration(unsigned int meshID=0);
        void doAllFaceIntegration(unsigned int meshID=0);
        
        void registerVTKWriteFunction(std::function<double(Base::Element*, const Geometry::PointReference&, size_t)> function, std::string name)
        {
            VTKDoubleWrite_.push_back({function,name});
        }
        
        void registerVTKWriteFunction(std::function<LinearAlgebra::NumericalVector(Base::Element*, const Geometry::PointReference&, size_t)> function, std::string name)
        {
            VTKVectorWrite_.push_back({function,name});
        }
        
        void registerVTKWriteFunction(std::function<LinearAlgebra::Matrix(Base::Element*, const Geometry::PointReference&, size_t)> function, std::string name)
        {
            VTKMatrixWrite_.push_back({function,name});
        }
        
        virtual void writeToTecplotFile(const ElementT*, const PointReferenceT&, std::ostream&){
            throw "If you want to call the function \'writeToTecplotFile\', please implement it";
        }


    private:
        
        //This is a function that checks the user defined initisation is fine.
        bool  checkInitialisation();
        
        void VTKWrite(Output::VTKTimeDependentWriter& out, double t)
        {
            //you would say this could be done more efficiently, but p.first has different types each time
            for(auto p:VTKDoubleWrite_)
            {
                out.write(p.first,p.second, t);
            }
            for(auto p:VTKVectorWrite_)
            {
                out.write(p.first,p.second, t);
            }
            for(auto p:VTKMatrixWrite_)
            {
                out.write(p.first,p.second, t);
            }
        }
        
    protected:
        
        double endTime_;
        double startTime_;
        double dt_;
        
    private:
        
        std::vector<std::pair<std::function<double(Base::Element*, const Geometry::PointReference&, size_t)>, std::string> > VTKDoubleWrite_;
        std::vector<std::pair<std::function<LinearAlgebra::NumericalVector(Base::Element*, const Geometry::PointReference&, size_t)>, std::string> > VTKVectorWrite_;
        std::vector<std::pair<std::function<LinearAlgebra::Matrix(Base::Element*, const Geometry::PointReference&, size_t)>, std::string> > VTKMatrixWrite_;
        
    };
}

#endif
