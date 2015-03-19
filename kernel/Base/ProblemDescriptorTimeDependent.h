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

#ifndef ProblemTimeHPP
#define ProblemTimeHPP

#include <functional>

#include "Base/Element.h"
#include "Base/HpgemAPIBase.h"
#include "Base/TimeIntegration/AllTimeIntegrators.h"
#include "Integration/ElementIntegrandBase.h"
#include "Integration/FaceIntegrandBase.h"
#include "Output/TecplotSingleElementWriter.h"
#include "Output/VTKTimeDependentWriter.h"

namespace Integration
{
    class FaceIntegral;
}

namespace Base
{
    ///\deprecated Please build your application in either HpgemAPIBase, HpgemAPISimplified or HpgemAPILinear
    class ProblemDescriptorTimeDependent : public HpgemAPIBase, Integration::ElementIntegrandBase<LinearAlgebra::Matrix>, Integration::FaceIntegrandBase<LinearAlgebra::Matrix>, Integration::FaceIntegrandBase<LinearAlgebra::NumericalVector>, Integration::ElementIntegrandBase<LinearAlgebra::NumericalVector>, public Output::TecplotSingleElementWriter
    {
    public:
        
        using ConstElementIterator = MeshManipulator::ConstElementIterator;
        using ElementIterator = MeshManipulator::ElementIterator;
        using ConstFaceIterator = MeshManipulator::ConstFaceIterator;
        using ElementT = Base::Element;
        using FaceT = Base::Face;
        using PointPhysicalT = Geometry::PointPhysical;
        using PointReferenceT = Geometry::PointReference;
        using FaceIntegralT = Integration::FaceIntegral;

        /// You need the basis functions before creating the mesh, because the mesh manipulator
        /// needs to create the elements, and they need some basis functions variables.
        ProblemDescriptorTimeDependent(std::size_t DIMension, std::size_t polynomialOrder, const Base::ButcherTableau* integrator);

        /// \brief Where the user creates a mesh
        bool virtual initialise() = 0;

        /// \brief User-defined element integrand for the left hand side
        virtual void elementIntegrand(const ElementT *element, const PointReferenceT& p, LinearAlgebra::Matrix& ret) = 0;

        /// \brief User-defined element integrand for the right hand side
        virtual void elementIntegrand(const ElementT *element, const PointReferenceT& p, LinearAlgebra::NumericalVector& ret) = 0;

        /// \brief User-defined face integrand for the left hand side
        virtual void faceIntegrand(const FaceT* face, const LinearAlgebra::NumericalVector& normal, const PointReferenceT& p, LinearAlgebra::Matrix& ret) = 0;

        /// \brief User-defined face integrand for the right hand side
        virtual void faceIntegrand(const FaceT *face, const LinearAlgebra::NumericalVector& normal, const PointReferenceT& p, LinearAlgebra::NumericalVector& ret) = 0;

        /// \brief User-defined initial conditions
        virtual double initialConditions(const PointPhysicalT& p) = 0;

        virtual void interpolate();

        /// \brief Everything that must be done before starting the time-integration.
        ///
        /// For example, user-defined integrals can be computed here. 
        virtual void beforeTimeIntegration()
        {
        }
        
        virtual void computeRhsLocal()
        {
            throw "If you want to call the function \'computeRhsLocal\', please implement it";
        }
        
        virtual void computeRhsFaces()
        {
            throw "If you want to call the function \'computeRhsFaces\', please implement it";
        }
        
        /**
         * \brief Executes one time step.
         */
        virtual void computeOneTimeStep();

        /**
         * \brief Performs (MPI) synchronisation between meshes
         */
        virtual void synchronize(std::size_t meshID = 0);

        /// \brief Does time integration.
        virtual bool solve();

        /// \brief Performs all the element integrations
        void doAllElementIntegration(std::size_t meshID = 0);

        /// \brief Performs all the face integrations
        void doAllFaceIntegration(std::size_t meshID = 0);

        void registerVTKWriteFunction(std::function<double(Base::Element*, const Geometry::PointReference&, std::size_t)> function, std::string name)
        {
            VTKDoubleWrite_.push_back( {function, name});
        }
        
        void registerVTKWriteFunction(std::function<LinearAlgebra::NumericalVector(Base::Element*, const Geometry::PointReference&, std::size_t)> function, std::string name)
        {
            VTKVectorWrite_.push_back( {function, name});
        }
        
        void registerVTKWriteFunction(std::function<LinearAlgebra::Matrix(Base::Element*, const Geometry::PointReference&, std::size_t)> function, std::string name)
        {
            VTKMatrixWrite_.push_back( {function, name});
        }
        
        virtual void writeToTecplotFile(const ElementT*, const PointReferenceT&, std::ostream&)
        {
            throw "If you want to call the function \'writeToTecplotFile\', please implement it";
        }
        
    protected:
        
        double endTime_;
        double startTime_;
        double dt_;

    private:
        
        ///Function that checks the user defined initialisation is fine.
        bool checkInitialisation();

        void VTKWrite(Output::VTKTimeDependentWriter& out, double t)
        {
            //you would say this could be done more efficiently, but p.first has different types each time
            for (auto p : VTKDoubleWrite_)
            {
                out.write(p.first, p.second, t);
            }
            for (auto p : VTKVectorWrite_)
            {
                out.write(p.first, p.second, t);
            }
            for (auto p : VTKMatrixWrite_)
            {
                out.write(p.first, p.second, t);
            }
        }
        
        std::vector<std::pair<std::function<double(Base::Element*, const Geometry::PointReference&, std::size_t)>, std::string> > VTKDoubleWrite_;
        std::vector<std::pair<std::function<LinearAlgebra::NumericalVector(Base::Element*, const Geometry::PointReference&, std::size_t)>, std::string> > VTKVectorWrite_;
        std::vector<std::pair<std::function<LinearAlgebra::Matrix(Base::Element*, const Geometry::PointReference&, std::size_t)>, std::string> > VTKMatrixWrite_;
        const ButcherTableau* integrator_;
        
    };
}

#endif
