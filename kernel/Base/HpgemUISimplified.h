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

#include "Base/HpgemAPIBase.h"
#include "Integration/ElementIntegrandBase.h"
#include "Integration/FaceIntegrandBase.h"
#include "Output/TecplotSingleElementWriter.h"
#include "Output/VTKTimeDependentWriter.h"
#include "Element.h"
#include <functional>

namespace Integration
{
    class FaceIntegral;
}

namespace Base
{
    /// \brief Interface for solving hyperbolic problems.
    /** \details To solve some linear time depent PDE you should do the following:
     * \li Create your own class that inherits this class.
     * \li Implement the function 'initialise' for creating the mesh.
     * \li Implement the function 'initialConditions' to define the initial conditions of your problem.
     * \li Implement the functions 'elementIntegrand' and 'faceIntegrand' for defining the integrands for element matrices and vectors and face matrices and vectors.
     * \li Implement the functions 'computeRhsLocal' and 'computeRhsFaces' for computing the right-hand-side corresponding to your time integration method.
     * \li Implement the function 'beforeTimeIntegration' when multiple element/face matrices/vectors are required.
     * \li Override the function 'solve' when using another time integration routine than forward Euler.
     * \li Implement the function 'writeToTecplotFile' to determine what data to write to the output file.
     */
    /** \details To solve the PDE do the following in the main routine:
     * \li Create an object of your own class. 
     * \li Define how the solution should be written in the VTK files using the function 'registerVTKWriteFunction'.
     * \li Call the function 'solve'.
     */
    /** \details For an example of using this interface see the application 'TutorialAdvection'.
     * \deprecated (use HpgemAPISimplified instead)
     * \todo clean up dependencies in the applications folder, then remove this class
     */
    class HpgemUISimplified : public HpgemAPIBase, Integration::ElementIntegrandBase<LinearAlgebra::Matrix>, Integration::FaceIntegrandBase<LinearAlgebra::Matrix>, Integration::FaceIntegrandBase<LinearAlgebra::NumericalVector>, Integration::ElementIntegrandBase<LinearAlgebra::NumericalVector>, public Output::TecplotSingleElementWriter
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
        HpgemUISimplified(std::size_t DIMension, std::size_t polynomialOrder = 2, std::size_t nrOfUnknowns = 1);
        
        HpgemUISimplified(const HpgemUISimplified &other) = delete;

        /// \brief Where the user creates a mesh
        bool virtual initialise() = 0;

        /// \brief User-defined initial conditions
        /// \details This function returns the initial condition at the physical point p.
        virtual double initialConditions(const PointPhysicalT& p) = 0;

        /// \brief Solve the problem.
        bool solve();

        /// \brief Everything that must be done before starting the time-integration.
        /// \details For example, user-defined integrals can be computed here, like stiffness matrices.
        virtual void beforeTimeIntegration()
        {
        }
        
        /// \brief Compute the right hand side corresponding to the element integrals.
        /// \details The right-hand side in this case corresponds to all integral terms that do not include the time derivative of the solution. For example, the right-hand side could be of the form \f[ Su+f \f], where \f$ S \f$ is the stiffness matrix, \f$ u \f$ are the coefficients of the basis functions of the solution and \f$ f \f$ is some external force.
        virtual void computeRhsLocal()
        {
#ifdef HPGEM_USE_MPI
            throw "If you want to call the function \'computeRhsLocal\', please implement it";
#endif
        }
        
        /// \brief Compute the right hand side corresponding to the face integrals.
        /// \details The right-hand side in this case corresponds to all integral terms that do not include the time derivative of the solution. For example, the right-hand side could be of the form \f$ Su+f \f$, where \f$ S \f$ is the stiffness matrix, \f$ u \f$ are the coefficients of the basis functions of the solution and \f$ f \f$ is some external force.
        virtual void computeRhsFaces()
        {
#ifdef HPGEM_USE_MPI
            throw "If you want to call the function \'computeRhsFaces\', please implement it";
#endif
        }
        
        /// \brief Interpolates the new solution from the right hand side computed before.
        virtual void interpolate();

        /**
         * \brief Performs (MPI) synchronisation between meshes
         */
        virtual void synchronize(std::size_t meshID = 0);

        /// \brief Performs all the element integrations.
        void doAllElementIntegration(std::size_t meshID = 0);

        /// \brief Performs all the face integrations
        void doAllFaceIntegration(std::size_t meshID = 0);

        /// \brief Define how the solution should be written in the VTK files.
        /// \details For an example of using this function, see for example the application 'TutorialAdvection' to find out how to use this function.
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
        
        /// \brief Write data to some output file that can be read by TecPlot.
        /// \details See some for example the application 'TutorialAdvection' to find out how to implement this function.
        virtual void writeToTecplotFile(const ElementT*, const PointReferenceT&, std::ostream&)
        {
            logger(ERROR, "If you want to call the function \'writeToTecplotFile\', please implement it");
        }
        
    protected:
        
        double endTime_;
        double startTime_;
        double dt_;

    private:
        
        //allow multiple templated functions with the same arguments, but different return types
        template<typename T>
        /// \brief User-defined element integrand for the right hand side
        T elementIntegrand(const ElementT* element, const PointReferenceT& p) = delete;
        
        template<typename T>
        /// \brief User-defined face integrand for the right hand side
        T faceIntegrand(const FaceT* face, const LinearAlgebra::NumericalVector& normal, const PointReferenceT& p) = delete;
        
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
        
    };
}

#endif
