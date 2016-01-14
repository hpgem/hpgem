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

#ifndef BaseLinearH
#define BaseLinearH

#include "Base/SerializationInclude.h"
#include "Base/ConfigurationData.h"
#include "Base/Element.h"
#include "Base/Face.h"
#include "Base/HpgemAPISimplified.h"
#include "Base/RectangularMeshDescriptor.h"
#include "Base/TimeIntegration/AllTimeIntegrators.h"
#include "Integration/ElementIntegrandBase.h"
#include "Integration/FaceIntegrandBase.h"
#include "Output/TecplotSingleElementWriter.h"
#include "Output/VTKTimeDependentWriter.h"
#include <functional>

namespace Base
{
    /// \brief Simplified Interface for solving linear PDE's.
    /** At the moment this class is well-suited for problems of the form \f[ l(\partial_t^k u) = f(t) + s(u) \f], where \f$ u\in R^{n_V} \f$ is a vector function, \f$ l:R^{n_V}\rightarrow R^{n_V} \f$ is a linear function, applied on the k-th order time-derivative of \f$ u \f$, and \f$ f(t) + s(u) \f$ is the right-hand side, where \f$ s(u) \f$ is some linear function of \f$ u \f$ (that can depend on arbitrary order spatial derivatives of \f$ u \f$) and \f$ f(t) \f$ is a source term. The resulting set of ODE's will have the form \f[ M\partial_t^ku = Su + f(t)\f], where \f$S\f$ is the stiffness matrix and \f$f(t)\f$ is the source term. If you do not want to store all element- and face matrices for the mass matrix and stiffness matrix, it is advised to use the superclass HpgemAPISimplified instead.
     */
    /** \details Let \f$ \{\phi_{i_B}^e\} \f$ be the set of DG basis functions, where \f$ e \f$ is an element and \f$ i_B \f$ is the index for a basis function corresponding to this element. The basis functions are such that \f$ \phi_{i_B}^e\f$ is non-zero only at element \f$ e \f$. The solution \f$ u \f$ is approximated as follows \f[ u_{i_V}(x,t)|_{x\in e} = \sum_{i_B} \bar{u}^e_{i_V,i_B}(t)\phi_{i_B}(x), \f] for \f$ i_V = 0 .. n_V-1 \f$, where \f$ \bar{u}^e\f$ are the solution coefficients corresponding to element \f$ e \f$.
     
     Let \f$ f \f$ be a face and \f$ i_S \f$ the index of a side of the face (either left or right. Let \f$ (f,i_S) \f$ denote the element at side \f$ i_S \f$ of face \f$ f \f$ (at the boundary we only have a left side). We can write the DG scheme as follows \f[ \sum_{j_V,j_B} M^e_{i_V,i_B;j_V,j_B} \partial_t \bar{u}^e_{j_V,j_B} = f^e_{i_V,i_B}(t) + f^f_{i_V,i_B}(t) + \sum_{j_V,j_B}S^e_{i_V,i_B;j_V.j_B}\bar{u}^e_{j_V.j_B} + \sum_{(f,i_S)=e, j_S,j_V,j_B} S^{f,i_S,j_S}_{i_V,i_B;j_V,j_B}\bar{u}^{(f,j_S)}_{j_V,j_B} \f] where \f$ M^e \f$ is the mass matrix at an element, \f$ S^e \f$ is the stiffness element matrix, \f$ S^f \f$ is the stiffness face matrix, \f$ f^e(t) \f$ is the source term corresponding to an element and \f$ f^f(t) \f$ is the source term corresponding to a boundary face.
     */
    /** \details To solve some linear time depent PDE with this class you should at least do the following:
     * \li Create your own class that inherits this class.
     * \li Implement the function 'createMeshDescription' to create a mesh description (e.g. domain, number of elements, etc.).
     * \li Implement the function 'getInitialConditions' to define the initial condition(s) of your problem.
     * \li Implement the function 'getSourceTerm' to define the source term (e.g. external force) if there is one.
     * \li Implement the function 'getSourceTermAtBoundary' to get the source term (e.g. because of a boundary force / boundary condition) if there is one.
     * \li Implement the function 'integrateSourceTermAtFace' for integrating the source term at a boundary face. One can also choose to implemen the function 'computeIntegrandSourceTermAtFace' for computing the integrand. The integration will be done by an automatic routine.
     * \li Implement the functions 'computeStiffnessMatrixAtElement' and 'computeStiffnessMatrixAtFace' for computing the stiffness matrix at an element or face. One can also choose to implement the functions 'computeIntegrandStiffnessMatrixAtElement' and 'computeIntegrandStiffnessMatrixAtFace' for computing the integrands at the face and element. The integration will be done by an automatic routine.
     */
    /** \details To solve the PDE do the following in the main routine:
     * \li Create an object of your own class, that inherits from this class and has the necessary functions implemented (see list above).
     * \li Call the function 'CreateMesh' to create the mesh.
     * \li Call the function 'setOutputNames' to set the names for the output files.
     * \li Call the function 'solve'.
     */
    /** \details Some other things you can do:
     * \li Implement the function 'getExactSolution' if you know the analytic solution and want to compute the error.
     * \li Implement the function 'integrateInitialSolutionAtElement' for integrating the initial solution at the element (by default this function computes the standard L2 inner product).
     * \li Implement the function 'computeMassMatrixAtElement' if you want to compute the mass matrix (by default a mass matrix is computed based on the L2 inner product).
     * \li Implement the function 'integrateSourceTermAtElement' to compute the source term at an element (if there is a source term, by default the L2 inner product is computed).
     * \li Implement the function 'integrateErrorAtElement' to compute the square of some user-defined norm of the error at an element (by default the L2-norm is computed).
     * \li Override the function 'writeToTecplotFile' to determine what data to write to the output file.
     * \li Override the function 'showProgress' to determine how you want to show the progress of the time integration routine.
     * \li Override the function 'solve' when using another time integration routine than a Runge-Kutta integration method.
     */
    /** \details For an example of using this interface see the application 'ExampleMultipleVariableProblem'.
     */

    template<std::size_t DIM>
    class HpgemAPILinear : public HpgemAPISimplified<DIM>
    {
    public:
        using typename HpgemAPIBase<DIM>::PointPhysicalT;
        using typename HpgemAPIBase<DIM>::PointReferenceT;
        using typename HpgemAPIBase<DIM>::PointReferenceOnFaceT;
        
        // Constructor
        HpgemAPILinear
        (
         const std::size_t numberOfUnknowns,
         const std::size_t polynomialOrder,
         const TimeIntegration::ButcherTableau * const ptrButcherTableau = TimeIntegration::AllTimeIntegrators::Instance().getRule(4, 4),
         const std::size_t numberOfTimeLevels = 0,
         const bool useSourceTerm = false,
         const bool useSourceTermAtBoundary = false
         );
        
        HpgemAPILinear
        (
         const std::size_t numberOfUnknowns,
         const std::size_t polynomialOrder,
         const std::size_t globalNummberOfTimeIntegrationVectors,
         const std::size_t numberOfTimeLevels = 0,
         const bool useSourceTerm = false,
         const bool useSourceTermAtBoundary = false
         );
        
        //If you want to implement the copy constructor and copy assignment operator,
        //make sure the copy constructor of HpgemAPIBase is implemented correctly first.
        HpgemAPILinear(const HpgemAPILinear &other) = delete;
        HpgemAPILinear& operator=(const HpgemAPILinear &other) = delete;
        
        /// \brief Create the mesh.
        virtual void createMesh(const std::size_t numberOfElementsPerDirection, const Base::MeshType meshType) override;
        
        /// \brief Compute the source term at a given physical point.
        virtual LinearAlgebra::MiddleSizeVector getSourceTerm(const PointPhysicalT &pPhys, const double &time, const std::size_t orderTimeDerivative)
        {
            logger(ERROR, "No source term implemented.");
            LinearAlgebra::MiddleSizeVector sourceTerm;
            return sourceTerm;
        }
        
        /// \brief Get the source term at the boundary at a given physical point.
        /// \details The source term at the boundary can be a result of certain boundary conditions (e.g. Neumann boundary conditions).
        virtual LinearAlgebra::MiddleSizeVector getSourceTermAtBoundary(const PointPhysicalT &pPhys, const double &time, const std::size_t orderTimeDerivative)
        {
            
            logger(ERROR, "No source term at the boundary implemented.");
            LinearAlgebra::MiddleSizeVector sourceTerm;
            return sourceTerm;
        }
        
        /// \brief Compute and store the mass matrices.
        virtual void createMassMatrices();
        
        /// \brief Solve the mass matrix equations.
        void solveMassMatrixEquations(const std::size_t timeIntegrationVectorId) override;
        
        /// \brief Compute the integrand for the stiffness matrix.
        virtual LinearAlgebra::MiddleSizeMatrix computeIntegrandStiffnessMatrixAtElement(Base::PhysicalElement<DIM> &element)
        {
            logger(ERROR, "No function for computing the integrand for the stiffness matrix at an element implemented.");
            LinearAlgebra::MiddleSizeMatrix integrandStiffnessMatrix;
            return integrandStiffnessMatrix;
        }
        
        /// \brief Compute the stiffness matrix corresponding to an element.
        virtual LinearAlgebra::MiddleSizeMatrix computeStiffnessMatrixAtElement(Base::Element *ptrElement);
        
        /// \brief Compute the integrand for the stiffness matrix.
        virtual Base::FaceMatrix computeIntegrandStiffnessMatrixAtFace(Base::PhysicalFace<DIM> &face)
        {
            logger(ERROR, "No function for computing the integrand for the stiffness matrix at a face implemented.");
            Base::FaceMatrix integrandStiffnessMatrix;
            return integrandStiffnessMatrix;
        }
        
        /// \brief Compute the stiffness matrix corresponding to a face.
        virtual Base::FaceMatrix computeStiffnessMatrixAtFace(Base::Face *ptrFace);
        
        /// \brief Compute and store stiffness matrices for computing the right hand side.
        virtual void createStiffnessMatrices();
        
        /// \brief Compute the integrand for the source term at the element.
        virtual LinearAlgebra::MiddleSizeVector computeIntegrandSourceTermAtElement(Base::PhysicalElement<DIM> &element, const double time, const std::size_t orderTimeDerivative);
        
        /// \brief Integrate the source term at a single element.
        virtual LinearAlgebra::MiddleSizeVector integrateSourceTermAtElement(Base::Element *ptrElement, const double time, const std::size_t orderTimeDerivative);
        
        /// \brief Compute the integrand for the source term at a face at the boundary.
        virtual LinearAlgebra::MiddleSizeVector computeIntegrandSourceTermAtFace(Base::PhysicalFace<DIM> &face, const double time, const std::size_t orderTimeDerivative)
        {
            logger(ERROR, "No function for computing the integrand for the source term at a face at the domain boundary implemented.");
            LinearAlgebra::MiddleSizeVector integrandSourceTerm;
            return integrandSourceTerm;
        }
        
        /// \brief Integrate the source term at a boundary face.
        virtual LinearAlgebra::MiddleSizeVector integrateSourceTermAtFace(Base::Face *ptrFace, const double time, const std::size_t orderTimeDerivative);
        
        /// \brief Multiply the stiffness matrices with the time integration vector with index 'inputVectorId' and store the result at time integration vector with index 'resultVectorId'.
        virtual void multiplyStiffnessMatrices(const std::size_t inputVectorId, const std::size_t resultVectorId);
        
        /// \brief Multiply the stiffness matrices with the linear combination of time integration vectors with indices 'inputVectorIds', with coefficients given in 'coefficientsInputVectors'. Store the result at time integration vector with index 'resultVectorId'.
        virtual void multiplyStiffnessMatrices(const std::vector<std::size_t> inputVectorIds, const std::vector<double> coefficientsInputVectors, const std::size_t resultVectorId);
        
        /// \brief Add the source term to the time integration vector with index 'resultVectorId'.
        virtual void addSourceTerm(const std::size_t resultVectorId, const double time, const std::size_t orderTimeDerivative);
        
        /// \brief Add the source term at the boundary of the domain to the time integration vector with index 'resultVectorId'.
        virtual void addSourceTermAtBoundary(const std::size_t resultVectorId, const double time, const std::size_t orderTimeDerivative);
        
        /// \brief Compute the right hand side for the DG function with coefficients given by the time integration vector with index 'inputVectorId' and store the result at the time integration vector with index 'resultVectorId'.
        void computeRightHandSide(const std::size_t inputVectorId, const std::size_t resultVectorId, const double time) override;
        
        /// \brief Compute the right hand side for the DG function with coefficients given by a linear combination of time integration vectors with indices 'inputVectorIds' and with coefficients given in 'coefficientsInputVectors'. Store the result at the time integration vector with index 'resultVectorId'.
        void computeRightHandSide(const std::vector<std::size_t> inputVectorId, const std::vector<double> coefficientsInputVectors, const std::size_t resultVectorId, const double time) override;
        
        /// \brief Create and Store things before solving the problem.
        void tasksBeforeSolving() override;
        

        template<class Archive>
        void serialize(Archive & ar, const unsigned int version) 
        {
            ///\todo serialize base classes

            ar & boost::serialization::base_object< Base::HpgemAPISimplified<DIM> >(*this);
        }
        
    protected:
        /// Boolean to indicate if there is a source term.
        const bool useSourceTerm_;
        
        /// Boolean to indicate if there is a source term on the domain boundary (e.g. because of a boundary force or boundary condition).
        const bool useSourceTermAtBoundary_;
        
        /// Index to indicate where the mass matrix is stored
        const std::size_t massMatrixID_;
        
        /// Index to indicate where the stiffness matrix for the elements is stored
        const std::size_t stiffnessElementMatrixID_;
        
        /// Index to indicate where the stiffness matrix for the elements is stored
        const std::size_t stiffnessFaceMatrixID_;
    };
    
    
}

template<class Archive, std::size_t DIM>
inline void save_construct_data(
    Archive & ar, const Base::HpgemAPILinear<DIM> * t, const unsigned int file_version)
{
    // save data required to construct instance
    ar << 100;
    ar << 1;
}

template<class Archive, std::size_t DIM>
inline void load_construct_data(
    Archive & ar, Base::HpgemAPILinear<DIM> * t, const unsigned int file_version)
{
    // retrieve data from archive required to construct new instance
    std::size_t numberOfUnknowns;
    std::size_t polynomialOrder;
    ar >> numberOfUnknowns;
    ar >> polynomialOrder;
    // invoke inplace constructor to initialize instance of my_class
    ::new(t)Base::HpgemAPILinear<DIM>(numberOfUnknowns, polynomialOrder);
}
        
#include "HpgemAPILinear_Impl.h"

#endif


