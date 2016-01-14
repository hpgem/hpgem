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

#include <cmath>
#include <functional>
#include <chrono>

#include "Base/SerializationInclude.h"
#include "Base/CommandLineOptions.h"
#include "Base/ConfigurationData.h"
#include "Base/Element.h"
#include "Base/Face.h"
#include "Base/HpgemAPILinear.h"
#include "Base/RectangularMeshDescriptor.h"
#include "Integration/ElementIntegral.h"
#include "Integration/FaceIntegral.h"
#include "Output/TecplotDiscontinuousSolutionWriter.h"
#include "Output/TecplotSingleElementWriter.h"
#include "Utilities/BasisFunctions2DH1ConformingTriangle.h"
#include "Logger.h"

const static std::size_t DIM = 2;

/// Linear advection equation du/dt + a[0] du/dx + a[1] du/dy = 0.
/// This class is meant for testing purposes.
//  Please verify that nobody went and tested a broken feature
//  Please keep the problem modelled here reasonably close to the linear advection problem
class AdvectionTest : public Base::HpgemAPILinear<DIM>
{
public:
    ///Constructor. Assign all private variables.
    AdvectionTest(std::size_t p) :
    HpgemAPILinear<DIM>(1, p)
    {
        //Choose the "direction" of the advection.
        //This cannot be implemented with iterators, and since the dimension is
        //not always 2, this is the most generic way to write it.
        for (std::size_t i = 0; i < DIM; ++i)
        {
            a[i] = .1 + 0.0 * i;
        }
    }
    
    /// Create a mesh description
    Base::RectangularMeshDescriptor<DIM> createMeshDescription(const std::size_t numberOfElementsPerDirection) override final
    {
        //describes a rectangular domain
        Base::RectangularMeshDescriptor<DIM> description;
        
        //this demo will use the square [0,1]^2
        for (std::size_t i = 0; i < DIM; ++i)
        {
            description.bottomLeft_[i] = 0;
            description.topRight_[i] = 1;
            //Define elements in each direction.
            description.numberOfElementsInDIM_[i] = numberOfElementsPerDirection;
            
            //Choose whether you want periodic boundary conditions or other (solid wall)
            //boundary conditions.
            description.boundaryConditions_[i] = Base::BoundaryType::PERIODIC;
        }
        
        return description;
    }
    
    ///Compute phi_i*(a.grad(phi_j)) on a reference point on an element for all
    ///basisfunctions phi_i and phi_j.
    ///You pass the reference point to the basisfunctions. Internally the basisfunctions will be mapped to the physical element
    ///so you wont have to do any transformations yourself
    LinearAlgebra::MiddleSizeMatrix computeIntegrandStiffnessMatrixAtElement(Base::PhysicalElement<DIM>& element) override final
    {
        std::size_t numberOfBasisFunctions = element.getElement()->getNumberOfBasisFunctions();
        LinearAlgebra::MiddleSizeMatrix&  result = element.getResultMatrix();
        for (std::size_t i = 0; i < numberOfBasisFunctions; ++i)
        {
            for (std::size_t j = 0; j < numberOfBasisFunctions; ++j)
            {
                result(j, i) = element.basisFunction(i) * (a * element.basisFunctionDeriv(j));
            }
        }
        
        return result;
    }
    
    /// \brief Compute the integrals of the left-hand side associated with faces.
    ///
    ///For every internal face, we want to compute the integral of the flux
    ///for all basisfunctions phi_i and phi_j that are non-zero on that face.
    ///For boundary faces, similar expressions can be obtained depending of the type of boundary condition.
    ///This function will compute these integrands for all basisfunctions phi_i and phi_j
    ///on a certain face at a reference point p. Then the integral can later be computed with appropriate (Gauss-)quadrature rules.
    ///The resulting matrix of values is then given in the matrix integrandVal, to which we passed a reference when calling it.
    ///Please note that you pass a reference point to the basisfunctions and the
    ///transformations are done internally. The class FaceMatrix consists of four element matrices for internal faces and one element matrix for faces on the boundary. Each element matrix corresponds to a pair of two adjacent elements of the face.
    Base::FaceMatrix computeIntegrandStiffnessMatrixAtFace(Base::PhysicalFace<DIM>& face) override final
    {
        //Get the number of basis functions, first of both sides of the face and
        //then only the basis functions associated with the left and right element.
        std::size_t numberOfBasisFunctions = face.getFace()->getNumberOfBasisFunctions();
        
        //Resize the result to the correct size and set all elements to 0.
        Base::FaceMatrix& integrandVal = face.getResultMatrix();
        integrandVal *= 0;
        
        //Check if the normal is in the same direction as the advection.
        //Note that normal does not have length 1!
        const double A = a * face.getUnitNormalVector();
        
        //Compute all entries of the integrand at this point:
        for (std::size_t i = 0; i < numberOfBasisFunctions; ++i)
        {
            Base::Side sideBasisFunction = face.getFace()->getSide(i);
            for (std::size_t j = 0; j < numberOfBasisFunctions; ++j)
            {
                //Give the terms of the upwind flux.
                //Advection in the same direction as outward normal of the left element:
                if ((A > 1e-12) && (sideBasisFunction == Base::Side::LEFT))
                {
                    integrandVal(j, i) = -(a * face.basisFunctionUnitNormal(j)) * face.basisFunction(i);
                }
                //Advection in the same direction as outward normal of right element:
                else if ((A < -1e-12) && (sideBasisFunction == Base::Side::RIGHT))
                {
                    integrandVal(j, i) = -(a * face.basisFunctionUnitNormal(j)) * face.basisFunction(i);
                }
                //Advection orthogonal to normal:
                else if (std::abs(A) < 1e-12)
                {
                    integrandVal(j, i) = -(a * face.basisFunctionUnitNormal(j)) * face.basisFunction(i) / 2.0;
                }
            }
        }
        
        return integrandVal;
    }
    
    /// Define a solution at time zero.
    double getSolutionAtTimeZero(const PointPhysicalT& point)
    {
        return (std::sin(2 * M_PI * point[0]) * std::sin(2 * M_PI * point[1]));
    }
    
    /// Define the exact solution. In this case that is \f$ u_0(\vec{x}-\vec{a}t) \f$, where \f$ u_0 \f$ is the solution at time zero.
    LinearAlgebra::MiddleSizeVector getExactSolution(const PointPhysicalT& point, const double &time, const std::size_t orderTimeDerivative) override final
    {
        LinearAlgebra::MiddleSizeVector exactSolution(1);
        if(orderTimeDerivative == 0)
        {
            PointPhysicalT displacement(-a*time);
            exactSolution(0) = getSolutionAtTimeZero(point + displacement);
            return exactSolution;
        }
        else
        {
            logger(ERROR, "No exact solution for order time derivative % implemented", orderTimeDerivative);
            exactSolution(0) = 0;
            return exactSolution;
        }
    }
    
    /// Define the initial conditions. In this case it is just the exact solution at the start time.
    LinearAlgebra::MiddleSizeVector getInitialSolution(const PointPhysicalT& point, const double &startTime, const std::size_t orderTimeDerivative) override final
    {
        return getExactSolution(point, startTime, orderTimeDerivative);
    }
    
    
        
    // When the class Archive corresponds to an output archive, the
    // & operator is defined similar to <<.  Likewise, when the class Archive
    // is a type of input archive the & operator is defined similar to >>.
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & boost::serialization::base_object<  Base::HpgemAPILinear<DIM> >(*this);
        ar & a;
    }
private:
    
    ///Advective vector
    LinearAlgebra::SmallVector<DIM> a;
};

auto& n = Base::register_argument<std::size_t>('n', "numelems", "Number of Elements", true);
auto& p = Base::register_argument<std::size_t>('p', "poly", "Polynomial order", true);

template<class Archive>
inline void save_construct_data(
    Archive & ar, const AdvectionTest * test, const unsigned int file_version)
{
    // save data required to construct instance (for now cheat and acquire the info from somewhere else)
    ar << p.getValue();
}

template<class Archive>
inline void load_construct_data(
    Archive & ar, AdvectionTest * t, const unsigned int file_version)
{
    // retrieve data from archive required to construct new instance
    std::size_t polynomialOrder;
    ar >> polynomialOrder;
    // invoke inplace constructor to initialize instance of my_class
    ::new(t)AdvectionTest(polynomialOrder);
}

///Make the problem and solve it.
int main(int argc, char **argv)
{
    Base::parse_options(argc, argv);
    const Base::MeshType meshType = Base::MeshType::RECTANGULAR;
    std::vector<std::string> variableNames;
    variableNames.push_back("u");
    AdvectionTest test(p.getValue());
    test.createMesh(n.getValue(), meshType);
    test.setOutputNames("output", "AdvectionTest", "AdvectionTest", variableNames);

    // save data to archive
    {
        std::ofstream outputFileStream("advectionFile");
        boost::archive::text_oarchive outputArchive(outputFileStream);
        //save a pointer so we don't need to figure out how to call the constructor when reloading
        outputArchive << &test;
    	// archive and stream closed when destructors are called
    }
    auto startTime = std::chrono::steady_clock::now();
    test.solve(Base::startTime.getValue(), Base::endTime.getValue(), Base::dt.getValue(), Base::numberOfSnapshots.getValue(), true);
    auto endTime = std::chrono::steady_clock::now();
    
    AdvectionTest* test2;
    {
        // create and open an archive for input
        std::ifstream inputFileStream("advectionFile");
        boost::archive::text_iarchive inputArchive(inputFileStream);
        // read class state from archive
        inputArchive >> test2;
        // archive and stream closed when destructors are called
    }
    
    test2->createMesh(n.getValue(), meshType);
    
    startTime = std::chrono::steady_clock::now();
    test2->solve(Base::startTime.getValue(), Base::endTime.getValue(), Base::dt.getValue(), Base::numberOfSnapshots.getValue(), true);
    endTime = std::chrono::steady_clock::now();

    logger(INFO, "Simulation took %ms.", std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count());

    //boost created some space for us, but it won't and shouldn't clean up
    delete test2;

    return 0;
}
