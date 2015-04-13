/*
 * Euler.h
 *
 *  Created on: Mar 20, 2015
 *      Author: marnix
 */

#ifndef EULER_H_
#define EULER_H_

#include "Base/HpgemAPISimplified.h"

class Euler : public Base::HpgemAPISimplified
{
public:
	Euler
	(
	 const std::size_t dimension,
	 const std::size_t numOfVariables,
	 const std::size_t polynomialOrder,
	 const Base::ButcherTableau * const ptrButcherTableau
	);

    /// \brief Create a domain
    Base::RectangularMeshDescriptor createMeshDescription(const std::size_t numOfElementPerDirection) override;

    /// \brief Compute the initial solution at a given point in space and time.
    LinearAlgebra::NumericalVector getInitialSolution(const PointPhysicalT &pPhys, const double &startTime, const std::size_t orderTimeDerivative = 0) override;

    /// *****************************************
    /// ***   Element integration functions   ***
    /// *****************************************

    /// Compute solution at an element
    LinearAlgebra::NumericalVector computeSolutionAtElement(const Base::Element *ptrElement, const LinearAlgebra::NumericalVector &solutionCoefficients, const Geometry::PointReference &pRef);

    /// Compute integrand of righthandside on an element
    LinearAlgebra::NumericalVector integrandRightHandSideOnRefElement(const Base::Element *ptrElement, const double &time, const Geometry::PointReference &pRef, const LinearAlgebra::NumericalVector &solutionCoefficients);

    /// \brief Compute the right hand side on an element
    LinearAlgebra::NumericalVector computeRightHandSideAtElement(Base::Element *ptrElement,	LinearAlgebra::NumericalVector &solutionCoefficients, const double time) override;

    /// *****************************************
    /// ***    face integration functions     ***
    /// *****************************************

    /// \brief Compute the Roe Riemann Flux.
    LinearAlgebra::NumericalVector RoeRiemannFluxFunction(const LinearAlgebra::NumericalVector &qReconstructionLeft, const LinearAlgebra::NumericalVector &qReconstructionRight, LinearAlgebra::NumericalVector &normal, const Base::Side &iSide);

    /// \brief Compute the integrand for the right hand side for the reference face corresponding to a boundary face.
    LinearAlgebra::NumericalVector integrandRightHandSideOnRefFace(const Base::Face *ptrFace, const double &time, const Geometry::PointReference &pRef, const LinearAlgebra::NumericalVector &solutionCoefficients);

    /// \brief Compute the integrand for the right hand side for the reference face corresponding to an internal face.
    LinearAlgebra::NumericalVector integrandRightHandSideOnRefFace(const Base::Face *ptrFace, const double &time, const Geometry::PointReference &pRef, const Base::Side &iSide, const LinearAlgebra::NumericalVector &solutionCoefficientsLeft, const LinearAlgebra::NumericalVector &solutionCoefficientsRight);

    /// \brief Compute the right-hand side corresponding to a boundary face
    LinearAlgebra::NumericalVector computeRightHandSideAtFace(Base::Face *ptrFace, LinearAlgebra::NumericalVector &solutionCoefficients, const double time) override final;

    /// \brief Compute the right-hand side corresponding to an internal face
    LinearAlgebra::NumericalVector computeRightHandSideAtFace(Base::Face *ptrFace, const Base::Side side, LinearAlgebra::NumericalVector &solutionCoefficientsLeft, LinearAlgebra::NumericalVector &solutionCoefficientsRight, const double time) override final;


    /// *****************************************
    /// ***    		Various Functions         ***
    /// *****************************************

    void showProgress(const double time, const std::size_t timeStepID);

private:
    /// Dimension of the domain
    const std::size_t DIM_;

    /// Specific heat ratio
    const double gamma_ = 1.4;

	/// Number of variables
	const std::size_t numOfVariables_;
};

#endif /* EULER_H_ */
