/* 
 * File:   SavageHutterBase.h
 * Author: irana
 *
 * Created on June 30, 2015, 9:14 AM
 */

#ifndef SAVAGEHUTTERBASE_H
#define	SAVAGEHUTTERBASE_H

#include "RightHandSideComputer.h"
#include "SlopeLimiter.h"
#include "HeightLimiter.h"
#include "Base/HpgemAPISimplified.h"

/// \param[in] numberOfVariables Number of variables in the PDE
/// \param[in] polynomialOrder Polynomial order of the basis functions
/// \param[in] useMatrixStorage Boolean to indicate if element and face matrices for the PDE should be stored
/// \param[in] ptrButcherTableau Pointer to a Butcher Tableau used to do the time integration with a Runge-Kutta scheme. By default this is a RK4 scheme.

struct SHConstructorStruct
{
    std::size_t numOfVariables;
    std::size_t polyOrder;
    std::size_t numElements;
    Base::MeshType meshType;
    Base::ButcherTableau * ptrButcherTableau;
};

class SavageHutterBase : public Base::HpgemAPISimplified<DIM>
{
public:
    SavageHutterBase(const SHConstructorStruct & inputValues);

    virtual ~SavageHutterBase()
    {
        delete rhsComputer_;
        delete slopeLimiter_;
        delete heightLimiter_;
    }

protected:

    /// \brief Create a domain
    Base::RectangularMeshDescriptor<DIM> createMeshDescription(const std::size_t numOfElementPerDirection) override final;
    
    /// Number of variables
    const std::size_t numOfVariables_;

    RightHandSideComputer* rhsComputer_;

    SlopeLimiter* slopeLimiter_;

    HeightLimiter* heightLimiter_;

    ///If the minimum height in an element is below this number, the element is considered to be dry.
    const double dryLimit_;

private:
    virtual SlopeLimiter * createSlopeLimiter(const SHConstructorStruct &inputValues) = 0;
    virtual HeightLimiter * createHeightLimiter(const SHConstructorStruct &inputValues) = 0;
    virtual RightHandSideComputer * createRightHandSideComputer(const SHConstructorStruct &inputValues) = 0;
    
    LinearAlgebra::MiddleSizeVector computeRightHandSideAtElement(Base::Element *ptrElement, 
    LinearAlgebra::MiddleSizeVector &solutionCoefficients, 
    const double time) override final;

    /// \brief Compute the right-hand side corresponding to an internal face
    LinearAlgebra::MiddleSizeVector computeRightHandSideAtFace(Base::Face *ptrFace,
        const Base::Side side,
        LinearAlgebra::MiddleSizeVector &solutionCoefficientsLeft,
        LinearAlgebra::MiddleSizeVector &solutionCoefficientsRight,
        const double time) override final;

    LinearAlgebra::MiddleSizeVector computeRightHandSideAtFace
    (
        Base::Face *ptrFace,
        LinearAlgebra::MiddleSizeVector &solutionCoefficients,
        const double time
        ) override final;

    void computeOneTimeStep(double &time, const double dt) override final;
    void limitSolutionOuterLoop();
    void limitSolutionInnerLoop();

    ///Compute the minimum of the height in the given element
    double getMinimumHeight(const Base::Element *element);

    void tasksBeforeSolving() override final
    {
        //todo: for one face integral you used referenceFaceIntegral (which does not scale with the magnitude of the normal) and for the other you used integrate (which does scale)
        //so it is not clear to me whether or not you need scaling. Please fix as needed
        faceIntegrator_.setTransformation(std::shared_ptr<Base::CoordinateTransformation<DIM> >(new Base::DoNotScaleIntegrands<DIM>(new Base::H1ConformingTransformation<DIM>())));
        Base::HpgemAPISimplified<DIM>::tasksBeforeSolving();
    }
    
    /// \details Make sure timeLevelResult is different from the timeLevelsIn.
    void computeRightHandSide(const std::vector<std::size_t> timeLevelsIn, const std::vector<double> coefficientsTimeLevels, const std::size_t timeLevelResult, const double time) override final;
    
    std::size_t temporaryTimeLevel_;
    
};

#endif	/* SAVAGEHUTTERBASE_H */

