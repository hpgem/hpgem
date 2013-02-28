#ifndef Expansion_hpp
#define Expansion_hpp

#include "../LinearAlgebra/NumericalVector.hpp"
#include "../Geometry/PointReference.hpp"
#include "../Base/BasisFunctionSet.hpp"
#include "../Base/Element_FAKE.hpp"
#include "../LinearAlgebra/Matrix.hpp"

namespace Base
{
    template <unsigned int dim>
    class Expansion
    {

    public:

        typedef BasisFunctionSet<dim>::PointReferenceT PointReferenceT;
        Expansion(const BasisFunctionSet<1>* BFSetPtr) :
            BFSetPtr_(BFSetPtr),
            coeff_(BFSetPtr->size()),
            size_(BFSetPtr->size())
        {}

        ~Expansion() {};

        void setBasisFunction(const BasisFunctionSet<1>* BFSetPtr)
        {
            BFSetPtr_ = BFSetPtr;
            size_ = BFSetPtr_->size();
            coeff_.resize(size_);
        }

        double EvalBasisFunction(unsigned int iBF, const PointReferenceT& p) const
        {
            return BFSetPtr_->Eval(iBF,p);
        }

        double EvalDerivBasisFunction(unsigned int iBF, const unsigned int jDir, const PointReferenceT& p) const
        {
            return BFSetPtr_->EvalDeriv(iBF,jDir,p);
        }

        void PhysGradientOfBasisFunction(
                const Element<dim>& el,
                const unsigned int iBF,
                const PointReferenceT& p,
                LinearAlgebra::NumericalVector& ret) const
        {
            ret.resize(dim);
            // get derivatives df
            // get Jacobian    jac
            // solve the linear system and return
        }

	double Eval(const PointReferenceT& p) const
	{
	  double ret(0.);
	  for (unsigned int i=0; i<size_; ++i)
	    ret += (coeff_[i] * BFSetPtr_->Eval(i,p));

	  return ret;
	}

	double EvalDeriv(const unsigned int jDir, const PointReferenceT& p) const
	{
	  double ret(0.);
	  for (unsigned int i=0; i<size_; ++i)
	    ret += (coeff_[i] * BFSetPtr_->EvalDeriv(i,jDir,p));

	  return ret;
	}

      private:
	BasisFunctionSet<1>* BFSetPtr_;
	LinearAlgebra::NumericalVector coeff_;
	unsigned int size_;
    };

};

#endif // Expansion_hpp
