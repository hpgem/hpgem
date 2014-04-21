#ifndef NumericalVectorHPP
#define NumericalVectorHPP

#include "GlobalNamespaceLinearAlgebra.hpp"
#include <valarray>
#include <iostream>


#define IAMNICO

namespace LinearAlgebra
{
    
    //This is dervied from valarray so import that information
    using std::valarray;
    using std::ostream;
    
    /// \class NumericalVector
    /// \brief This is a vector of doubles
    ///
    /// \details
    /// This implements a vector of doubles and all the standard operators for it.
    /// Note it is encapulating a valarray for its data storage.
    class NumericalVector
    {
        
    public:
        
        NumericalVector();
        
        NumericalVector(int m);
        
        NumericalVector(const NumericalVector& other);
        
        NumericalVector(const double array[], int size);
        
        void resize(unsigned int size);

        NumericalVector& operator= (const NumericalVector& right);
        
        NumericalVector operator+ (const NumericalVector& right);
        
        NumericalVector operator+ (const NumericalVector& right) const;
        
        NumericalVector operator- (const NumericalVector& right);
        
        NumericalVector operator- (const NumericalVector& right) const;
        
        NumericalVector operator* (const double& right);
        
        NumericalVector operator* (const double& right) const;
        
        double operator* (const NumericalVector& right) const;

        NumericalVector& operator/= (const double& right);
        
        NumericalVector operator/ (const double& right);
        
        void axpy(double a, const NumericalVector& x);
        
        bool operator== (const NumericalVector& right);

        bool operator== (const NumericalVector& right) const;
        
        bool operator< (const NumericalVector& right) const;

        NumericalVector& operator+= (const NumericalVector& right);
        
        NumericalVector& operator-= (const NumericalVector& right);
        
        NumericalVector& operator*= (const double& right);
        
        double& operator[] (const unsigned int n);
        
        const double&  operator[] (const unsigned int n) const {return data_[n];}
        
        double& operator() (const unsigned int n) {return data_[n];}
        
        const double&  operator() (const unsigned int n) const {return data_[n];}
        
        int size() const {return data_.size();}
        
        int size() {return data_.size();}
        
        friend NumericalVector operator*(const double& left, const NumericalVector& right);
        
        friend NumericalVector   operator-(const NumericalVector& right);
 
        friend ostream& operator<<(ostream& os, const NumericalVector& A);
        
   
    private:

        valarray<double> data_;
        
    };
    
}
#endif
