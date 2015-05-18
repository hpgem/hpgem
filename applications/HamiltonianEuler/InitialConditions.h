#ifndef INITIAL_CONDITIONS_HH
#define INITIAL_CONDITIONS_HH

#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <vector>
#include <sstream>
#include <cstring>
#include <cmath>
#include <complex>
using namespace std;

#include "Geometry/PointPhysical.h"
#include "Geometry/PointReference.h"
#include "Base/Element.h"
#include  "LinearAlgebra/MiddleSizeVector.h"
#include "Integration/ElementIntegrandBase.h"

class ExactSolutionBase
{
public:
    using PointPhysicalT = Geometry::PointPhysical;
public:
    ExactSolutionBase()
    {
    }
    ~ExactSolutionBase()
    {
    }
    
    virtual double getU(const PointPhysicalT& pPhys, double t) const=0;
    virtual double getV(const PointPhysicalT& pPhys, double t) const=0;
    virtual double getW(const PointPhysicalT& pPhys, double t) const=0;
    virtual double getP(const PointPhysicalT& pPhys, double t) const=0;

public:
    static const double Pi;
};

// class ExternalSolution:public ExactSolutionBase
// {
// 	public:
//     	Compressible3DOneThirdPeriodic(const std::string& file)
//     		:baseFileName_(file)
//     		{}
// 
// 		double getU(const PointPhysicalT& pPhys, double t)const;
// 		double getV(const PointPhysicalT& pPhys, double t)const;
// 		double getW(const PointPhysicalT& pPhys, double t)const;
// 		double getP(const PointPhysicalT& pPhys, double t)const;
// 	private:
// 		std::string 	baseFileName_;
// 		
// };

class Compressible3DOneThirdPeriodic : public ExactSolutionBase
{
public:
    using ComplexNumber = complex<double>;

public:
    Compressible3DOneThirdPeriodic(double lx, double ly, double lz);

    double getU(const PointPhysicalT& pPhys, double t) const;
    double getV(const PointPhysicalT& pPhys, double t) const;
    double getW(const PointPhysicalT& pPhys, double t) const;
    double getP(const PointPhysicalT& pPhys, double t) const;

private:
    const ComplexNumber iComplex_;
    double lx_;
    double ly_;
    double lz_;
    double k_;
    double l_;
    double sigma_;
    
};

class Compressible3DPeriodic : public ExactSolutionBase
{
    
public:
    Compressible3DPeriodic();

    double getU(const PointPhysicalT& pPhys, double t) const;
    double getV(const PointPhysicalT& pPhys, double t) const;
    double getW(const PointPhysicalT& pPhys, double t) const;
    double getP(const PointPhysicalT& pPhys, double t) const;
};

class Compressible3DWalls : public ExactSolutionBase
{
    
public:
    Compressible3DWalls();

    double getU(const PointPhysicalT& pPhys, double t) const;
    double getV(const PointPhysicalT& pPhys, double t) const;
    double getW(const PointPhysicalT& pPhys, double t) const;
    double getP(const PointPhysicalT& pPhys, double t) const;
};

class Incompressible3DPeriodic : public ExactSolutionBase
{
    
public:
    Incompressible3DPeriodic();

    double getU(const PointPhysicalT& pPhys, double t) const;
    double getV(const PointPhysicalT& pPhys, double t) const;
    double getW(const PointPhysicalT& pPhys, double t) const;
    double getP(const PointPhysicalT& pPhys, double t) const;
};

class InitialVelocityConstructorTaylor : public ExactSolutionBase
{
public:
    using ComplexNumber = complex<double>;

    using Coefficients = std::vector<ComplexNumber>;
    using Vector = std::vector<double>;
    using IstreamIterator = istream_iterator<double>;

    using ExactSolutionBase::PointPhysicalT;
    using ExactSolutionBase::Pi;

public:
    InitialVelocityConstructorTaylor();
    InitialVelocityConstructorTaylor(double lx, double ly, double lz);

    virtual double getU(const PointPhysicalT& pPhys, double t) const;
    virtual double getV(const PointPhysicalT& pPhys, double t) const;
    virtual double getW(const PointPhysicalT& pPhys, double t) const;
    virtual double getP(const PointPhysicalT& pPhys, double t) const;

private:
    double lx_;
    double ly_;
    double lz_;
    double sigma_;
    double a_;
    int size_;

    Coefficients muj_;
    Coefficients vCoeff_;
    Vector b_;
    const ComplexNumber iComplex_;
};

// class InitialVelocityConstructor
// {
// 
//   typedef std::vector<std::pair<int, int> >  Modes;
//   typedef std::vector<double>                Coefficients;
//   typedef std::istream_iterator<double>      IstreamIterator;
// 
//   public:
// 
//     InitialVelocityConstructor(double lx,double ly, double lz):
//     lx_(lx),
//     ly_(ly),
//     lz_(lz),
//     sigma_(0.657014913905734),//sigma eigenfrequency calculated in Matlab code
//     size_(66)//N=10 in Matlab code
//     {
//         int i=0;
// 		cout <<"Lx= "<<lx_<<"Ly= "<<ly_<<"Lz= "<<lz_<<endl;
// //*************creation of antiSymmetric modes*******************
//         for (int n=0;n<=size_;++n)
//         {
//             for (int j=0;j<=n;++j)
//             {
// //         antiSymmetricAlpha(i,1)      = 2*(n-j)+1;
// //         antiSymmetricAlpha(i,2)      = 2*j;
//         
//                 antiSymmetric_.push_back(std::make_pair(2*(n-j)+1, 2*j));
//         
// //         antiSymmetricAlphaPrime(i,2) = 2*(n-j)+1;
// //         antiSymmetricAlphaPrime(i,1) = 2*j;
//         
//                 antiSymmetricPrime_.push_back(std::make_pair(2*j, 2*(n-j)+1));
//             }
//         }      
// //*************creation of antiSymmetric modes*******************
// 
// //*************read coefficients from Matlab output files*******************
//         ifstream p("pCoeff.txt");
//         ifstream pPrime("pPrimeCoeff.txt");
//         ifstream q("qCoeff.txt");
//         ifstream qPrime("qPrimeCoeff.txt");;
// 
//         ifstream r("rCoeff.txt");
//         ifstream rPrime("rPrimeCoeff.txt");;
// 
// 		std::copy(IstreamIterator(p),      IstreamIterator(), std::back_inserter(pAlpha_));
// 		std::copy(IstreamIterator(pPrime), IstreamIterator(), std::back_inserter(pPrimeAlpha_));
// 		std::copy(IstreamIterator(q),      IstreamIterator(), std::back_inserter(qAlpha_));
// 		std::copy(IstreamIterator(qPrime), IstreamIterator(), std::back_inserter(qPrimeAlpha_));
// 		std::copy(IstreamIterator(r),      IstreamIterator(), std::back_inserter(rAlpha_));
// 		std::copy(IstreamIterator(rPrime), IstreamIterator(), std::back_inserter(rPrimeAlpha_));
// 
//         cout <<"Check the coefficients of P"<<endl;
//         copy (pAlpha_.begin(), pAlpha_.end(), std::ostream_iterator<double> (cout, " "));
//         cout <<"*****"<<endl;
// //*************read coefficients from Matlab output files*******************
// // 		createCashedCoefficientsa();
//     }
// // 	void createCashedCoefficientsa()
// // 	{
// // 		
// // 	}
// // 	
//     double getU(const PhysSpacePoint<DIM>& pPhys, double t)const
//     {
//         double u = 0;
//         double x = pPhys[0];
//         double y = pPhys[1];
//         double z = pPhys[2];
//         
//         int    l, k, lPrime, kPrime;
//             
//         for (int i =0; i<size_;++i)
//         {
// 
//         k      = antiSymmetric_[i].first;
//         l      = antiSymmetric_[i].second;
// 
//         kPrime = antiSymmetricPrime_[i].first;
//         lPrime = antiSymmetricPrime_[i].second;
// 
// 
// 
//         u += (pAlpha_[i]*sin(sigma_*t)*(epsilon(k, l)*oneOverNu(k, l)*((k)/lx_)*sin((k*pi*x)/lx_)*cos((l*pi*y)/ly_)));
//         u += (pPrimeAlpha_[i]*cos(sigma_*t)*(epsilon(kPrime, lPrime)*oneOverNu(kPrime, lPrime)*((kPrime)/lx_)*sin((kPrime*pi*x)/lx_)*cos((lPrime*pi*y)/ly_)));
//         u += -(qAlpha_[i]*sin(sigma_*t)*(2*oneOverNu(k, l)*((l)/ly_)*sin((k*pi*x)/lx_)*cos((l*pi*y)/ly_)));
//         u += -(qPrimeAlpha_[i]*cos(sigma_*t)*(2*oneOverNu(kPrime, lPrime)*((lPrime)/ly_)*sin((kPrime*pi*x)/lx_)*cos((lPrime*pi*y)/ly_)));
//         }
// 		u = u * lz_ * cos(pi*z/lz_);
// 		if (u>1000)
// 		{
// 			cout <<"U="<<u<<endl;
// 			cout <<x<<" "<<y<<" "<<z<<endl;
// 		}
//         return u;
//     }
// 
//     double getV(const PhysSpacePoint<DIM>& pPhys, double t)const
//     {
//         double v = 0;
//         double x = pPhys[0];
//         double y = pPhys[1];
//         double z = pPhys[2];
//     
//         int    l, k, lPrime, kPrime;
//     
//         for (int i =0; i<size_;++i)
//         {
//             k      = antiSymmetric_[i].first;
//             l      = antiSymmetric_[i].second;
//             kPrime = antiSymmetricPrime_[i].first;
//             lPrime = antiSymmetricPrime_[i].second;
//     
//             v += (pAlpha_[i]*sin(sigma_*t)*(epsilon(k, l)*oneOverNu(k, l)*((l)/ly_)*cos((k*pi*x)/lx_)*sin((l*pi*y)/ly_)));
//             v += (pPrimeAlpha_[i]*cos(sigma_*t)*(epsilon(kPrime, lPrime)*oneOverNu(kPrime, lPrime)*((lPrime)/ly_)*cos((kPrime*pi*x)/lx_)*sin((lPrime*pi*y)/ly_)));
//             v += (qAlpha_[i]*sin(sigma_*t)*(2*oneOverNu(k, l)*((k)/lx_)*cos((k*pi*x)/lx_)*sin((l*pi*y)/ly_)));
//             v += (qPrimeAlpha_[i]*cos(sigma_*t)*(2*oneOverNu(kPrime, lPrime)*((kPrime)/lx_)*cos((kPrime*pi*x)/lx_)*sin((lPrime*pi*y)/ly_)));
// 			
//         }
// 		v = v * lz_ *cos(pi*z/lz_);
//         return v;
//     }
// 
//     double getW(const PhysSpacePoint<DIM>& pPhys, double t)const
//     {
//         double w = 0;
//         double x = pPhys[0];
//         double y = pPhys[1];
//         double z = pPhys[2];
//     
//         int    l, k, lPrime, kPrime;
//     
//         for (int i =0; i<size_;++i)
//         {
//             k      = antiSymmetric_[i].first;
//             l      = antiSymmetric_[i].second;
//             kPrime = antiSymmetricPrime_[i].first;
//             lPrime = antiSymmetricPrime_[i].second;
//     
//             w += -(sigma_*rAlpha_[i]*sin(sigma_*t)*(epsilon(k, l)*cos((k*pi*x)/lx_)*sin((l*pi*y)/ly_)));
//             w += (sigma_*rPrimeAlpha_[i]*cos(sigma_*t)*(epsilon(kPrime, lPrime)*cos((kPrime*pi*x)/lx_)*cos((lPrime*pi*y)/ly_)));
// 			
//         }
// 		w = w *lz_/pi*sin(pi*z/lz_);
//         return w;
//     }
// 	
// 	double getS1(const PhysSpacePoint<DIM>& pPhys, double t)const
// 	{
// 		double s1 = 0;
// 		double x = pPhys[0];
// 		double y = pPhys[1];
// 		double z = pPhys[2];
// 		int    l, k, lPrime, kPrime;
//     
// 		for (int i =0; i<size_;++i)
// 		{
// 			k      = antiSymmetric_[i].first;
// 			l      = antiSymmetric_[i].second;
// 			kPrime = antiSymmetricPrime_[i].first;
// 			lPrime = antiSymmetricPrime_[i].second;
//     
// 			s1 += (pAlpha_[i]*sigma_*cos(sigma_*t)*(epsilon(k, l)*oneOverNu(k, l)*((k)/lx_)*sin((k*pi*x)/lx_)*cos((l*pi*y)/ly_)));
// 			s1 += -(pPrimeAlpha_[i]*sigma_*sin(sigma_*t)*(epsilon(kPrime, lPrime)*oneOverNu(kPrime, lPrime)*((kPrime)/lx_)*sin((kPrime*pi*x)/lx_)*cos((lPrime*pi*y)/ly_)));
// 			s1 += -(qAlpha_[i]*sigma_*cos(sigma_*t)*(2*oneOverNu(k, l)*((l)/ly_)*sin((k*pi*x)/lx_)*cos((l*pi*y)/ly_)));
// 			s1 += (qPrimeAlpha_[i]*sigma_*sin(sigma_*t)*(2*oneOverNu(kPrime, lPrime)*((lPrime)/ly_)*sin((kPrime*pi*x)/lx_)*cos((lPrime*pi*y)/ly_)));
// 			
// 			s1 += -(pAlpha_[i]*sin(sigma_*t)*(epsilon(k, l)*oneOverNu(k, l)*((l)/ly_)*cos((k*pi*x)/lx_)*sin((l*pi*y)/ly_)));
// 			s1 += -(pPrimeAlpha_[i]*cos(sigma_*t)*(epsilon(kPrime, lPrime)*oneOverNu(kPrime, lPrime)*((lPrime)/ly_)*cos((kPrime*pi*x)/lx_)*sin((lPrime*pi*y)/ly_)));
// 			s1 += -(qAlpha_[i]*sin(sigma_*t)*(2*oneOverNu(k, l)*((k)/lx_)*cos((k*pi*x)/lx_)*sin((l*pi*y)/ly_)));
// 			s1 += -(qPrimeAlpha_[i]*cos(sigma_*t)*(2*oneOverNu(kPrime, lPrime)*((kPrime)/lx_)*cos((kPrime*pi*x)/lx_)*sin((lPrime*pi*y)/ly_)));
// 			
// 			s1 += lz_/pi*sigma_*sigma_*rAlpha_[i]*cos(sigma_*t)*(epsilon(k, l)*k/lx_*sin((k*pi*x)/lx_)*cos((l*pi*y)/ly_));
// 			s1 += lz_/pi*sigma_*sigma_*rPrimeAlpha_[i]*sin(sigma_*t)*(epsilon(kPrime, lPrime)*kPrime/lx_*sin((kPrime*pi*x)/lx_)*cos((lPrime*pi*y)/ly_));
// 			
// 			
// 		}
// 		s1 = s1 *lz_/**cos(pi*z/lz_)*/;
// 		return s1;
// 	}			
// 
// 	double getXDerS1(const PhysSpacePoint<DIM>& pPhys, double t)const
// 	{
// 		double sx1 = 0;
// 		double x = pPhys[0];
// 		double y = pPhys[1];
// 		double z = pPhys[2];
// 		int    l, k, lPrime, kPrime;
//     
// 		for (int i = 0; i < size_;++i)
// 		{
// 			k      = antiSymmetric_[i].first;
// 			l      = antiSymmetric_[i].second;
// 			kPrime = antiSymmetricPrime_[i].first;
// 			lPrime = antiSymmetricPrime_[i].second;
//     
// 			sx1 += (pAlpha_[i]*sigma_*cos(sigma_*t)*(epsilon(k, l)*oneOverNu(k, l)*((k*pi)/lx_)*((k)/lx_)*cos((k*pi*x)/lx_)*cos((l*pi*y)/ly_)));
// 			sx1 += -(pPrimeAlpha_[i]*sigma_*sin(sigma_*t)*(epsilon(kPrime, lPrime)*oneOverNu(kPrime, lPrime)*((kPrime)/lx_)*((kPrime*pi)/lx_)*cos((kPrime*pi*x)/lx_)*cos((lPrime*pi*y)/ly_)));
// 			sx1 += -(qAlpha_[i]*sigma_*cos(sigma_*t)*(2*oneOverNu(k, l)*((l)/ly_)*((k*pi)/lx_)*cos((k*pi*x)/lx_)*cos((l*pi*y)/ly_)));
// 			sx1 += (qPrimeAlpha_[i]*sigma_*sin(sigma_*t)*(2*oneOverNu(kPrime, lPrime)*((lPrime)/ly_)*((kPrime*pi)/lx_)*cos((kPrime*pi*x)/lx_)*cos((lPrime*pi*y)/ly_)));
// 			
// 			sx1 += (pAlpha_[i]*sin(sigma_*t)*(epsilon(k, l)*oneOverNu(k, l)*((l)/ly_)*((k*pi)/lx_)*sin((k*pi*x)/lx_)*sin((l*pi*y)/ly_)));
// 			sx1 += (pPrimeAlpha_[i]*cos(sigma_*t)*(epsilon(kPrime, lPrime)*oneOverNu(kPrime, lPrime)*((lPrime)/ly_)*((kPrime*pi)/lx_)*sin((kPrime*pi*x)/lx_)*sin((lPrime*pi*y)/ly_)));
// 			sx1 += (qAlpha_[i]*sin(sigma_*t)*(2*oneOverNu(k, l)*((k)/lx_)*(k*pi)/lx_)*sin((k*pi*x)/lx_)*sin((l*pi*y)/ly_));
// 			sx1 += (qPrimeAlpha_[i]*cos(sigma_*t)*(2*oneOverNu(kPrime, lPrime)*((kPrime)/lx_)*((kPrime*pi)/lx_)*sin((kPrime*pi*x)/lx_)*sin((lPrime*pi*y)/ly_)));
// 			
// 			sx1 += lz_/pi*sigma_*sigma_*rAlpha_[i]*cos(sigma_*t)*(epsilon(k, l)*k/lx_*((k*pi)/lx_)*cos((k*pi*x)/lx_)*cos((l*pi*y)/ly_));
// 			
// 			sx1 += lz_/pi*sigma_*sigma_*rPrimeAlpha_[i]*sin(sigma_*t)*epsilon(kPrime, lPrime)*kPrime/lx_*((kPrime*pi*x)/lx_)*cos((kPrime*pi*x)/lx_)*cos((lPrime*pi*y)/ly_);
// 			
// 			
// 		}
// 		sx1 = sx1 *lz_/**cos(pi*z/lz_)*/;
// 		return sx1;
// 	}
// 
// 	double getS2(const PhysSpacePoint<DIM>& pPhys, double t)const
// 	{
// 		double s2 = 0;
// 		double x = pPhys[0];
// 		double y = pPhys[1];
// 		double z = pPhys[2];
// 		int    l, k, lPrime, kPrime;
//     
// 		for (int i =0; i<size_;++i)
// 		{
// 			k      = antiSymmetric_[i].first;
// 			l      = antiSymmetric_[i].second;
// 			kPrime = antiSymmetricPrime_[i].first;
// 			lPrime = antiSymmetricPrime_[i].second;
//     
// 			s2 += (pAlpha_[i]*sigma_*cos(sigma_*t)*(epsilon(k, l)*oneOverNu(k, l)*((l)/ly_)*cos((k*pi*x)/lx_)*sin((l*pi*y)/ly_)));
// 			s2 += -(pPrimeAlpha_[i]*sigma_*sin(sigma_*t)*(epsilon(kPrime, lPrime)*oneOverNu(kPrime, lPrime)*((lPrime)/ly_)*cos((kPrime*pi*x)/lx_)*sin((lPrime*pi*y)/ly_)));
// 			s2 += (qAlpha_[i]*sigma_*cos(sigma_*t)*(2*oneOverNu(k, l)*((kPrime)/lx_)*cos((k*pi*x)/lx_)*cos((l*pi*y)/ly_)));
// 			s2 += -(qPrimeAlpha_[i]*sigma_*sin(sigma_*t)*(2*oneOverNu(kPrime, lPrime)*((kPrime)/lx_)*cos((kPrime*pi*x)/lx_)*sin((lPrime*pi*y)/ly_)));
// 			
// 			s2 += (pAlpha_[i]*sin(sigma_*t)*(epsilon(k, l)*oneOverNu(k, l)*((k)/lx_)*sin((k*pi*x)/lx_)*cos((l*pi*y)/ly_)));
// 			s2 += (pPrimeAlpha_[i]*cos(sigma_*t)*(epsilon(kPrime, lPrime)*oneOverNu(kPrime, lPrime)*((kPrime)/lx_)*sin((kPrime*pi*x)/lx_)*cos((lPrime*pi*y)/ly_)));
// 			s2 += -(qAlpha_[i]*sin(sigma_*t)*(2*oneOverNu(k, l)*((l)/ly_)*sin((k*pi*x)/lx_)*cos((l*pi*y)/ly_)));
// 			s2 += -(qPrimeAlpha_[i]*cos(sigma_*t)*(2*oneOverNu(kPrime, lPrime)*((lPrime)/ly_)*sin((kPrime*pi*x)/lx_)*cos((lPrime*pi*y)/ly_)));
// 			
// 			s2 += lz_/pi*sigma_*sigma_*rAlpha_[i]*cos(sigma_*t)*(epsilon(k, l)*l/ly_*cos((k*pi*x)/lx_)*sin((l*pi*y)/ly_));
// 			s2 += lz_/pi*sigma_*sigma_*rPrimeAlpha_[i]*sin(sigma_*t)*(epsilon(kPrime, lPrime)*lPrime/ly_*cos((kPrime*pi*x)/lx_)*sin((lPrime*pi*y)/ly_));
// 			
// 			
// 		}
// 		s2 = s2 *lz_/**cos(pi*z/lz_)*/;
// 		return s2;
// 	}
// 	
// 	double getYDerS2(const PhysSpacePoint<DIM>& pPhys, double t)const
// 	{
// 		double s2 = 0;
// 		double x = pPhys[0];
// 		double y = pPhys[1];
// 		double z = pPhys[2];
// 		int    l, k, lPrime, kPrime;
// 
// 		for (int i =0; i<size_;++i)
// 		{
// 			k      = antiSymmetric_[i].first;
// 			l      = antiSymmetric_[i].second;
// 			kPrime = antiSymmetricPrime_[i].first;
// 			lPrime = antiSymmetricPrime_[i].second;
// 
// 			s2 += (pAlpha_[i]*sigma_*cos(sigma_*t)*(epsilon(k, l)*oneOverNu(k, l)*((l)/ly_)*((l*pi)/ly_)*cos((k*pi*x)/lx_)*cos((l*pi*y)/ly_)));
// 			s2 += -(pPrimeAlpha_[i]*sigma_*sin(sigma_*t)*(epsilon(kPrime, lPrime)*oneOverNu(kPrime, lPrime)*((lPrime)/ly_)*((lPrime*pi*y)/ly_)*cos((kPrime*pi*x)/lx_)*cos((lPrime*pi*y)/ly_)));
// 			s2 += -(qAlpha_[i]*sigma_*cos(sigma_*t)*(2*oneOverNu(k, l)*((kPrime)/lx_)*((l*pi*y)/ly_)*cos((k*pi*x)/lx_)*sin((l*pi*y)/ly_)));
// 			s2 += -(qPrimeAlpha_[i]*sigma_*sin(sigma_*t)*(2*oneOverNu(kPrime, lPrime)*((kPrime)/lx_)*((lPrime*pi*y)/ly_)*cos((kPrime*pi*x)/lx_)*cos((lPrime*pi*y)/ly_)));
// 			
// 			s2 += -(pAlpha_[i]*sin(sigma_*t)*(epsilon(k, l)*oneOverNu(k, l)*((k)/lx_)*((l*pi*y)/ly_)*sin((k*pi*x)/lx_)*sin((l*pi*y)/ly_)));
// 			s2 += -(pPrimeAlpha_[i]*cos(sigma_*t)*(epsilon(kPrime, lPrime)*oneOverNu(kPrime, lPrime)*((kPrime)/lx_)*((lPrime*pi*y)/ly_)*sin((kPrime*pi*x)/lx_)*sin((lPrime*pi*y)/ly_)));
// 			s2 += (qAlpha_[i]*sin(sigma_*t)*(2*oneOverNu(k, l)*((l)/ly_)*((lPrime*pi*y)/ly_)*sin((k*pi*x)/lx_)*sin((l*pi*y)/ly_)));
// 			s2 += (qPrimeAlpha_[i]*cos(sigma_*t)*(2*oneOverNu(kPrime, lPrime)*((lPrime)/ly_)*((lPrime*pi*y)/ly_)*sin((kPrime*pi*x)/lx_)*sin((lPrime*pi*y)/ly_)));
// 
// 			s2 += lz_/pi*sigma_*sigma_*rAlpha_[i]*cos(sigma_*t)*(epsilon(k, l)*l/ly_*((lPrime*pi*y)/ly_)*cos((k*pi*x)/lx_)*cos((l*pi*y)/ly_));
// 			s2 += lz_/pi*sigma_*sigma_*rPrimeAlpha_[i]*sin(sigma_*t)*(epsilon(kPrime, lPrime)*lPrime/ly_*((lPrime*pi*y)/ly_)*cos((kPrime*pi*x)/lx_)*cos((lPrime*pi*y)/ly_));
// 		}
// 		s2 = s2 *lz_*cos(pi*z/lz_);
// 		return s2;
// 	}
// 	
//   private:
//     double epsilon(int k,int l)const
//     {
//       if ((k*l)==0)
//       {
//         return std::sqrt(2);
//       }
//       else
//       {
//         return 2;
//       }
//     }
//     
//     double oneOverNu(int k,int l)const
//     {
//       return 1/(pi*std::sqrt((k*k)/(lx_*lx_)+(l*l)/(ly_*ly_)));
//     }
//     
// 
//   private:
//     double         sigma_;     
//     int            size_;
//     double         lx_;
//     double         ly_;
//     double         lz_;
// 	
// // 	double         uCosCoeff_;
// // 	double         uSinCoeff_;
// 	
// // 	double         vCosCoeff_;
// // 	double         vSinCoeff_;
// 	   
// // 	double         wCosCoeff_;
// // 	double         wSinCoeff_;
// 	
// 	
//     Modes          antiSymmetricPrime_;
//     Modes          antiSymmetric_;
//        
//     Coefficients   pAlpha_;
//     Coefficients   pPrimeAlpha_;
//     Coefficients   qAlpha_;
//     Coefficients   qPrimeAlpha_;
//     Coefficients   rAlpha_;
//     Coefficients   rPrimeAlpha_;
// };
// initial condition for u

class InitCond
{
public:
    
    using PointPhysicalT = ExactSolutionBase::PointPhysicalT;
    using ElementT = Base::Element;
    using PointReferenceT = Geometry::PointReference;
    using ReturnType = LinearAlgebra::MiddleSizeVector;
public:
    
    InitCond(const ExactSolutionBase* init)
            : velocity_(init)
    {
    }
    
    virtual void operator()(const ElementT* element, const PointReferenceT& p, LinearAlgebra::MiddleSizeVector& r) const = 0;

protected:
    const ExactSolutionBase* const velocity_;
};

class InitCondU : public InitCond, public Integration::ElementIntegrandBase<LinearAlgebra::MiddleSizeVector>
{
public:
    using ReturnType = LinearAlgebra::MiddleSizeVector;
public:
    InitCondU(const ExactSolutionBase* init);

    void operator()(const ElementT* element, const PointReferenceT& p, LinearAlgebra::MiddleSizeVector& r) const;
    void elementIntegrand(const ElementT* element, const PointReferenceT& p, LinearAlgebra::MiddleSizeVector& r)
    {
        operator()(element, p, r);
    }
};

// initial condition for v
class InitCondV : public InitCond, public Integration::ElementIntegrandBase<LinearAlgebra::MiddleSizeVector>
{
public:
    
    InitCondV(const ExactSolutionBase* init);

    void operator()(const ElementT* element, const PointReferenceT& p, LinearAlgebra::MiddleSizeVector& r) const;
    void elementIntegrand(const ElementT* element, const PointReferenceT& p, LinearAlgebra::MiddleSizeVector& r)
    {
        operator()(element, p, r);
    }
};
// 

// initial condition for w
class InitCondW : public InitCond, public Integration::ElementIntegrandBase<LinearAlgebra::MiddleSizeVector>
{
public:
    InitCondW(const ExactSolutionBase* init);

    void operator()(const ElementT* element, const PointReferenceT& p, LinearAlgebra::MiddleSizeVector& r) const;
    void elementIntegrand(const ElementT* element, const PointReferenceT& p, LinearAlgebra::MiddleSizeVector& r)
    {
        operator()(element, p, r);
    }
    
};

// initial condition for Lambda, in the compressible set-up used as densitu /rho
class InitCondLambda : public InitCond, public Integration::ElementIntegrandBase<LinearAlgebra::MiddleSizeVector>
{
public:
    
    InitCondLambda(const ExactSolutionBase* init);

    void operator()(const ElementT* element, const PointReferenceT& p, LinearAlgebra::MiddleSizeVector& r) const;
    void elementIntegrand(const ElementT* element, const PointReferenceT& p, LinearAlgebra::MiddleSizeVector& r)
    {
        operator()(element, p, r);
    }
    
};

#endif
//------------------------------------------------------------------------------
