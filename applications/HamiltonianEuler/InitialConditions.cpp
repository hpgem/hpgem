/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2014, University of Twente
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

 1. Redistributions of source code must retain the above copyright notice, this
 list of conditions and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.

 3. Neither the name of the copyright holder nor the names of its contributors
 may be used to endorse or promote products derived from this software without
 specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "InitialConditions.h"
#include "LinearAlgebra/MiddleSizeMatrix.h"

const double ExactSolutionBase::Pi = 3.14159265;

Compressible3DOneThirdPeriodic::Compressible3DOneThirdPeriodic(double lx,
                                                               double ly,
                                                               double lz)
    : ExactSolutionBase(), iComplex_(0, -1), lx_(lx), ly_(ly), lz_(lz) {
    k_ = 1 * 2 * Pi / lx_;
    l_ = 1 * 2 * Pi / ly_;
    sigma_ = 1 / (std::sqrt(k_ * k_ + l_ * l_ + 1));

    cout << "k=" << k_ << ", l=" << l_ << ", sigma=" << sigma_ << endl;
}

double Compressible3DOneThirdPeriodic::getU(const PointPhysicalT& pPhys,
                                            double t) const {
    // 	ComplexNumber u = 0;
    // 	double xx = pPhys[0];
    // 	double yy = pPhys[1];
    // 	double zz = pPhys[2];
    //
    // 	u =
    // iComplex_*k_*(sigma_*sigma_*sigma_)/(1-sigma_*sigma_)*(1+(l_/(sigma_*k_))*(l_/(sigma_*k_)))*sin(k_*xx);
    // 	u =
    // u*(cos(l_*yy)-iComplex_*sin(l_*yy))*cos(Pi*zz/lz_)*(-cos(sigma_*t)+iComplex_*sin(sigma_*t));
    // 	return real(u);

    double u = 0;
    double xx = pPhys[0];
    double yy = pPhys[1];
    double zz = pPhys[2];

    u = k_ * (sigma_ * sigma_ * sigma_) / (1 - sigma_ * sigma_) *
        (1 + (l_ * l_) / (sigma_ * k_ * sigma_ * k_)) * sin(k_ * xx) *
        (-sin(l_ * yy - sigma_ * t)) * cos(2 * Pi * zz / lz_);

    return u;
}

double Compressible3DOneThirdPeriodic::getV(const PointPhysicalT& pPhys,
                                            double t) const {
    // 	ComplexNumber v = 0;
    // 	double xx = pPhys[0];
    // 	double yy = pPhys[1];
    // 	double zz = pPhys[2];
    //
    // 	v
    // =(-l_*sigma_*cos(k_*xx)+1/k_*sin(k_*xx))*(cos(l_*yy)-iComplex_*sin(l_*yy))*cos(Pi*zz/lz_)*(-cos(sigma_*t)+iComplex_*sin(sigma_*t));
    // 	return real(v);

    double v = 0;
    double xx = pPhys[0];
    double yy = pPhys[1];
    double zz = pPhys[2];

    v = (-l_ * sigma_ * cos(k_ * xx) + 1 / k_ * sin(k_ * xx)) *
        (cos(l_ * yy - sigma_ * t)) * cos(2 * Pi * zz / lz_);
    return v;
}

double Compressible3DOneThirdPeriodic::getW(const PointPhysicalT& pPhys,
                                            double t) const {
    // 	ComplexNumber w = 0;
    // 	double xx = pPhys[0];
    // 	double yy = pPhys[1];
    // 	double zz = pPhys[2];
    //
    // 	w =
    // (-iComplex_*sigma_)*(cos(k_*xx)+l_/(k_*sigma_)*sin(k_*xx))*(cos(l_*yy)-iComplex_*sin(l_*yy))*sin(Pi*zz/lz_)*(-cos(sigma_*t)+iComplex_*sin(sigma_*t));
    // 	return real(w);

    double w = 0;
    double xx = pPhys[0];
    double yy = pPhys[1];
    double zz = pPhys[2];

    w = (sigma_) * (cos(k_ * xx) + l_ / (k_ * sigma_) * sin(k_ * xx)) *
        (sin(l_ * yy - sigma_ * t)) * sin(2 * Pi * zz / lz_);
    return w;
}

double Compressible3DOneThirdPeriodic::getP(const PointPhysicalT& pPhys,
                                            double t) const {

    // 	ComplexNumber p = 0;
    // 	double xx = pPhys[0];
    // 	double yy = pPhys[1];
    // 	double zz = pPhys[2];
    //
    // 	p =
    // (-iComplex_*sigma_)*(-iComplex_*sigma_)*(cos(k_*xx)+l_/(k_*sigma_)*sin(k_*xx))*(cos(l_*yy)-iComplex_*sin(l_*yy))*cos(Pi*zz/lz_)*(-cos(sigma_*t)+iComplex_*sin(sigma_*t));
    // 	return real(p);

    double p = 0;
    double xx = pPhys[0];
    double yy = pPhys[1];
    double zz = pPhys[2];

    p = -(sigma_) * (sigma_) *
        (cos(k_ * xx) + l_ / (k_ * sigma_) * sin(k_ * xx)) *
        (cos(l_ * yy - sigma_ * t)) * cos(2 * Pi * zz / lz_);
    return p;
}

Compressible3DPeriodic::Compressible3DPeriodic() : ExactSolutionBase() {}

double Compressible3DPeriodic::getU(const PointPhysicalT& pPhys,
                                    double t) const {
    //	return 1;
    return std::cos(2 * Pi * (t - pPhys[0]));
}

double Compressible3DPeriodic::getV(const PointPhysicalT& pPhys,
                                    double t) const {
    //    return 2;
    return std::cos(2 * Pi * (t - pPhys[1]));
}

double Compressible3DPeriodic::getW(const PointPhysicalT& pPhys,
                                    double t) const {
    //    return 3;
    return std::cos(2 * Pi * (t - pPhys[2]));
}

double Compressible3DPeriodic::getP(const PointPhysicalT& pPhys,
                                    double t) const {
    //    return 4;
    return std::cos(2 * Pi * (t - pPhys[0])) +
           std::cos(2 * Pi * (t - pPhys[1])) +
           std::cos(2 * Pi * (t - pPhys[2]));
}

Compressible3DWalls::Compressible3DWalls() : ExactSolutionBase() {}

double Compressible3DWalls::getU(const PointPhysicalT& pPhys, double t) const {
    return -cos(2 * Pi * (t)) * std::sin(2 * Pi * (pPhys[0]));
}

double Compressible3DWalls::getV(const PointPhysicalT& pPhys, double t) const {
    return -cos(2 * Pi * (t)) * std::sin(2 * Pi * (pPhys[1]));
}

double Compressible3DWalls::getW(const PointPhysicalT& pPhys, double t) const {
    return -cos(2 * Pi * (t)) * std::sin(2 * Pi * (pPhys[2]));
}

double Compressible3DWalls::getP(const PointPhysicalT& pPhys, double t) const {
    return cos(2 * Pi * (pPhys[0])) * sin(2 * Pi * (t)) +
           cos(2 * Pi * (pPhys[1])) * sin(2 * Pi * (t)) +
           cos(2 * Pi * (pPhys[2])) * sin(2 * Pi * (t));
}

Incompressible3DPeriodic::Incompressible3DPeriodic() : ExactSolutionBase() {}

double Incompressible3DPeriodic::getU(const PointPhysicalT& pPhys,
                                      double t) const {
    return (std::sqrt(3) * cos(2 * Pi * (pPhys[0] + pPhys[1] + pPhys[2]) +
                               std::sqrt(3) * t / 3) +
            3 * sin(2 * Pi * (pPhys[0] + pPhys[1] + pPhys[2]) +
                    std::sqrt(3) * t / 3)) /
           Pi / 2;
}

double Incompressible3DPeriodic::getV(const PointPhysicalT& pPhys,
                                      double t) const {
    return (std::sqrt(3) * cos(2 * Pi * (pPhys[0] + pPhys[1] + pPhys[2]) +
                               std::sqrt(3) * t / 3) -
            3 * sin(2 * Pi * (pPhys[0] + pPhys[1] + pPhys[2]) +
                    std::sqrt(3) * t / 3)) /
           Pi / 2;
}

double Incompressible3DPeriodic::getW(const PointPhysicalT& pPhys,
                                      double t) const {
    return -(2 * std::sqrt(3) *
             cos(2 * Pi * (pPhys[0] + pPhys[1] + pPhys[2]) +
                 std::sqrt(3) * t / 3)) /
           Pi / 2;
}

double Incompressible3DPeriodic::getP(const PointPhysicalT& pPhys,
                                      double t) const {
    return 2 *
           cos(2 * Pi * (pPhys[0] + pPhys[1] + pPhys[2]) +
               std::sqrt(3) * t / 3) /
           (4 * Pi * Pi);
}

InitialVelocityConstructorTaylor::InitialVelocityConstructorTaylor()
    : ExactSolutionBase() {}

InitialVelocityConstructorTaylor::InitialVelocityConstructorTaylor(double lx,
                                                                   double ly,
                                                                   double lz)
    : ExactSolutionBase(),
      lx_(lx),
      ly_(ly),
      lz_(lz),
      sigma_(0.6570),  // sigma eigenfrequency calculated in Matlab code
      a_(lx_ / (Pi * sigma_)),
      size_(50),  // N=10 in Matlab code
      muj_(size_ + 1),
      vCoeff_(size_ + 1),
      iComplex_(0, -1) {
    try {

        cout << "Lx= " << lx_ << " Ly= " << ly_ << " Lz= " << lz_ << endl;
        ifstream bCoeff("bCoeff.txt");  // do some checking with existence!!!

        std::copy(IstreamIterator(bCoeff), IstreamIterator(),
                  std::back_inserter(b_));

        //************creation of Symmetric modes*******************
        ComplexNumber a;
        for (int j = 1; j <= size_; ++j) {
            // a = Pi/lx_*std::sqrt((lx_/Pi)*(lx_/Pi) + j*j - a_*a_);
            muj_[j] = Pi / lx_ *
                      std::sqrt(static_cast<ComplexNumber>(
                          (lx_ / Pi) * (lx_ / Pi) + j * j - a_ * a_));

            std::cerr << endl;
            std::cerr << "muj[j]" << muj_[j] << endl;
            std::cerr << "b[j]" << b_[j - 1] << endl;

            if (j % 2 == 0) {
                vCoeff_[j] = -iComplex_ * sinh(ly_ / 2) /
                             sinh(muj_[j] * (ly_ / 2)) * a_ * Pi / lx_ *
                             cos(a_ * Pi / 2) * b_[j - 1];
            } else {
                vCoeff_[j] = -(cosh(ly_ / 2) / cosh(muj_[j] * (ly_ / 2))) *
                             (a_ * Pi) / lx_ * sin(a_ * Pi / 2) * b_[j - 1];
            }
            std::cerr << "vCoeff_[j]" << vCoeff_[j] << endl;
        }
    } catch (const char* c) {
        cout << c;
    }
}

double InitialVelocityConstructorTaylor::getU(const PointPhysicalT& pPhys,
                                              double t) const {
    ComplexNumber u = 0;
    double xx = pPhys[0];
    double yy = pPhys[1];
    double zz = pPhys[2];

    for (int j = 1; j <= size_; ++j) {
        if (j % 2 == 0) {
            u = u +
                (iComplex_ * vCoeff_[j] * (a_ / j - j / a_) *
                 pow(-1, (double)j / 2) * sin(j * (Pi / lx_ * xx - Pi / 2)) *
                 sinh(muj_[j] * (yy - ly_ / 2)));

        } else {
            u = u + (iComplex_ * vCoeff_[j] * (a_ / j - j / a_) *
                     pow(-1, (double)(j - 1) / 2) *
                     cos(j * (Pi / lx_ * xx - Pi / 2)) *
                     cosh(muj_[j] * (yy - ly_ / 2)));
        }
    }
    u = u * cos(Pi * zz / lz_) *
        (cos(sigma_ * t) - iComplex_ * sin(sigma_ * t));
    // cout <<"Constructed u" <<real(u)<<endl;
    return real(u);
}

double InitialVelocityConstructorTaylor::getV(const PointPhysicalT& pPhys,
                                              double t) const {
    ComplexNumber v = 0;
    ComplexNumber bla = 0;
    double xx = pPhys[0];
    double yy = pPhys[1];
    double zz = pPhys[2];
    v = cosh(yy - ly_ / 2) * cos(a_ * (Pi / lx_ * xx - Pi / 2)) -
        iComplex_ * sinh(yy - ly_ / 2) * sin(a_ * (Pi / lx_ * xx - Pi / 2));

    for (int j = 1; j <= size_; ++j) {
        if (j % 2 == 0) {
            v = v +
                (vCoeff_[j] * lx_ / Pi * pow(-1, (double)j / 2) *
                 (iComplex_ * muj_[j] / a_ * cos(j * (Pi / lx_ * xx - Pi / 2)) *
                      cosh(muj_[j] * (yy - ly_ / 2)) +
                  sin(j * (Pi / lx_ * xx - Pi / 2)) *
                      sinh(muj_[j] * (yy - ly_ / 2)) / (ComplexNumber)j));

        } else { /*if (j==3)
                {
                v = v+ (vCoeff_[j]*lx_/Pi*pow(-1,(j -
                1)/2)*(-iComplex_*muj_[j]/a_*sin(j*(Pi/lx_*xx-Pi/2))*sinh(muj_[j]*(yy-ly_/2))
                +
                (cos(j*(Pi/lx_*xx-Pi/2))*cosh(muj_[j]*(yy-ly_/2)))/(ComplexNumber)j));
                }*/
            v = v + (vCoeff_[j] * lx_ / Pi * pow(-1, (double)(j - 1) / 2) *
                     (-iComplex_ * muj_[j] / a_ *
                          sin(j * (Pi / lx_ * xx - Pi / 2)) *
                          sinh(muj_[j] * (yy - ly_ / 2)) +
                      (cos(j * (Pi / lx_ * xx - Pi / 2)) *
                       cosh(muj_[j] * (yy - ly_ / 2))) /
                          (ComplexNumber)j));
        }
    }
    v = real(v) * cos(Pi * zz / lz_) *
        (cos(sigma_ * t) - iComplex_ * sin(sigma_ * t));
    return real(v);
}

double InitialVelocityConstructorTaylor::getW(const PointPhysicalT& pPhys,
                                              double t) const {
    ComplexNumber w = 0;
    double xx = pPhys[0];
    double yy = pPhys[1];
    double zz = pPhys[2];

    w = -iComplex_ * sigma_ *
        (-iComplex_ * (1 / sigma_) *
         (sinh(yy - ly_ / 2) * cos(a_ * (Pi / lx_ * xx - Pi / 2)) -
          iComplex_ * cosh(yy - ly_ / 2) * sin(a_ * (Pi / lx_ * xx - Pi / 2))));

    for (int j = 1; j <= size_; ++j) {
        if (j % 2 == 0) {
            w = w - iComplex_ * sigma_ * (-iComplex_) * (1 / sigma_) *
                        vCoeff_[j] * lx_ / Pi * pow(-1, (double)j / 2) *
                        (iComplex_ * (muj_[j] / a_) *
                             cos(j * (Pi / lx_ * xx - Pi / 2)) * muj_[j] *
                             sinh(muj_[j] * (yy - ly_ / 2)) +
                         (sin(j * (Pi / lx_ * xx - Pi / 2)) * muj_[j] *
                          cosh(muj_[j] * (yy - ly_ / 2))) /
                             (ComplexNumber)j);
            w = w - iComplex_ * sigma_ * (-iComplex_) * (1 / sigma_) *
                        iComplex_ * vCoeff_[j] * (a_ / j - j / a_) *
                        pow(-1, (double)j / 2) * (double)j * (Pi / lx_) *
                        cos(j * (Pi / lx_ * xx - Pi / 2)) *
                        sinh(muj_[j] * (yy - ly_ / 2));

        } else {
            w = w - iComplex_ * sigma_ *
                        (-iComplex_ * (1 / sigma_) * vCoeff_[j] * lx_ / Pi *
                         pow(-1, (double)(j - 1) / 2) *
                         (-iComplex_ * muj_[j] / a_ *
                              sin(j * (Pi / lx_ * xx - Pi / 2)) * muj_[j] *
                              cosh(muj_[j] * (yy - ly_ / 2)) +
                          (cos(j * (Pi / lx_ * xx - Pi / 2)) * muj_[j] *
                           sinh(muj_[j] * (yy - ly_ / 2))) /
                              (ComplexNumber)j));
            w = w - iComplex_ * sigma_ *
                        (-iComplex_ * (1 / sigma_) * iComplex_ * vCoeff_[j] *
                         (a_ / j - j / a_) * pow(-1, (double)(j - 1) / 2) *
                         (-j * Pi / lx_) * sin(j * (Pi / lx_ * xx - Pi / 2)) *
                         cosh(muj_[j] * (yy - ly_ / 2)));
        }
    }
    w = w * sin(Pi * zz / lz_) *
        (cos(sigma_ * t) - iComplex_ * sin(sigma_ * t));
    return real(w);
}

double InitialVelocityConstructorTaylor::getP(const PointPhysicalT& pPhys,
                                              double t) const {
    ComplexNumber p = 0;
    double xx = pPhys[0];
    double yy = pPhys[1];
    double zz = pPhys[2];

    p = (-iComplex_ * sigma_) * (-iComplex_ * sigma_) *
        (-iComplex_ * (1 / sigma_) *
         (sinh(yy - ly_ / 2) * cos(a_ * (Pi / lx_ * xx - Pi / 2)) -
          iComplex_ * cosh(yy - ly_ / 2) * sin(a_ * (Pi / lx_ * xx - Pi / 2))));

    for (int j = 1; j <= size_; ++j) {
        if (j % 2 == 0) {
            p = p - (-iComplex_ * sigma_) * iComplex_ * sigma_ * (-iComplex_) *
                        (1 / sigma_) * vCoeff_[j] * lx_ / Pi *
                        pow(-1, (double)j / 2) *
                        (iComplex_ * (muj_[j] / a_) *
                             cos(j * (Pi / lx_ * xx - Pi / 2)) * muj_[j] *
                             sinh(muj_[j] * (yy - ly_ / 2)) +
                         (sin(j * (Pi / lx_ * xx - Pi / 2)) * muj_[j] *
                          cosh(muj_[j] * (yy - ly_ / 2))) /
                             (ComplexNumber)(ComplexNumber)j);
            p = p - (-iComplex_ * sigma_) * iComplex_ * sigma_ * (-iComplex_) *
                        (1 / sigma_) * iComplex_ * vCoeff_[j] *
                        (a_ / j - j / a_) * pow(-1, (double)j / 2) * (double)j *
                        (Pi / lx_) * cos(j * (Pi / lx_ * xx - Pi / 2)) *
                        sinh(muj_[j] * (yy - ly_ / 2));

        } else {
            p = p - (-iComplex_ * sigma_) * iComplex_ * sigma_ *
                        (-iComplex_ * (1 / sigma_) * vCoeff_[j] * lx_ / Pi *
                         pow(-1, (double)(j - 1) / 2) *
                         (-iComplex_ * muj_[j] / a_ *
                              sin(j * (Pi / lx_ * xx - Pi / 2)) * muj_[j] *
                              cosh(muj_[j] * (yy - ly_ / 2)) +
                          (cos(j * (Pi / lx_ * xx - Pi / 2)) * muj_[j] *
                           sinh(muj_[j] * (yy - ly_ / 2))) /
                              (ComplexNumber)j));
            p = p - (-iComplex_ * sigma_) * iComplex_ * sigma_ *
                        (-iComplex_ * (1 / sigma_) * iComplex_ * vCoeff_[j] *
                         (a_ / j - j / a_) * pow(-1, (double)(j - 1) / 2) *
                         (-j * Pi / lx_) * sin(j * (Pi / lx_ * xx - Pi / 2)) *
                         cosh(muj_[j] * (yy - ly_ / 2)));
        }
    }
    p = p * cos(Pi * zz / lz_) *
        (cos(sigma_ * t) - iComplex_ * sin(sigma_ * t));
    return real(p);
}

InitCondU::InitCondU(const ExactSolutionBase* init) : InitCond(init) {}

void InitCondU::operator()(const ElementT* element, const PointReferenceT& pRef,
                           LinearAlgebra::MiddleSizeVector& r) const {

    PointPhysicalT pPhys =
        element->referenceToPhysical(pRef);  // ...transform the point.
                                             //
    unsigned int numberOfDegreesOfFreedom = element->getNrOfBasisFunctions();
    //
    double uVal = InitCond::velocity_->getU(pPhys, 0);
    //

    // cout << "**************************"<<numberOfDegreesOfFreedom<<endl;
    for (unsigned int j = 0; j < numberOfDegreesOfFreedom; ++j) {
        //
        r[j] = element->basisFunction(j, pRef) * uVal;
    }
}

// initial condition for v
InitCondV::InitCondV(const ExactSolutionBase* init) : InitCond(init) {}

void InitCondV::operator()(const ElementT* element, const PointReferenceT& pRef,
                           LinearAlgebra::MiddleSizeVector& r) const {
    PointPhysicalT pPhys =
        element->referenceToPhysical(pRef);  // ...transform the point.

    unsigned int numberOfDegreesOfFreedom = element->getNrOfBasisFunctions();

    double vVal = InitCond::velocity_->getV(pPhys, 0);

    for (unsigned int j = 0; j < numberOfDegreesOfFreedom; ++j) {
        r(j) = element->basisFunction(j, pRef) * vVal;
    }
}

// initial condition for w
InitCondW::InitCondW(const ExactSolutionBase* init) : InitCond(init) {}

void InitCondW::operator()(const ElementT* element, const PointReferenceT& pRef,
                           LinearAlgebra::MiddleSizeVector& r) const {

    PointPhysicalT pPhys =
        element->referenceToPhysical(pRef);  // ...transform the point.

    unsigned int numberOfDegreesOfFreedom = element->getNrOfBasisFunctions();

    double wVal = InitCond::velocity_->getW(pPhys, 0);

    for (unsigned int j = 0; j < numberOfDegreesOfFreedom; ++j) {
        r(j) = element->basisFunction(j, pRef) * wVal;
    }
}

InitCondLambda::InitCondLambda(const ExactSolutionBase* init)
    : InitCond(init) {}

void InitCondLambda::operator()(const ElementT* element,
                                const PointReferenceT& pRef,
                                LinearAlgebra::MiddleSizeVector& r) const {
    PointPhysicalT pPhys =
        element->referenceToPhysical(pRef);  // ...transform the point.

    unsigned int numberOfDegreesOfFreedom = element->getNrOfBasisFunctions();

    double pVal = InitCond::velocity_->getP(pPhys, 0);

    for (unsigned int j = 0; j < numberOfDegreesOfFreedom; ++j) {
        r(j) = element->basisFunction(j, pRef) * pVal;
    }
}
