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

#include "BasisFunctionCollection_Curl.h"
#include "LinearAlgebra/NumericalVector.h"
#include "Base/Element.h"
#include "Base/BasisFunctionSet.h"
#include "Geometry/PointReference.h"
#include "Geometry/PointPhysical.h"
#include "Geometry/ReferenceGeometry.h"
#include "Base/ElementCacheData.h"

void OuterProduct(const LinearAlgebra::NumericalVector& a, const LinearAlgebra::NumericalVector& b, LinearAlgebra::NumericalVector& ret)
{
    ret.resize(3);
    //direct computation using the definition
    ret[0] = a[1] * b[2] - a[2] * b[1];
    ret[1] = a[2] * b[0] - a[0] * b[2];
    ret[2] = a[0] * b[1] - a[1] * b[0];
}

double Base::LegendrePolynomial(int degree, double x)
{
    switch (degree)
    {
        case 0:
            return 1;
            break;
        case 1:
            return x;
            break;
        default:
            return double(2 * degree - 1) / double(degree) * x * LegendrePolynomial(degree - 1, x) - double(degree - 1) / double(degree) * LegendrePolynomial(degree - 2, x);
            break;
    }
}

double Base::LegendrePolynomialDerivative(int degree, double x)
{
    switch (degree)
    {
        case 0:
            return 0;
            break;
        case 1:
            return 1;
            break;
        default:
            return double(2 * degree - 1) / double(degree) * (x * LegendrePolynomialDerivative(degree - 1, x) + LegendrePolynomial(degree - 1, x)) - double(degree - 1) / double(degree) * LegendrePolynomialDerivative(degree - 2, x);
            break;
    }
}

Base::Basis_Curl_Bari::Basis_Curl_Bari(int vertex)
        : BaseBasisFunction(), VertexNr(vertex)
{
}
double Base::Basis_Curl_Bari::eval(const Base::BaseBasisFunction::PointReferenceT& p) const
{
    return (VertexNr == 0) ? (1 - p[0] - p[1] - p[2]) : (p[VertexNr - 1]);
}
double Base::Basis_Curl_Bari::evalDeriv0(const Base::BaseBasisFunction::PointReferenceT& p) const
{
    return (VertexNr == 1) - (VertexNr == 0);
}
double Base::Basis_Curl_Bari::evalDeriv1(const Base::BaseBasisFunction::PointReferenceT& p) const
{
    return (VertexNr == 2) - (VertexNr == 0);
}
double Base::Basis_Curl_Bari::evalDeriv2(const Base::BaseBasisFunction::PointReferenceT& p) const
{
    return (VertexNr == 3) - (VertexNr == 0);
}

void Base::Basis_Curl_Bari::evalDeriv(const Base::BaseBasisFunction::PointReferenceT& p, LinearAlgebra::NumericalVector& ret)
{
    ret.resize(3);
    if (VertexNr == 0)
    {
        ret[0] = -1;
        ret[1] = -1;
        ret[2] = -1;
    }
    else
    {
        //clear the return vector so we don't return trash
        ret[0] = 0;
        ret[1] = 0;
        ret[2] = 0;
        ret[VertexNr - 1] = 1;
    }
}

Base::Basis_Curl_Bari Base::Basis_Curl_Bari::barycentricFunctions[] = {Basis_Curl_Bari(0), Basis_Curl_Bari(1), Basis_Curl_Bari(2), Basis_Curl_Bari(3)};

double Base::threeDBasisFunction::eval(const Base::BaseBasisFunction::PointReferenceT& p) const
{
    std::cout << "This type of basisfunctions expects to return vectors, but the return type of this function expects scalars";
    exit(-1);
}
double Base::threeDBasisFunction::evalDeriv0(const Base::BaseBasisFunction::PointReferenceT& p) const
{
    std::cout << "The derivative of a vector-function is not implemented. Perhaps you meant evalCurl?";
    exit(-1);
}
double Base::threeDBasisFunction::evalDeriv1(const Base::BaseBasisFunction::PointReferenceT& p) const
{
    std::cout << "The derivative of a vector-function is not implemented. Perhaps you meant evalCurl?";
    exit(-1);
}
double Base::threeDBasisFunction::evalDeriv2(const Base::BaseBasisFunction::PointReferenceT& p) const
{
    std::cout << "The derivative of a vector-function is not implemented. Perhaps you meant evalCurl?";
    exit(-1);
}

Base::Basis_Curl_Edge::Basis_Curl_Edge(int degree, int localFirstVertex, int localSecondVertex)
        : deg(degree), o(localFirstVertex), i(localSecondVertex)
{
}

void Base::Basis_Curl_Edge::eval(const Base::BaseBasisFunction::PointReferenceT& p, LinearAlgebra::NumericalVector& ret) const
{
    LinearAlgebra::NumericalVector dummy(3); //dummy vectors are used to store parial results
    Basis_Curl_Bari::barycentricFunctions[o].evalDeriv(p, ret);
    Basis_Curl_Bari::barycentricFunctions[i].evalDeriv(p, dummy);
    ret *= Basis_Curl_Bari::barycentricFunctions[i].eval(p); //use only *= and += and so on near numerical vectors or you will be rapidly creating
    dummy *= Basis_Curl_Bari::barycentricFunctions[o].eval(p); //and destroying a gazilion of them, which is a waste of time
    //switch between the special cases in the definition
    switch (deg)
    {
        case 0:
            ret -= dummy;
            break;
        case 1:
            ret += dummy;
            break;
        default:
            //function values at points i and o
            double valI(Basis_Curl_Bari::barycentricFunctions[i].eval(p)), valO(Basis_Curl_Bari::barycentricFunctions[o].eval(p));
            ret -= dummy; //degree 0 value
            dummy *= 2;
            dummy += ret; //degree 1 value
            ret *= -double(deg - 1) / double(deg) * LegendrePolynomial(deg - 2, valI - valO);
            dummy *= double(2 * deg - 1) / double(deg) * LegendrePolynomial(deg - 1, valI - valO);
            ret += dummy;
    }
}

void Base::Basis_Curl_Edge::evalCurl(const Base::BaseBasisFunction::PointReferenceT& p, LinearAlgebra::NumericalVector& ret) const
{
    LinearAlgebra::NumericalVector dummy(3), dummy2(3);
    ret.resize(3);
    switch (deg)
    {
        case 0:
            Basis_Curl_Bari::barycentricFunctions[o].evalDeriv(p, dummy);
            Basis_Curl_Bari::barycentricFunctions[i].evalDeriv(p, dummy2);
            OuterProduct(dummy2, dummy, ret);
            ret *= 2;
            break;
        case 1:
            ret[0] = 0;
            ret[1] = 0;
            ret[2] = 0;
            break;
        default:
            Basis_Curl_Bari::barycentricFunctions[o].evalDeriv(p, dummy);
            Basis_Curl_Bari::barycentricFunctions[i].evalDeriv(p, dummy2);
            LinearAlgebra::NumericalVector dummy3(3), dummy4(dummy2);
            dummy4 -= dummy; //dummy4=nambla(labda_i-labda_o)
            double valI(Basis_Curl_Bari::barycentricFunctions[i].eval(p)), valO(Basis_Curl_Bari::barycentricFunctions[o].eval(p));
            OuterProduct(dummy2, dummy, ret);
            ret *= -double(deg - 1) / double(deg) * LegendrePolynomial(deg - 2, valI - valO) * 2;
            dummy *= valI;
            dummy2 *= valO;
            dummy += dummy2; //dummy=phi_1
            OuterProduct(dummy4, dummy, dummy3);
            dummy3 *= double(2 * deg - 1) / double(deg) * LegendrePolynomialDerivative(deg - 1, valI - valO);
            ret += dummy3;
            dummy2 *= 2;
            dummy -= dummy2; //dummy=phi_0
            OuterProduct(dummy4, dummy, dummy3);
            dummy3 *= double(deg - 1) / double(deg) * LegendrePolynomialDerivative(deg - 2, valI - valO);
            ret -= dummy3;
    }
}

void Base::Basis_Curl_Edge::getReasonableNode(const Base::Element& element, Geometry::PointPhysical node)
{
    PointReferenceT leftnode(3), rightnode(3);
    element.getReferenceGeometry()->getNode(i, leftnode);
    element.getReferenceGeometry()->getNode(o, rightnode);
    element.referenceToPhysical((leftnode + rightnode) * .5, node);
}

Base::Basis_Curl_Edge_Face::Basis_Curl_Edge_Face(int degree, int localOpposingVertex, int localSpecialVertex)
        : deg(degree), c(localSpecialVertex)
{
    a = 0;
    //find the edge (a,b) this functions is based on
    while (c == a || localOpposingVertex == a)
    {
        a++;
    }
    b = a + 1;
    while (c == b || localOpposingVertex == b)
    {
        b++;
    }
}

void Base::Basis_Curl_Edge_Face::eval(const Base::BaseBasisFunction::PointReferenceT& p, LinearAlgebra::NumericalVector& ret) const
{
    Basis_Curl_Bari::barycentricFunctions[c].evalDeriv(p, ret);
    double valA(Basis_Curl_Bari::barycentricFunctions[a].eval(p)), valB(Basis_Curl_Bari::barycentricFunctions[b].eval(p));
    ret *= valA * valB * LegendrePolynomial(deg - 2, valB - valA);
}

void Base::Basis_Curl_Edge_Face::evalCurl(const Base::BaseBasisFunction::PointReferenceT& p, LinearAlgebra::NumericalVector& ret) const
{
    LinearAlgebra::NumericalVector dummy(3), dummy2(3), dummy3(3);
    ret.resize(3);
    double valA(Basis_Curl_Bari::barycentricFunctions[a].eval(p)), valB(Basis_Curl_Bari::barycentricFunctions[b].eval(p));
    Basis_Curl_Bari::barycentricFunctions[a].evalDeriv(p, dummy);
    Basis_Curl_Bari::barycentricFunctions[b].evalDeriv(p, dummy2);
    Basis_Curl_Bari::barycentricFunctions[c].evalDeriv(p, dummy3);
    OuterProduct(dummy, dummy3, ret);
    OuterProduct(dummy2, dummy3, dummy);
    ret *= valB * LegendrePolynomial(deg - 2, valB - valA) - valA * valB * LegendrePolynomialDerivative(deg - 2, valB - valA);
    dummy *= valA * LegendrePolynomial(deg - 2, valB - valA) + valA * valB * LegendrePolynomialDerivative(deg - 2, valB - valA);
    ret += dummy;
}

void Base::Basis_Curl_Edge_Face::getReasonableNode(const Base::Element& element, Geometry::PointPhysical node)
{
    PointReferenceT leftnode(3), rightnode(3);
    element.getReferenceGeometry()->getNode(a, leftnode);
    element.getReferenceGeometry()->getNode(b, rightnode);
    element.referenceToPhysical((leftnode + rightnode) * .5, node);
}

Base::Basis_Curl_Face::Basis_Curl_Face(int degree1, int degree2, int localOpposingVertex, int direction)
        : deg1(degree1), deg2(degree2)
{
    a = 0;
    //construct the face this function lives on
    if (localOpposingVertex == a)
    {
        a++;
    }
    b = a + 1;
    if (localOpposingVertex == b)
    {
        b++;
    }
    c = b + 1;
    if (localOpposingVertex == c)
    {
        c++;
    }
    //lets the tangent always point to vertex b in the rest of this class
    if (direction == 2)
    { //this also swaps the use of l and m in the ainsworth definition
        int temp = c; //but that should be inconsequential
        c = b;
        b = temp;
    }
}

void Base::Basis_Curl_Face::eval(const Base::BaseBasisFunction::PointReferenceT& p, LinearAlgebra::NumericalVector& ret) const
{
    ret.resize(3);
    double valA(Basis_Curl_Bari::barycentricFunctions[a].eval(p)), valB(Basis_Curl_Bari::barycentricFunctions[b].eval(p)), valC(Basis_Curl_Bari::barycentricFunctions[c].eval(p));
    ret[0] = (b == 1 - a == 1);
    ret[1] = (b == 2);
    ret[2] = (b == 3);
    ret *= valA * valB * valC * LegendrePolynomial(deg1, valB - valA) * LegendrePolynomial(deg2, valC - valA);
}

void Base::Basis_Curl_Face::evalCurl(const Base::BaseBasisFunction::PointReferenceT& p, LinearAlgebra::NumericalVector& ret) const
{
    double valA(Basis_Curl_Bari::barycentricFunctions[a].eval(p)), valB(Basis_Curl_Bari::barycentricFunctions[b].eval(p)), valC(Basis_Curl_Bari::barycentricFunctions[c].eval(p));
    LinearAlgebra::NumericalVector dummy(3), dummy2(3), dummy3(3);
    ret.resize(3);
    dummy[0] = (b == 1 - a == 1);
    dummy[1] = (b == 2);
    dummy[2] = (b == 3);
    Basis_Curl_Bari::barycentricFunctions[a].evalDeriv(p, dummy2);
    OuterProduct(dummy2, dummy, ret);
    ret *= valB * valC * LegendrePolynomial(deg1, valB - valA) * LegendrePolynomial(deg2, valC - valA) - valA * valB * valC * (LegendrePolynomialDerivative(deg1, valB - valA) * LegendrePolynomial(deg2, valC - valA) + LegendrePolynomial(deg1, valB - valA) * LegendrePolynomialDerivative(deg2, valC - valA));
    Basis_Curl_Bari::barycentricFunctions[b].evalDeriv(p, dummy2);
    OuterProduct(dummy2, dummy, dummy3);
    dummy3 *= valA * valC * LegendrePolynomial(deg1, valB - valA) * LegendrePolynomial(deg2, valC - valA) + valA * valB * valC * LegendrePolynomialDerivative(deg1, valB - valA) * LegendrePolynomial(deg2, valC - valA);
    ret += dummy3;
    Basis_Curl_Bari::barycentricFunctions[c].evalDeriv(p, dummy2);
    OuterProduct(dummy2, dummy, dummy3);
    dummy3 *= valA * valB * LegendrePolynomial(deg1, valB - valA) * LegendrePolynomial(deg2, valC - valA) + valA * valB * valC * LegendrePolynomial(deg1, valB - valA) * LegendrePolynomialDerivative(deg2, valC - valA);
    ret += dummy3;
}

void Base::Basis_Curl_Face::getReasonableNode(const Base::Element& element, Geometry::PointPhysical node)
{
    PointReferenceT leftnode(3), rightnode(3), Cnode(3);
    element.getReferenceGeometry()->getNode(a, leftnode);
    element.getReferenceGeometry()->getNode(b, rightnode);
    element.getReferenceGeometry()->getNode(c, Cnode);
    element.referenceToPhysical((leftnode + rightnode + Cnode) * .33333, node);
}

Base::Basis_Curl_Face_interior::Basis_Curl_Face_interior(int degree1, int degree2, int localOpposingVertex)
        : deg1(degree1), deg2(degree2), d(localOpposingVertex)
{
    a = 0;
    //find the face...
    if (d == a)
    {
        a++;
    }
    b = a + 1;
    if (d == b)
    {
        b++;
    }
    c = b + 1;
    if (d == c)
    {
        c++;
    }
}

void Base::Basis_Curl_Face_interior::eval(const Base::BaseBasisFunction::PointReferenceT& p, LinearAlgebra::NumericalVector& ret) const
{
    double valA(Basis_Curl_Bari::barycentricFunctions[a].eval(p)), valB(Basis_Curl_Bari::barycentricFunctions[b].eval(p)), valC(Basis_Curl_Bari::barycentricFunctions[c].eval(p));
    Basis_Curl_Bari::barycentricFunctions[d].evalDeriv(p, ret);
    ret *= valA * valB * valC * LegendrePolynomial(deg1, valB - valA) * LegendrePolynomial(deg2, valC - valA);
}

void Base::Basis_Curl_Face_interior::evalCurl(const Base::BaseBasisFunction::PointReferenceT& p, LinearAlgebra::NumericalVector& ret) const
{
    double valA(Basis_Curl_Bari::barycentricFunctions[a].eval(p)), valB(Basis_Curl_Bari::barycentricFunctions[b].eval(p)), valC(Basis_Curl_Bari::barycentricFunctions[c].eval(p));
    LinearAlgebra::NumericalVector dummy(3), dummy2(3), dummy3(3);
    ret.resize(3);
    Basis_Curl_Bari::barycentricFunctions[d].evalDeriv(p, dummy);
    Basis_Curl_Bari::barycentricFunctions[a].evalDeriv(p, dummy2);
    OuterProduct(dummy2, dummy, ret);
    ret *= valB * valC * LegendrePolynomial(deg1, valB - valA) * LegendrePolynomial(deg2, valC - valA) - valA * valB * valC * (LegendrePolynomialDerivative(deg1, valB - valA) * LegendrePolynomial(deg2, valC - valA) + LegendrePolynomial(deg1, valB - valA) * LegendrePolynomialDerivative(deg2, valC - valA));
    Basis_Curl_Bari::barycentricFunctions[b].evalDeriv(p, dummy2);
    OuterProduct(dummy2, dummy, dummy3);
    dummy3 *= valA * valC * LegendrePolynomial(deg1, valB - valA) * LegendrePolynomial(deg2, valC - valA) + valA * valB * valC * LegendrePolynomialDerivative(deg1, valB - valA) * LegendrePolynomial(deg2, valC - valA);
    ret += dummy3;
    Basis_Curl_Bari::barycentricFunctions[c].evalDeriv(p, dummy2);
    OuterProduct(dummy2, dummy, dummy3);
    dummy3 *= valA * valB * LegendrePolynomial(deg1, valB - valA) * LegendrePolynomial(deg2, valC - valA) + valA * valB * valC * LegendrePolynomial(deg1, valB - valA) * LegendrePolynomialDerivative(deg2, valC - valA);
    ret += dummy3;
}

void Base::Basis_Curl_Face_interior::getReasonableNode(const Base::Element& element, Geometry::PointPhysical node)
{
    PointReferenceT leftnode(3), rightnode(3), Cnode(3);
    element.getReferenceGeometry()->getNode(a, leftnode);
    element.getReferenceGeometry()->getNode(b, rightnode);
    element.getReferenceGeometry()->getNode(c, Cnode);
    element.referenceToPhysical((leftnode + rightnode + Cnode) * .3333333333, node);
}

Base::Basis_Curl_interior::Basis_Curl_interior(int degree1, int degree2, int degree3, int direction)
        : deg1(degree1), deg2(degree2), deg3(degree3), direction(direction)
{
}

void Base::Basis_Curl_interior::eval(const Base::BaseBasisFunction::PointReferenceT& p, LinearAlgebra::NumericalVector& ret) const
{
    ret.resize(3);
    ret[0] = 0;
    ret[1] = 0;
    ret[2] = 0;
    ret[direction - 1] = 1;
    double val0(Basis_Curl_Bari::barycentricFunctions[0].eval(p)), val1(Basis_Curl_Bari::barycentricFunctions[1].eval(p)), val2(Basis_Curl_Bari::barycentricFunctions[2].eval(p)), val3(Basis_Curl_Bari::barycentricFunctions[3].eval(p));
    ret *= val0 * val1 * val2 * val3 * LegendrePolynomial(deg1, val1 - val0) * LegendrePolynomial(deg2, val2 - val0) * LegendrePolynomial(deg3, val3 - val0);
}

void Base::Basis_Curl_interior::evalCurl(const Base::BaseBasisFunction::PointReferenceT& p, LinearAlgebra::NumericalVector& ret) const
{
    double val0(Basis_Curl_Bari::barycentricFunctions[0].eval(p)), val1(Basis_Curl_Bari::barycentricFunctions[1].eval(p)), val2(Basis_Curl_Bari::barycentricFunctions[2].eval(p)), val3(Basis_Curl_Bari::barycentricFunctions[3].eval(p));
    LinearAlgebra::NumericalVector dummy(3), dummy2(3), dummy3(3);
    dummy[0] = 0;
    dummy[1] = 0;
    dummy[2] = 0;
    dummy[direction - 1] = 1;
    Basis_Curl_Bari::barycentricFunctions[0].evalDeriv(p, dummy2);
    ret.resize(3);
    OuterProduct(dummy2, dummy, ret);
    ret *= val1 * val2 * val3 * LegendrePolynomial(deg1, val1 - val0) * LegendrePolynomial(deg2, val2 - val0) * LegendrePolynomial(deg3, val3 - val0) - val0 * val1 * val2 * val3 * (LegendrePolynomialDerivative(deg1, val1 - val0) * LegendrePolynomial(deg2, val2 - val0) * LegendrePolynomial(deg3, val3 - val0) + LegendrePolynomial(deg1, val1 - val0) * LegendrePolynomialDerivative(deg2, val2 - val0) * LegendrePolynomial(deg3, val3 - val0) + LegendrePolynomial(deg1, val1 - val0) * LegendrePolynomial(deg2, val2 - val0) * LegendrePolynomialDerivative(deg3, val3 - val0));
    Basis_Curl_Bari::barycentricFunctions[1].evalDeriv(p, dummy2);
    OuterProduct(dummy2, dummy, dummy3);
    dummy3 *= val0 * val2 * val3 * LegendrePolynomial(deg1, val1 - val0) * LegendrePolynomial(deg2, val2 - val0) * LegendrePolynomial(deg3, val3 - val0) + val0 * val1 * val2 * val3 * LegendrePolynomialDerivative(deg1, val1 - val0) * LegendrePolynomial(deg2, val2 - val0) * LegendrePolynomial(deg3, val3 - val0);
    ret += dummy3;
    Basis_Curl_Bari::barycentricFunctions[2].evalDeriv(p, dummy2);
    OuterProduct(dummy2, dummy, dummy3);
    dummy3 *= val0 * val1 * val3 * LegendrePolynomial(deg1, val1 - val0) * LegendrePolynomial(deg2, val2 - val0) * LegendrePolynomial(deg3, val3 - val0) + val0 * val1 * val2 * val3 * LegendrePolynomial(deg1, val1 - val0) * LegendrePolynomialDerivative(deg2, val2 - val0) * LegendrePolynomial(deg3, val3 - val0);
    ret += dummy3;
    Basis_Curl_Bari::barycentricFunctions[3].evalDeriv(p, dummy2);
    OuterProduct(dummy2, dummy, dummy3);
    dummy3 *= val0 * val1 * val2 * LegendrePolynomial(deg1, val1 - val0) * LegendrePolynomial(deg2, val2 - val0) * LegendrePolynomial(deg3, val3 - val0) + val0 * val1 * val2 * val3 * LegendrePolynomial(deg1, val1 - val0) * LegendrePolynomial(deg2, val2 - val0) * LegendrePolynomialDerivative(deg3, val3 - val0);
    ret += dummy3;
}

void Base::Basis_Curl_interior::getReasonableNode(const Base::Element& element, Geometry::PointPhysical node)
{
    PointReferenceT center(3);
    element.getReferenceGeometry()->getCenter(center);
    element.referenceToPhysical(center, node);
}

MyMeshManipulator::MyMeshManipulator(const Base::ConfigurationData* data, bool xPer, bool yPer, bool zPer, std::size_t order, std::size_t idRangeBegin, int nrOfElementMatrixes, int nrOfElementVectors, int nrOfFaceMatrixes, int nrOfFaceVectors)
        : Base::MeshManipulator(data, xPer, yPer, zPer, order, idRangeBegin, nrOfElementMatrixes, nrOfElementVectors, nrOfFaceMatrixes, nrOfFaceVectors)
{
    //std::cout<<nrOfElementVectors<<std::endl;
    createBasisFunctions(order);
    
}

void MyMeshManipulator::createBasisFunctions(unsigned int order)
{
    Base::BasisFunctionSet* bFset = new Base::BasisFunctionSet(order);
    for (int p = 0; p <= order; ++p)
    {
        //constructor takes first the degree of the function then the two vertices on the edge
        bFset->addBasisFunction(new Base::Basis_Curl_Edge(p, 0, 1));
        bFset->addBasisFunction(new Base::Basis_Curl_Edge(p, 0, 2));
        bFset->addBasisFunction(new Base::Basis_Curl_Edge(p, 0, 3));
        bFset->addBasisFunction(new Base::Basis_Curl_Edge(p, 2, 3));
        bFset->addBasisFunction(new Base::Basis_Curl_Edge(p, 1, 3));
        bFset->addBasisFunction(new Base::Basis_Curl_Edge(p, 1, 2));
        //	    std::cout<<"constructed edge functions with p="<<p<<std::endl;
    }
    for (int p = 2; p <= order; ++p)
    {
        //constructor takes first the degree of the function the the vertex opposing the face, then the vertex not on the edge
        bFset->addBasisFunction(new Base::Basis_Curl_Edge_Face(p, 3, 0));
        bFset->addBasisFunction(new Base::Basis_Curl_Edge_Face(p, 3, 1));
        bFset->addBasisFunction(new Base::Basis_Curl_Edge_Face(p, 3, 2));
        bFset->addBasisFunction(new Base::Basis_Curl_Edge_Face(p, 2, 0));
        bFset->addBasisFunction(new Base::Basis_Curl_Edge_Face(p, 2, 1));
        bFset->addBasisFunction(new Base::Basis_Curl_Edge_Face(p, 2, 3));
        bFset->addBasisFunction(new Base::Basis_Curl_Edge_Face(p, 1, 0));
        bFset->addBasisFunction(new Base::Basis_Curl_Edge_Face(p, 1, 2));
        bFset->addBasisFunction(new Base::Basis_Curl_Edge_Face(p, 1, 3));
        bFset->addBasisFunction(new Base::Basis_Curl_Edge_Face(p, 0, 1));
        bFset->addBasisFunction(new Base::Basis_Curl_Edge_Face(p, 0, 2));
        bFset->addBasisFunction(new Base::Basis_Curl_Edge_Face(p, 0, 3));
        //	    std::cout<<"constructed edge based face functions with p="<<p<<std::endl;
    }
    for (int l = 0; l + 3 <= order; ++l)
    {
        for (int m = 0; m + l + 3 <= order; ++m)
        {
            //constructor takes first the degrees of the legendre polinomials, then the vertex opposing the face, then a 1 or 2 to fix the direction of the bubble function
            bFset->addBasisFunction(new Base::Basis_Curl_Face(l, m, 3, 1));
            bFset->addBasisFunction(new Base::Basis_Curl_Face(l, m, 3, 2));
            bFset->addBasisFunction(new Base::Basis_Curl_Face(l, m, 2, 1));
            bFset->addBasisFunction(new Base::Basis_Curl_Face(l, m, 2, 2));
            bFset->addBasisFunction(new Base::Basis_Curl_Face(l, m, 1, 1));
            bFset->addBasisFunction(new Base::Basis_Curl_Face(l, m, 1, 2));
            bFset->addBasisFunction(new Base::Basis_Curl_Face(l, m, 0, 1));
            bFset->addBasisFunction(new Base::Basis_Curl_Face(l, m, 0, 2));
            //constructor takes first the degrees of the legendre polinomials, then the vertex opposing the face
            bFset->addBasisFunction(new Base::Basis_Curl_Face_interior(l, m, 3));
            bFset->addBasisFunction(new Base::Basis_Curl_Face_interior(l, m, 2));
            bFset->addBasisFunction(new Base::Basis_Curl_Face_interior(l, m, 1));
            bFset->addBasisFunction(new Base::Basis_Curl_Face_interior(l, m, 0));
            //	    std::cout<<"constructed face functions and face based interior functions with l="<<l<<" and m="<<m<<std::endl;
        }
    }
    for (int l = 0; l + 4 <= order; ++l)
    {
        for (int m = 0; m + l + 4 <= order; ++m)
        {
            for (int n = 0; n + m + l + 4 <= order; ++n)
            {
                //constructor takes first the degrees of the legendre polinomials, then a 1, 2 or 3 to fix the direction of the bubble function
                bFset->addBasisFunction(new Base::Basis_Curl_interior(l, m, n, 1));
                bFset->addBasisFunction(new Base::Basis_Curl_interior(l, m, n, 2));
                bFset->addBasisFunction(new Base::Basis_Curl_interior(l, m, n, 3));
                //	    std::cout<<"constructed interior functions with l="<<l<<" and m="<<m<<" and n="<<n<<std::endl;
            }
        }
    }
    setDefaultBasisFunctionSet(bFset);
}
