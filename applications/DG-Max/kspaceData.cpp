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

#define USE_MATH_DEFINES
#include "kspaceData.h"
#include "Base/L2Norm.h"
#include "Geometry/PointPhysical.h"
#include "Logger.h"
/*
KspaceData::KspaceData(int pointsPerDirection)
        : minimumsize_(999999), current_(0)
{
    //cout<<"\n\n\n=========================================================================\n\n";
    
    //this constructor sorts the k-points for easy indexing, not efficient eigenvalue computations
    LinearAlgebra::NumericalVector k(3);
    deltak_.push_back(k);
    double stepPerPoint = M_PI / double(pointsPerDirection - 1);
    for (int i = 0; i < pointsPerDirection; ++i)
    {
        for (int j = 0; j <= i; ++j)
        {
            for (int l = 0; l <= j; ++l)
            {
                k[0] = double(i) * stepPerPoint;
                k[1] = double(j) * stepPerPoint;
                k[2] = double(l) * stepPerPoint;
                kpoints_.push_back(k);
            }
        }
    }
    
    for (int i = 1; i < kpoints_.size(); ++i)
    {
        deltak_.push_back(kpoints_[i] - kpoints_[i - 1]);
    }
    
    //create the tetrahedra per orientation
    for (int i = 0; i < pointsPerDirection - 1; ++i)
    {
        for (int j = 0; j <= i; ++j)
        {
            for (int l = 0; l <= j; ++l)
            {
                std::vector<int> localElement(4);
                localElement[0] = l + (j + 1) * j / 2 + i * (i + 1) * (i + 2) / 6;
                localElement[1] = localElement[0] + (i + 1) * (i + 2) / 2;
                localElement[2] = localElement[1] + j + 1;
                localElement[3] = localElement[2] + 1;
                elements_.push_back(localElement);
                //cout<<kpoints_[localElement[0]]<<endl<<kpoints_[localElement[1]]<<endl<<kpoints_[localElement[2]]<<endl<<kpoints_[localElement[3]]<<endl;
            }
        }
    }
    
    for (int i = 0; i < pointsPerDirection - 1; ++i)
    {
        for (int j = 0; j < i; ++j)
        {
            for (int l = 0; l <= j; ++l)
            {
                std::vector<int> localElement(4);
                localElement[0] = l + (j + 1) * j / 2 + i * (i + 1) * (i + 2) / 6;
                localElement[1] = localElement[0] + j + 1;
                localElement[2] = localElement[1] + (i + 1) * (i + 2) / 2;
                localElement[3] = localElement[2] + 1;
                elements_.push_back(localElement);
                //cout<<kpoints_[localElement[0]]<<endl<<kpoints_[localElement[1]]<<endl<<kpoints_[localElement[2]]<<endl<<kpoints_[localElement[3]]<<endl;
                localElement[2] = localElement[1] + 1;
                localElement[3] = localElement[2] + (i + 1) * (i + 2) / 2;
                elements_.push_back(localElement);
                //cout<<kpoints_[localElement[0]]<<endl<<kpoints_[localElement[1]]<<endl<<kpoints_[localElement[2]]<<endl<<kpoints_[localElement[3]]<<endl;
            }
        }
    }
    
    for (int i = 0; i < pointsPerDirection - 1; ++i)
    {
        for (int j = 0; j <= i; ++j)
        {
            for (int l = 0; l < j; ++l)
            {
                std::vector<int> localElement(4);
                localElement[0] = l + (j + 1) * j / 2 + i * (i + 1) * (i + 2) / 6;
                localElement[1] = localElement[0] + (i + 1) * (i + 2) / 2;
                localElement[2] = localElement[1] + 1;
                localElement[3] = localElement[2] + j + 1;
                elements_.push_back(localElement);
                //cout<<kpoints_[localElement[0]]<<endl<<kpoints_[localElement[1]]<<endl<<kpoints_[localElement[2]]<<endl<<kpoints_[localElement[3]]<<endl;
                localElement[1] = localElement[0] + 1;
                localElement[2] = localElement[1] + (i + 1) * (i + 2) / 2;
                localElement[3] = localElement[2] + j + 1;
                elements_.push_back(localElement);
                //cout<<kpoints_[localElement[0]]<<endl<<kpoints_[localElement[1]]<<endl<<kpoints_[localElement[2]]<<endl<<kpoints_[localElement[3]]<<endl;
            }
        }
    }
    
    for (int i = 0; i < pointsPerDirection - 1; ++i)
    {
        for (int j = 0; j < i; ++j)
        {
            for (int l = 0; l < j; ++l)
            {
                std::vector<int> localElement(4);
                localElement[0] = l + (j + 1) * j / 2 + i * (i + 1) * (i + 2) / 6;
                localElement[1] = localElement[0] + 1;
                localElement[2] = localElement[1] + j + 1;
                localElement[3] = localElement[2] + (i + 1) * (i + 2) / 2;
                elements_.push_back(localElement);
                //cout<<kpoints_[localElement[0]]<<endl<<kpoints_[localElement[1]]<<endl<<kpoints_[localElement[2]]<<endl<<kpoints_[localElement[3]]<<endl;
            }
        }
    }
}

bool KspaceData::hasNextPoint()
{
    return current_ < kpoints_.size() - 1;
}

LinearAlgebra::NumericalVector& KspaceData::nextPoint()
{
    //make sure the user is really done at this k-point
    logger.assert(omegaAtKpoints_.size() > current_, "");
    std::cout << current_++ << std::endl;
    if (functionValuesAtKpoints_.size() < current_)
    {
        LinearAlgebra::NumericalVector one(1);
        one[0] = 48;
        std::vector<LinearAlgebra::NumericalVector> data(minimumsize_, one);
        functionValuesAtKpoints_.push_back(data);
    }
    return deltak_[current_];
}

void KspaceData::setOmega(std::vector<double>& omega)
{
    //the user is allowed the freedom to think of their own intelligent ordering for the k-points so dont chack that it is consistent
    if (omega.size() < minimumsize_)
    {
        minimumsize_ = omega.size();
    }
    omegaAtKpoints_.push_back(omega);
}

void KspaceData::setFunctionValues(std::vector<LinearAlgebra::NumericalVector>& functionValues)
{
    if (functionValues.size() < minimumsize_)
    {
        minimumsize_ = functionValues.size();
    }
    functionValuesAtKpoints_.push_back(functionValues);
}

void KspaceData::getIntegral(double omega, LinearAlgebra::NumericalVector& result)
{
    logger.assert(omegaAtKpoints_.size() == kpoints_.size(), "");
    if (functionValuesAtKpoints_.size() <= current_)
    {
        LinearAlgebra::NumericalVector one(1);
        one[0] = 48;
        std::vector<LinearAlgebra::NumericalVector> data(minimumsize_, one);
        functionValuesAtKpoints_.push_back(data);
    }
    for (int i = 0; i < elements_.size(); ++i)
    {
        for (int k = 0; k < minimumsize_; ++k)
        {
            double localResult = 0;
            std::vector<LinearAlgebra::NumericalVector> kpointsWithHigherOmega, kpointsWithLowerOmega, fAtHigherPoints, fAtLowerPoints;
            std::vector<double> omegaAtHigherPoints, omegaAtLowerPoints;
            LinearAlgebra::NumericalVector deltakLocal[3], reciprocalk[3], localFunctionValue(result.size());
            double deltaOmega[3];
            
            //find the location of the omega is constant plane wrt. the corners of the tetrahedron
            for (int j = 0; j < 4; ++j)
            {
                if (omegaAtKpoints_[elements_[i][j]][k] > omega)
                {
                    kpointsWithHigherOmega.push_back(kpoints_[elements_[i][j]]);
                    omegaAtHigherPoints.push_back(omegaAtKpoints_[elements_[i][j]][k]);
                    fAtHigherPoints.push_back(functionValuesAtKpoints_[elements_[i][j]][k]);
                }
                else
                {
                    kpointsWithLowerOmega.push_back(kpoints_[elements_[i][j]]);
                    omegaAtLowerPoints.push_back(omegaAtKpoints_[elements_[i][j]][k]);
                    fAtLowerPoints.push_back(functionValuesAtKpoints_[elements_[i][j]][k]);
                }
            }
            
            if (kpointsWithHigherOmega.size() == 1)
            {
                //compute 1/|nabla_k omega|
                for (int j = 0; j < 3; ++j)
                {
                    deltaOmega[j] = omegaAtHigherPoints[0] - omegaAtLowerPoints[j];
                    deltakLocal[j] = kpointsWithHigherOmega[0] - kpointsWithLowerOmega[j];
                } //this bit is the same in all 3 cases
                OuterProduct(deltakLocal[0], deltakLocal[1], reciprocalk[2]);
                OuterProduct(deltakLocal[1], deltakLocal[2], reciprocalk[0]);
                OuterProduct(deltakLocal[2], deltakLocal[0], reciprocalk[1]);
                
                localResult = Base::L2Norm(reciprocalk[0] * deltaOmega[0] + reciprocalk[1] * deltaOmega[1] + reciprocalk[2] * deltaOmega[2]);
                localResult = fabs(reciprocalk[0][0] * deltakLocal[0][0] + reciprocalk[0][1] * deltakLocal[0][1] + reciprocalk[0][2] * deltakLocal[0][2]) / localResult;
                
                //find the three corners of the intersection
                LinearAlgebra::NumericalVector corner[3], area;
                for (int j = 0; j < 3; ++j)
                {
                    corner[j] = (omegaAtHigherPoints[0] - omega) / (omegaAtHigherPoints[0] - omegaAtLowerPoints[j]) * kpointsWithLowerOmega[j] + (omega - omegaAtLowerPoints[j]) / (omegaAtHigherPoints[0] - omegaAtLowerPoints[j]) * kpointsWithHigherOmega[0];
                }
                OuterProduct(corner[1] - corner[0], corner[2] - corner[0], area);
                
                for (int j = 0; j < 3; ++j)
                {
                    localFunctionValue += (omegaAtHigherPoints[0] - omega) / (omegaAtHigherPoints[0] - omegaAtLowerPoints[j]) * fAtLowerPoints[j] + (omega - omegaAtLowerPoints[j]) / (omegaAtHigherPoints[0] - omegaAtLowerPoints[j]) * fAtHigherPoints[0];
                }
                
                //and compute the local contribution to the DOS
                localResult *= Base::L2Norm(area) / 6;
                
                localFunctionValue *= localResult;
            }
            else if (kpointsWithHigherOmega.size() == 2)
            {
                //compute 1/|nabla_k omega|
                deltaOmega[0] = omegaAtHigherPoints[0] - omegaAtHigherPoints[1];
                deltakLocal[0] = kpointsWithHigherOmega[0] - kpointsWithHigherOmega[1];
                for (int j = 0; j < 2; ++j)
                {
                    deltaOmega[j + 1] = omegaAtHigherPoints[0] - omegaAtLowerPoints[j];
                    deltakLocal[j + 1] = kpointsWithHigherOmega[0] - kpointsWithLowerOmega[j];
                } //this bit is the same in all 3 cases
                OuterProduct(deltakLocal[0], deltakLocal[1], reciprocalk[2]);
                OuterProduct(deltakLocal[1], deltakLocal[2], reciprocalk[0]);
                OuterProduct(deltakLocal[2], deltakLocal[0], reciprocalk[1]);
                
                localResult = Base::L2Norm(reciprocalk[0] * deltaOmega[0] + reciprocalk[1] * deltaOmega[1] + reciprocalk[2] * deltaOmega[2]);
                localResult = fabs(reciprocalk[0][0] * deltakLocal[0][0] + reciprocalk[0][1] * deltakLocal[0][1] + reciprocalk[0][2] * deltakLocal[0][2]) / localResult;
                
                //find the four corners of the intersection
                LinearAlgebra::NumericalVector corner[4], area, f[4];
                for (int j = 0; j < 2; ++j)
                {
                    corner[2 * j] = (omegaAtHigherPoints[0] - omega) / (omegaAtHigherPoints[0] - omegaAtLowerPoints[j]) * kpointsWithLowerOmega[j] + (omega - omegaAtLowerPoints[j]) / (omegaAtHigherPoints[0] - omegaAtLowerPoints[j]) * kpointsWithHigherOmega[0];
                    corner[2 * j + 1] = (omegaAtHigherPoints[1] - omega) / (omegaAtHigherPoints[1] - omegaAtLowerPoints[j]) * kpointsWithLowerOmega[j] + (omega - omegaAtLowerPoints[j]) / (omegaAtHigherPoints[1] - omegaAtLowerPoints[j]) * kpointsWithHigherOmega[1];
                    
                    f[2 * j] = (omegaAtHigherPoints[0] - omega) / (omegaAtHigherPoints[0] - omegaAtLowerPoints[j]) * fAtLowerPoints[j] + (omega - omegaAtLowerPoints[j]) / (omegaAtHigherPoints[0] - omegaAtLowerPoints[j]) * fAtHigherPoints[0];
                    f[2 * j + 1] = (omegaAtHigherPoints[1] - omega) / (omegaAtHigherPoints[1] - omegaAtLowerPoints[j]) * fAtLowerPoints[j] + (omega - omegaAtLowerPoints[j]) / (omegaAtHigherPoints[1] - omegaAtLowerPoints[j]) * fAtHigherPoints[1];
                }
                
                //split the quadrilateral in triangles so the barycentre is a bit easier to find when the LDOS is computed
                OuterProduct(corner[1] - corner[0], corner[2] - corner[0], area);
                LinearAlgebra::NumericalVector temp = Base::L2Norm(area) * (f[0] + f[1] + f[2]) * localResult / 6;
                OuterProduct(corner[2] - corner[1], corner[3] - corner[1], area);
                localResult *= (Base::L2Norm(area)) / 6;
                localFunctionValue = (f[1] + f[2] + f[3]) * localResult;
                localFunctionValue += temp;
            }
            else if (kpointsWithHigherOmega.size() == 3)
            {
                //compute 1/|nabla_k omega|
                for (int j = 0; j < 3; ++j)
                {
                    deltaOmega[j] = omegaAtLowerPoints[0] - omegaAtHigherPoints[j];
                    deltakLocal[j] = kpointsWithLowerOmega[0] - kpointsWithHigherOmega[j];
                } //this bit is the same in all 3 cases
                OuterProduct(deltakLocal[0], deltakLocal[1], reciprocalk[2]);
                OuterProduct(deltakLocal[1], deltakLocal[2], reciprocalk[0]);
                OuterProduct(deltakLocal[2], deltakLocal[0], reciprocalk[1]);
                
                localResult = Base::L2Norm(reciprocalk[0] * deltaOmega[0] + reciprocalk[1] * deltaOmega[1] + reciprocalk[2] * deltaOmega[2]);
                localResult = fabs(reciprocalk[0][0] * deltakLocal[0][0] + reciprocalk[0][1] * deltakLocal[0][1] + reciprocalk[0][2] * deltakLocal[0][2]) / localResult;
                
                //find the three corners of the intersection
                LinearAlgebra::NumericalVector corner[3], area;
                for (int j = 0; j < 3; ++j)
                {
                    corner[j] = (omega - omegaAtLowerPoints[0]) / (omegaAtHigherPoints[j] - omegaAtLowerPoints[0]) * kpointsWithHigherOmega[j] + (omegaAtHigherPoints[j] - omega) / (omegaAtHigherPoints[j] - omegaAtLowerPoints[0]) * kpointsWithLowerOmega[0];
                }
                OuterProduct(corner[1] - corner[0], corner[2] - corner[0], area);
                
                for (int j = 0; j < 3; ++j)
                {
                    localFunctionValue += (omega - omegaAtLowerPoints[0]) / (omegaAtHigherPoints[j] - omegaAtLowerPoints[0]) * fAtHigherPoints[j] + (omegaAtHigherPoints[j] - omega) / (omegaAtHigherPoints[j] - omegaAtLowerPoints[0]) * fAtLowerPoints[0];
                }
                
                //and compute the local contribution to the DOS
                localResult *= Base::L2Norm(area) / 6;
                
                localFunctionValue *= localResult;
            } //else there is no intersection so contribution 0 for this element
            result += localFunctionValue;
        }
    }
    result /= double(8) * M_PI * M_PI * M_PI;
}
*/
