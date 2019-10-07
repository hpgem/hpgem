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

#include "ElementInfos.h"
#include "Geometry/PointReference.h"
#include "Geometry/ReferenceGeometry.h"
#include "Geometry/PointPhysical.h"

template<std::size_t DIM>
double computeEpsilon(const Base::Element& element, std::size_t structureType);



ElementInfos::ElementInfos(double epsilon)
    : epsilon_ (epsilon)
{}

template<std::size_t DIM>
ElementInfos* ElementInfos::createStructure(
        const Base::Element &element,
        std::function<double(const Geometry::PointPhysical<DIM> &)> epsilon)
{
    logger.assert_debug(element.getReferenceGeometry()->getDimension() == DIM,
            "Incorrect dimension");
    const Geometry::PointReference<DIM>& p = element.getReferenceGeometry()->getCenter();
    Geometry::PointPhysical<DIM> pPhys = element.referenceToPhysical(p);
    return new ElementInfos(epsilon(pPhys));
}

template<std::size_t DIM>
double jelmerStructure(const Geometry::PointPhysical<DIM>& pPhys, std::size_t structureType)
{
    if(structureType == 0)
    {// Vacuum Case
        return 1;
    }
    else if (structureType == 1)
    {//Bragg Stack
        if(pPhys[0]<0.5)
        {
            return 13;
        }
        else
        {
            return 1;
        }
    }
    else if (structureType == 2)
    {//Cylinder Case with radius 0.2a
        if((pPhys[0]-0.5)*(pPhys[0]-0.5) + (pPhys[1]-0.5)*(pPhys[1]-0.5) <= 0.2 * 0.2)
        {
            return 13;
            //std::cout << pPhys[0] << " " << pPhys[1] << " " << pPhys[2]  << " " << epsilon_ << "\n";
        }
        else
        {
            return 1;
            //std::cout << pPhys[0] << " " << pPhys[1] << " " << pPhys[2]  << " " << epsilon_ << "\n";
        }
    }
    else if (structureType == 3)
    {//Cube in Cuboid Case with width of pilars of 0.1a
        if(pPhys[0] < 0.1 || pPhys[0] > 0.9 || pPhys[1] < 0.1 || pPhys[1] > 0.9)
        {
            return 1;
            //std::cout << pPhys[0] << " " << pPhys[1] << " " << pPhys[2]  << " " << epsilon_ << "\n";

        }
        else
        {
            return 13;
            //std::cout << pPhys[0] << " " << pPhys[1] << " " << pPhys[2]  << " " << epsilon_ << "\n";

        }
    }

    else if (structureType == 4)
    {//Inverse Woodpile
        //here y has length 10, whereas the lenght of x and y are length(y)/sqrt(2) = 7.07.
        // The first three circles, 2 halves and 1 whole circle, are in x,y plane.
        // The second set of circles, 4 quarters and 1 whole, is in the y,z plane. Diameter of cylinders are 0.19a.
        /*
        if((pPhys[0]-0.3535)*(pPhys[0]-0.3535) + (pPhys[1]-0.75)*(pPhys[1]-0.75) <= 0.19*0.19    ||
           (pPhys[0]-0.707)*(pPhys[0]-0.707) + (pPhys[1]-0.25)*(pPhys[1]-0.25) <= 0.19*0.19      ||
           (pPhys[0])*(pPhys[0]) + (pPhys[1]-0.25)*(pPhys[1]-0.25) <= 0.19*0.19                  ||

           (pPhys[2]-0.3535)*(pPhys[2]-0.3535) + (pPhys[1]-0.5)*(pPhys[1]-0.5) <= 0.19*0.19      ||
           (pPhys[2])*(pPhys[2]) + (pPhys[1])*(pPhys[1]) <= 0.19*0.19                            ||
           (pPhys[2])*(pPhys[2]) + (pPhys[1]-1)*(pPhys[1]-1) <= 0.19*0.19                        ||
           (pPhys[2]-0.707)*(pPhys[2]-0.707) + (pPhys[1])*(pPhys[1])                             ||
           (pPhys[2]-0.707)*(pPhys[2]-0.707) + (pPhys[1]-1)*(pPhys[1]-1) <= 0.19*0.19)
        {
            epsilon_ = 13;
            std::cout << pPhys[0] << " " << pPhys[1] << " " << pPhys[2]  << " " << epsilon_ << "\n";
        }
        else
        {
            epsilon_ = 1;
            std::cout << pPhys[0] << " " << pPhys[1] << " " << pPhys[2]  << " " << epsilon_ << "\n";
        }
         */

        if((pPhys[2]-0.3535)*(pPhys[2]-0.3535) + (pPhys[1]-0.5)*(pPhys[1]-0.5) <= 0.19*0.19       ||
           (pPhys[2]-0.707)*(pPhys[2]-0.707) + (pPhys[1])*(pPhys[1]) <= 0.19*0.19                    ||
           (pPhys[2]-0.707)*(pPhys[2]-0.707) + (pPhys[1]-1)*(pPhys[1]-1) <= 0.19*0.19                ||
           (pPhys[2])*(pPhys[2]) + (pPhys[1])*(pPhys[1]) <= 0.19*0.19                                ||
           (pPhys[2])*(pPhys[2]) + (pPhys[1]-1)*(pPhys[1]-1) <= 0.19*0.19                            ||

           (pPhys[0]-0.3535)*(pPhys[0]-0.3535) + (pPhys[1]-0.75)*(pPhys[1]-0.75) <= 0.19*0.19        ||
           (pPhys[0])*(pPhys[0]) + (pPhys[1]-0.25)*(pPhys[1]-0.25) <= 0.19*0.19                      ||
           (pPhys[0]-0.707)*(pPhys[0]-0.707) + (pPhys[1]-0.25)*(pPhys[1]-0.25) <= 0.19*0.19)
        {
            return 1;
            //std::cout << pPhys[0] << " " << pPhys[1] << " " << pPhys[2]  << " " << epsilon_ << "\n";
        }
        else
        {
            return 13;
            //std::cout << pPhys[0] << " " << pPhys[1] << " " << pPhys[2]  << " " << epsilon_ << "\n";
        }

    }
    else
    {
        std::cout << "Incorrect value for SetEpsilon" << "\n";
        return 1.0;
    }
}

// Instantiate
template
ElementInfos* ElementInfos::createStructure(const Base::Element&,
        std::function<double(const Geometry::PointPhysical<2> &)> epsilon);
template
ElementInfos* ElementInfos::createStructure(
        const Base::Element&,
        std::function<double(const Geometry::PointPhysical<3> &)> epsilon);

template
double jelmerStructure(const Geometry::PointPhysical<2>&, std::size_t);
template
double jelmerStructure(const Geometry::PointPhysical<3>&, std::size_t);
