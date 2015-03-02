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
#ifndef GaussQuadratureRulesForPyramid_hpp
#define GaussQuadratureRulesForPyramid_hpp
//---------------------------------------------------------------------------
#include "Integration/QuadratureRules/GaussQuadratureRule.hpp"
#include <vector>

//---------------------------------------------------------------------------
namespace QuadratureRules
{

    using Geometry::PointReference;
    using Geometry::ReferenceGeometry;

//---------------------------------------------------------------------------
    class Pyramid_1_1: public GaussQuadratureRule
    {
    public:
        using PointReferenceT = PointReference;
        using ReferenceGeometryT = ReferenceGeometry;
    public:
        static Pyramid_1_1& Instance()
            {
                static Pyramid_1_1 theInstance;
                return theInstance;
            }

        virtual std::string             getName() const;
        virtual std::size_t            order() const;
        virtual std::size_t            dimension() const;
        virtual std::size_t            nrOfPoints() const;
        virtual double                  weight(std::size_t i) const;
        virtual const PointReferenceT&                    getPoint(std::size_t i) const;
        virtual ReferenceGeometryT*     forReferenceGeometry() const;

    private:
        Pyramid_1_1();
        Pyramid_1_1(const Pyramid_1_1&);
        virtual ~Pyramid_1_1();
    private:

        const std::string               name_;
        double                          weight_[4];
        ReferenceGeometryT* const       refGeoPtr_;
        std::vector<PointReferenceT>                   gp_;
    };

//---------------------------------------------------------------------------
    class Pyramid_3_1: public GaussQuadratureRule
    {
    public:
        using PointReferenceT = PointReference;
        using ReferenceGeometryT = ReferenceGeometry;
    public:
        static Pyramid_3_1& Instance()
            {
                static Pyramid_3_1 theInstance;
                return theInstance;
            }

        virtual std::string             getName() const;
        virtual std::size_t            order() const;
        virtual std::size_t            dimension() const;
        virtual std::size_t            nrOfPoints() const;
        virtual double                  weight(std::size_t i) const;
        virtual const PointReferenceT&                    getPoint(std::size_t i) const;
        virtual ReferenceGeometryT*     forReferenceGeometry() const;

    private:
        Pyramid_3_1();
        Pyramid_3_1(const Pyramid_3_1&);
        virtual ~Pyramid_3_1();
    private:
        const std::string               name_;
        double                          weight_[16];
        ReferenceGeometryT* const       refGeoPtr_;
        std::vector<PointReferenceT>                   gp_;
    };

//---------------------------------------------------------------------------
    class Pyramid_5_1: public GaussQuadratureRule
    {
     public:
        using PointReferenceT = PointReference;
        using ReferenceGeometryT = ReferenceGeometry;
    public:
        static Pyramid_5_1& Instance()
            {
                static Pyramid_5_1 theInstance;
                return theInstance;
            }

        virtual std::string             getName() const;
        virtual std::size_t            order() const;
        virtual std::size_t            dimension() const;
        virtual std::size_t            nrOfPoints() const;
        virtual double                  weight(std::size_t i) const;
        virtual const PointReferenceT&                    getPoint(std::size_t i) const;
        virtual ReferenceGeometryT*     forReferenceGeometry() const;

    private:
        Pyramid_5_1();
        Pyramid_5_1(const Pyramid_5_1&);
        virtual ~Pyramid_5_1();
    private:
        const std::string               name_;
        double                          weight_[36];
        ReferenceGeometryT* const       refGeoPtr_;
        std::vector<PointReferenceT>                   gp_;
    };

//---------------------------------------------------------------------------
    class Pyramid_7_1: public GaussQuadratureRule
    {
     public:
        using PointReferenceT = PointReference;
        using ReferenceGeometryT = ReferenceGeometry;
    public:
        static Pyramid_7_1& Instance()
            {
                static Pyramid_7_1 theInstance;
                return theInstance;
            }

        virtual std::string             getName() const;
        virtual std::size_t            order() const;
        virtual std::size_t            dimension() const;
        virtual std::size_t            nrOfPoints() const;
        virtual double                  weight(std::size_t i) const;
        virtual const PointReferenceT&                    getPoint(std::size_t i) const;
        virtual ReferenceGeometryT*     forReferenceGeometry() const;

    private:
        Pyramid_7_1();
        Pyramid_7_1(const Pyramid_7_1&);
        virtual ~Pyramid_7_1();
    private:
        const std::string               name_;
        double                          weight_[48];
        ReferenceGeometryT* const       refGeoPtr_;
        std::vector<PointReferenceT>                   gp_;
    };

//---------------------------------------------------------------------------
} // close namespace QuadratureRules
#endif
