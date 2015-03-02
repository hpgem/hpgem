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
//------------------------------------------------------------------------------
#ifndef GaussQuadratureRulesForTriangle_hpp
#define GaussQuadratureRulesForTriangle_hpp
//---------------------------------------------------------------------------
#include "Integration/QuadratureRules/GaussQuadratureRule.hpp"
#include <vector>

//---------------------------------------------------------------------------
namespace QuadratureRules
{
    using Geometry::PointReference;
    using Geometry::ReferenceGeometry;

//---------------------------------------------------------------------------
    class Tn2_1_1: public GaussQuadratureRule
    {
    public:
        using ReferenceGeometryT = ReferenceGeometry;
        using PointReferenceT = PointReference;
    public:
        static Tn2_1_1& Instance()
            {
                static Tn2_1_1 theInstance;
                return theInstance;
            }

        virtual std::string             getName() const;
        virtual std::size_t            order() const;
        virtual std::size_t            dimension() const;
        virtual std::size_t            nrOfPoints() const;
        virtual double                  weight(std::size_t i) const;
        virtual void                    getPoint(std::size_t i, PointReferenceT& p) const;
        virtual ReferenceGeometryT*     forReferenceGeometry() const;

    private:
        Tn2_1_1();
        Tn2_1_1(const Tn2_1_1&);
        virtual ~Tn2_1_1();
    private:
        const std::string               name_;
        double                          weight_[1];
        ReferenceGeometryT* const       refGeoPtr_;
        std::vector<PointReferenceT>                 gp_;
    };

//---------------------------------------------------------------------------
    class Tn2_2_1: public GaussQuadratureRule
    {
    public:
        using ReferenceGeometryT = ReferenceGeometry;
        using PointReferenceT = PointReference;
    public:
        static Tn2_2_1& Instance()
            {
                static Tn2_2_1 theInstance;
                return theInstance;
            }

        virtual std::string             getName() const;
        virtual std::size_t            order() const;
        virtual std::size_t            dimension() const;
        virtual std::size_t            nrOfPoints() const;
        virtual double                  weight(std::size_t i) const;
        virtual void                    getPoint(std::size_t i, PointReferenceT& p) const;
        virtual ReferenceGeometryT*     forReferenceGeometry() const;

    private:
        Tn2_2_1();
        Tn2_2_1(const Tn2_2_1&);
        virtual ~Tn2_2_1();
    private:
        const std::string               name_;
        double                          weight_[3];
        ReferenceGeometryT* const       refGeoPtr_;
        std::vector<PointReferenceT>                   gp_;
    };

//---------------------------------------------------------------------------
    class Tn2_3_1: public GaussQuadratureRule
    {
    public:
        using ReferenceGeometryT = ReferenceGeometry;
        using PointReferenceT = PointReference;
    public:
        static Tn2_3_1& Instance()
            {
                static Tn2_3_1 theInstance;
                return theInstance;
            }

        virtual std::string             getName() const;
        virtual std::size_t            order() const;
        virtual std::size_t            dimension() const;
        virtual std::size_t            nrOfPoints() const;
        virtual double                  weight(std::size_t i) const;
        virtual void                    getPoint(std::size_t i, PointReferenceT& p) const;
        virtual ReferenceGeometryT*     forReferenceGeometry() const;

    private:
        Tn2_3_1();
        Tn2_3_1(const Tn2_3_1&);
        virtual ~Tn2_3_1();
    private:
        const std::string               name_;
        double                          weight_[4];
        ReferenceGeometryT* const       refGeoPtr_;
        std::vector<PointReferenceT>                   gp_;
    };

//---------------------------------------------------------------------------
    class Tn2_4_1: public GaussQuadratureRule
    {
    public:
        using ReferenceGeometryT = ReferenceGeometry;
        using PointReferenceT = PointReference;
    public:
        static Tn2_4_1& Instance()
            {
                static Tn2_4_1 theInstance;
                return theInstance;
            }

        virtual std::string             getName() const;
        virtual std::size_t            order() const;
        virtual std::size_t            dimension() const;
        virtual std::size_t            nrOfPoints() const;
        virtual double                  weight(std::size_t i) const;
        virtual void                    getPoint(std::size_t i, PointReferenceT& p) const;
        virtual ReferenceGeometryT*     forReferenceGeometry() const;

    private:
        Tn2_4_1();
        Tn2_4_1(const Tn2_4_1&);
        virtual ~Tn2_4_1();
    private:
        const std::string               name_;
        double                          weight_[6];
        ReferenceGeometryT* const       refGeoPtr_;
        std::vector<PointReferenceT>                   gp_;
    };

//---------------------------------------------------------------------------
    class T2_5_1: public GaussQuadratureRule
    {
    public:
        using ReferenceGeometryT = ReferenceGeometry;
        using PointReferenceT = PointReference;
    public:
        static T2_5_1& Instance()
            {
                static T2_5_1 theInstance;
                return theInstance;
            }

        virtual std::string             getName() const;
        virtual std::size_t            order() const;
        virtual std::size_t            dimension() const;
        virtual std::size_t            nrOfPoints() const;
        virtual double                  weight(std::size_t i) const;
        virtual void                    getPoint(std::size_t i, PointReferenceT& p) const;
        virtual ReferenceGeometryT*     forReferenceGeometry() const;

    private:
        T2_5_1();
        T2_5_1(const T2_5_1&);
        virtual ~T2_5_1();

        const std::string               name_;
        double                          weight_[7];
        ReferenceGeometryT* const       refGeoPtr_;
        std::vector<PointReferenceT>                   gp_;
    };

//---------------------------------------------------------------------------
    class T2_6_1: public GaussQuadratureRule
    {
    public:
        using ReferenceGeometryT = ReferenceGeometry;
        using PointReferenceT = PointReference;
    public:
        static T2_6_1& Instance()
            {
                static T2_6_1 theInstance;
                return theInstance;
            }

        virtual std::string             getName() const;
        virtual std::size_t            order() const;
        virtual std::size_t            dimension() const;
        virtual std::size_t            nrOfPoints() const;
        virtual double                  weight(std::size_t i) const;
        virtual void                    getPoint(std::size_t i, PointReferenceT& p) const;
        virtual ReferenceGeometryT*     forReferenceGeometry() const;

    private:
        T2_6_1();
        T2_6_1(const T2_6_1&);
        virtual ~T2_6_1();
    private:
        const std::string               name_;
        double                          weight_[12];
        ReferenceGeometryT* const       refGeoPtr_;
        std::vector<PointReferenceT>                   gp_;
    };

//---------------------------------------------------------------------------
    class T2_7_1: public GaussQuadratureRule
    {
    public:
        using ReferenceGeometryT = ReferenceGeometry;
        using PointReferenceT = PointReference;
        public:
        static T2_7_1& Instance()
            {
                static T2_7_1 theInstance;
                return theInstance;
            }

        virtual std::string             getName() const;
        virtual std::size_t            order() const;
        virtual std::size_t            dimension() const;
        virtual std::size_t            nrOfPoints() const;
        virtual double                  weight(std::size_t i) const;
        virtual void                    getPoint(std::size_t i, PointReferenceT& p) const;
        virtual ReferenceGeometryT*     forReferenceGeometry() const;

    private:
        T2_7_1();
        T2_7_1(const T2_7_1&);
        virtual ~T2_7_1();
    private:
        const std::string               name_;
        double                          weight_[13];
        ReferenceGeometryT* const       refGeoPtr_;
        std::vector<PointReferenceT>                   gp_;
    };

//---------------------------------------------------------------------------
    class T2_8_1: public GaussQuadratureRule
    {
    public:
        using ReferenceGeometryT = ReferenceGeometry;
        using PointReferenceT = PointReference;
    public:
        static T2_8_1& Instance()
            {
                static T2_8_1 theInstance;
                return theInstance;
            }

        virtual std::string             getName() const;
        virtual std::size_t            order() const;
        virtual std::size_t            dimension() const;
        virtual std::size_t            nrOfPoints() const;
        virtual double                  weight(std::size_t i) const;
        virtual void                    getPoint(std::size_t i, PointReferenceT& p) const;
        virtual ReferenceGeometryT*     forReferenceGeometry() const;

    private:
        T2_8_1();
        T2_8_1(const T2_8_1&);
        virtual ~T2_8_1();
    private:
        const std::string               name_;
        double                          weight_[16];
        ReferenceGeometryT* const       refGeoPtr_;
        std::vector<PointReferenceT>                   gp_;
    };

//---------------------------------------------------------------------------
    class T2_9_1: public GaussQuadratureRule
    {
    public:
        using ReferenceGeometryT = ReferenceGeometry;
        using PointReferenceT = PointReference;
    public:
        static T2_9_1& Instance()
            {
                static T2_9_1 theInstance;
                return theInstance;
            }

        virtual std::string             getName() const;
        virtual std::size_t            order() const;
        virtual std::size_t            dimension() const;
        virtual std::size_t            nrOfPoints() const;
        virtual double                  weight(std::size_t i) const;
        virtual void                    getPoint(std::size_t i, PointReferenceT& p) const;
        virtual ReferenceGeometryT*     forReferenceGeometry() const;

    private:
        T2_9_1();
        T2_9_1(const T2_9_1&);
        virtual ~T2_9_1();
    private:
        const std::string               name_;
        double                          weight_[19];
        ReferenceGeometryT* const       refGeoPtr_;
        std::vector<PointReferenceT>                   gp_;
    };

//---------------------------------------------------------------------------
    class T2_10_1: public GaussQuadratureRule
    {
    public:
        using ReferenceGeometryT = ReferenceGeometry;
        using PointReferenceT = PointReference;
    public:
        static T2_10_1& Instance()
            {
                static T2_10_1 theInstance;
                return theInstance;
            }

        virtual std::string             getName() const;
        virtual std::size_t            order() const;
        virtual std::size_t            dimension() const;
        virtual std::size_t            nrOfPoints() const;
        virtual double                  weight(std::size_t i) const;
        virtual void                    getPoint(std::size_t i, PointReferenceT& p) const;
        virtual ReferenceGeometryT*     forReferenceGeometry() const;

    private:
        T2_10_1();
        T2_10_1(const T2_10_1&);
        virtual ~T2_10_1();
    private:
        const std::string               name_;
        double                          weight_[25];
        ReferenceGeometryT* const       refGeoPtr_;
        std::vector<PointReferenceT>                   gp_;
    };

//---------------------------------------------------------------------------
    class T2_11_1: public GaussQuadratureRule
    {
    public:
        using ReferenceGeometryT = ReferenceGeometry;
        using PointReferenceT = PointReference;
    public:
        static T2_11_1& Instance()
            {
                static T2_11_1 theInstance;
                return theInstance;
            }

        virtual std::string             getName() const;
        virtual std::size_t            order() const;
        virtual std::size_t            dimension() const;
        virtual std::size_t            nrOfPoints() const;
        virtual double                  weight(std::size_t i) const;
        virtual void                    getPoint(std::size_t i, PointReferenceT& p) const;
        virtual ReferenceGeometryT*     forReferenceGeometry() const;

    private:
        T2_11_1();
        T2_11_1(const T2_11_1&);
        virtual ~T2_11_1();
    private:
        const std::string               name_;
        double                          weight_[28];
        ReferenceGeometryT* const       refGeoPtr_;
        std::vector<PointReferenceT>                   gp_;
    };

//---------------------------------------------------------------------------
} // close namespace QuadratureRules
#endif
