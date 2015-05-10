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
#ifndef GaussQuadratureRulesForTetrahedron_h
#define GaussQuadratureRulesForTetrahedron_h
//---------------------------------------------------------------------------
#include "Integration/QuadratureRules/GaussQuadratureRule.h"
#include <vector>

//---------------------------------------------------------------------------
namespace QuadratureRules
{
    using Geometry::PointReference;
    using Geometry::ReferenceGeometry;
    
//---------------------------------------------------------------------------
    class Tn3_1_1 : public GaussQuadratureRule
    {
    public:
        using ReferenceGeometryT = ReferenceGeometry;
        using PointReferenceT = PointReference;
    public:
        static Tn3_1_1& Instance()
        {
            static Tn3_1_1 theInstance;
            return theInstance;
        }
        
        virtual std::string getName() const;
        virtual std::size_t order() const;
        virtual std::size_t dimension() const;
        virtual std::size_t nrOfPoints() const;
        virtual double weight(std::size_t i) const;
        virtual const PointReferenceT& getPoint(std::size_t i) const;
        virtual ReferenceGeometryT* forReferenceGeometry() const;

    private:
        Tn3_1_1();
        Tn3_1_1(const Tn3_1_1&);
        virtual ~Tn3_1_1();
    private:
        const std::string name_;
        double weight_[1];
        ReferenceGeometryT* const refGeoPtr_;
        std::vector<PointReferenceT> gp_;
    };
    
//---------------------------------------------------------------------------
    class Tn3_2_1 : public GaussQuadratureRule
    {
    public:
        using ReferenceGeometryT = ReferenceGeometry;
        using PointReferenceT = PointReference;
    public:
        static Tn3_2_1& Instance()
        {
            static Tn3_2_1 theInstance;
            return theInstance;
        }
        
        virtual std::string getName() const;
        virtual std::size_t order() const;
        virtual std::size_t dimension() const;
        virtual std::size_t nrOfPoints() const;
        virtual double weight(std::size_t i) const;
        virtual const PointReferenceT& getPoint(std::size_t i) const;
        virtual ReferenceGeometryT* forReferenceGeometry() const;

    private:
        Tn3_2_1();
        Tn3_2_1(const Tn3_2_1&);
        virtual ~Tn3_2_1();
    private:
        const std::string name_;
        double weight_[4];
        ReferenceGeometryT* const refGeoPtr_;
        std::vector<PointReferenceT> gp_;
    };
    
//---------------------------------------------------------------------------
    class Tn3_3_1 : public GaussQuadratureRule
    {
    public:
        using ReferenceGeometryT = ReferenceGeometry;
        using PointReferenceT = PointReference;
    public:
        static Tn3_3_1& Instance()
        {
            static Tn3_3_1 theInstance;
            return theInstance;
        }
        
        virtual std::string getName() const;
        virtual std::size_t order() const;
        virtual std::size_t dimension() const;
        virtual std::size_t nrOfPoints() const;
        virtual double weight(std::size_t i) const;
        virtual const PointReferenceT& getPoint(std::size_t i) const;
        virtual ReferenceGeometryT* forReferenceGeometry() const;

    private:
        Tn3_3_1();
        Tn3_3_1(const Tn3_3_1&);
        virtual ~Tn3_3_1();
    private:
        const std::string name_;
        double weight_[5];
        ReferenceGeometryT* const refGeoPtr_;
        std::vector<PointReferenceT> gp_;
    };
    
//---------------------------------------------------------------------------
    class Tn3_4_1 : public GaussQuadratureRule
    {
    public:
        using ReferenceGeometryT = ReferenceGeometry;
        using PointReferenceT = PointReference;
    public:
        static Tn3_4_1& Instance()
        {
            static Tn3_4_1 theInstance;
            return theInstance;
        }
        
        virtual std::string getName() const;
        virtual std::size_t order() const;
        virtual std::size_t dimension() const;
        virtual std::size_t nrOfPoints() const;
        virtual double weight(std::size_t i) const;
        virtual const PointReferenceT& getPoint(std::size_t i) const;
        virtual ReferenceGeometryT* forReferenceGeometry() const;

    private:
        Tn3_4_1();
        Tn3_4_1(const Tn3_4_1&);
        virtual ~Tn3_4_1();
    private:
        const std::string name_;
        double weight_[11];
        ReferenceGeometryT* const refGeoPtr_;
        std::vector<PointReferenceT> gp_;
    };
    
//---------------------------------------------------------------------------
    class T3_5_1 : public GaussQuadratureRule
    {
    public:
        using ReferenceGeometryT = ReferenceGeometry;
        using PointReferenceT = PointReference;
    public:
        static T3_5_1& Instance()
        {
            static T3_5_1 theInstance;
            return theInstance;
        }
        
        virtual std::string getName() const;
        virtual std::size_t order() const;
        virtual std::size_t dimension() const;
        virtual std::size_t nrOfPoints() const;
        virtual double weight(std::size_t i) const;
        virtual const PointReferenceT& getPoint(std::size_t i) const;
        virtual ReferenceGeometryT* forReferenceGeometry() const;

    private:
        T3_5_1();
        T3_5_1(const T3_5_1&);
        virtual ~T3_5_1();
    private:
        const std::string name_;
        double weight_[14];
        ReferenceGeometryT* const refGeoPtr_;
        std::vector<PointReferenceT> gp_;
    };
    
//---------------------------------------------------------------------------
    class T3_6_1 : public GaussQuadratureRule
    {
    public:
        using ReferenceGeometryT = ReferenceGeometry;
        using PointReferenceT = PointReference;
    public:
        static T3_6_1& Instance()
        {
            static T3_6_1 theInstance;
            return theInstance;
        }
        
        virtual std::string getName() const;
        virtual std::size_t order() const;
        virtual std::size_t dimension() const;
        virtual std::size_t nrOfPoints() const;
        virtual double weight(std::size_t i) const;
        virtual const PointReferenceT& getPoint(std::size_t i) const;
        virtual ReferenceGeometryT* forReferenceGeometry() const;

    private:
        T3_6_1();
        T3_6_1(const T3_6_1&);
        virtual ~T3_6_1();
    private:
        const std::string name_;
        double weight_[24];
        ReferenceGeometryT* const refGeoPtr_;
        std::vector<PointReferenceT> gp_;
    };
    
//---------------------------------------------------------------------------
    class T3_7_1 : public GaussQuadratureRule
    {
    public:
        using ReferenceGeometryT = ReferenceGeometry;
        using PointReferenceT = PointReference;
    public:
        static T3_7_1& Instance()
        {
            static T3_7_1 theInstance;
            return theInstance;
        }
        
        virtual std::string getName() const;
        virtual std::size_t order() const;
        virtual std::size_t dimension() const;
        virtual std::size_t nrOfPoints() const;
        virtual double weight(std::size_t i) const;
        virtual const PointReferenceT& getPoint(std::size_t i) const;
        virtual ReferenceGeometryT* forReferenceGeometry() const;

    private:
        T3_7_1();
        T3_7_1(const T3_7_1&);
        virtual ~T3_7_1();
    private:
        const std::string name_;
        double weight_[31];
        ReferenceGeometryT* const refGeoPtr_;
        std::vector<PointReferenceT> gp_;
    };
    
//---------------------------------------------------------------------------
    class T3_8_1 : public GaussQuadratureRule
    {
    public:
        using ReferenceGeometryT = ReferenceGeometry;
        using PointReferenceT = PointReference;
    public:
        static T3_8_1& Instance()
        {
            static T3_8_1 theInstance;
            return theInstance;
        }
        
        virtual std::string getName() const;
        virtual std::size_t order() const;
        virtual std::size_t dimension() const;
        virtual std::size_t nrOfPoints() const;
        virtual double weight(std::size_t i) const;
        virtual const PointReferenceT& getPoint(std::size_t i) const;
        virtual ReferenceGeometryT* forReferenceGeometry() const;

    private:
        T3_8_1();
        T3_8_1(const T3_8_1&);
        virtual ~T3_8_1();
    private:
        const std::string name_;
        double weight_[43];
        ReferenceGeometryT* const refGeoPtr_;
        std::vector<PointReferenceT> gp_;
    };
    
//---------------------------------------------------------------------------
    class T3_9_1 : public GaussQuadratureRule
    {
    public:
        using ReferenceGeometryT = ReferenceGeometry;
        using PointReferenceT = PointReference;
    public:
        static T3_9_1& Instance()
        {
            static T3_9_1 theInstance;
            return theInstance;
        }
        
        virtual std::string getName() const;
        virtual std::size_t order() const;
        virtual std::size_t dimension() const;
        virtual std::size_t nrOfPoints() const;
        virtual double weight(std::size_t i) const;
        virtual const PointReferenceT& getPoint(std::size_t i) const;
        virtual ReferenceGeometryT* forReferenceGeometry() const;

    private:
        T3_9_1();
        T3_9_1(const T3_9_1&);
        virtual ~T3_9_1();
    private:
        const std::string name_;
        double weight_[53];
        ReferenceGeometryT* const refGeoPtr_;
        std::vector<PointReferenceT> gp_;
    };
    
//---------------------------------------------------------------------------
    class T3_10_1 : public GaussQuadratureRule
    {
    public:
        using ReferenceGeometryT = ReferenceGeometry;
        using PointReferenceT = PointReference;
    public:
        static T3_10_1& Instance()
        {
            static T3_10_1 theInstance;
            return theInstance;
        }
        
        virtual std::string getName() const;
        virtual std::size_t order() const;
        virtual std::size_t dimension() const;
        virtual std::size_t nrOfPoints() const;
        virtual double weight(std::size_t i) const;
        virtual const PointReferenceT& getPoint(std::size_t i) const;
        virtual ReferenceGeometryT* forReferenceGeometry() const;

    private:
        T3_10_1();
        T3_10_1(const T3_10_1&);
        virtual ~T3_10_1();
    private:
        const std::string name_;
        double weight_[126];
        ReferenceGeometryT* const refGeoPtr_;
        std::vector<PointReferenceT> gp_;
    };

//---------------------------------------------------------------------------
}// close namespace QuadratureRules
#endif
