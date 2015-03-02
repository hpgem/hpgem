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
#ifndef GaussQuadratureRulesForSquare_hpp
#define GaussQuadratureRulesForSquare_hpp
//---------------------------------------------------------------------------
#include "Integration/QuadratureRules/GaussQuadratureRule.hpp"
#include <vector>

//---------------------------------------------------------------------------
namespace QuadratureRules
{
    using Geometry::PointReference;
    using Geometry::ReferenceGeometry;
    
//---------------------------------------------------------------------------
    class Cn2_1_1:public GaussQuadratureRule
    {
    public:
        typedef ReferenceGeometry    ReferenceGeometryT;
        typedef PointReference       PointReferenceT;
        
    public:
        static Cn2_1_1& Instance()
        {
            static Cn2_1_1 theInstance;
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
        Cn2_1_1();
        Cn2_1_1(const Cn2_1_1&);
        virtual ~Cn2_1_1();

    private:
        const std::string               name_;
        double                          weight_[1];
        ReferenceGeometryT* const       refGeoPtr_;
        std::vector<PointReferenceT>                   gp_;
    };

//---------------------------------------------------------------------------
    class Cn2_3_4: public GaussQuadratureRule
    {
    public:
        typedef ReferenceGeometry    ReferenceGeometryT;
        typedef PointReference       PointReferenceT;
    public:
        static Cn2_3_4& Instance()
        {
            static Cn2_3_4 theInstance;
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
        Cn2_3_4();
        Cn2_3_4(const Cn2_3_4&);
        virtual ~Cn2_3_4();
    private:
        const std::string               name_;
        double                          weight_[4];
        ReferenceGeometryT* const       refGeoPtr_;
        std::vector<PointReferenceT>                   gp_;
    };

//---------------------------------------------------------------------------
    class Cn2_5_9: public GaussQuadratureRule
    {
    public:
        typedef ReferenceGeometry    ReferenceGeometryT;
        typedef PointReference       PointReferenceT;
    public:
        static Cn2_5_9& Instance()
            {
                static Cn2_5_9 theInstance;
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
        Cn2_5_9();
        Cn2_5_9(const Cn2_5_9&);
        virtual ~Cn2_5_9();
    private:
        const std::string               name_;
        double                          weight_[9];
        ReferenceGeometryT* const       refGeoPtr_;
        std::vector<PointReferenceT>                   gp_;
    };

//---------------------------------------------------------------------------
    class C2_7_4:public GaussQuadratureRule
    {
    public:
        typedef ReferenceGeometry    ReferenceGeometryT;
        typedef PointReference       PointReferenceT;
    public:
        static C2_7_4& Instance()
        {
            static C2_7_4 theInstance;
            return theInstance;
        }

        virtual std::string getName() const;
        
        virtual std::size_t            order() const;
        
        virtual std::size_t            dimension() const;
        
        virtual std::size_t            nrOfPoints() const;
        
        virtual double                  weight(std::size_t i) const;
        
        virtual void                    getPoint(std::size_t i, PointReferenceT& p) const;
        
        virtual ReferenceGeometryT*     forReferenceGeometry() const;

    private:
        C2_7_4();
        C2_7_4(const C2_7_4&);
        virtual ~C2_7_4();
    private:
        const std::string               name_;
        double                          weight_[16];
        ReferenceGeometryT* const       refGeoPtr_;
        std::vector<PointReferenceT>                   gp_;
    };
//---------------------------------------------------------------------------
    class C2_9_5:public GaussQuadratureRule
    {
    public:
        typedef ReferenceGeometry    ReferenceGeometryT;
        typedef PointReference       PointReferenceT;
    public:
        static C2_9_5& Instance()
        {
            static C2_9_5 theInstance;
            return theInstance;
        }

        virtual std::string getName() const;

        virtual std::size_t            order() const;

        virtual std::size_t            dimension() const;

        virtual std::size_t            nrOfPoints() const;

        virtual double                  weight(std::size_t i) const;

        virtual void                    getPoint(std::size_t i, PointReferenceT& p) const;

        virtual ReferenceGeometryT*     forReferenceGeometry() const;

    private:
        C2_9_5();
        C2_9_5(const C2_9_5&);
        virtual ~C2_9_5();
    private:
        const std::string               name_;
        double                          weight_[25];
        ReferenceGeometryT* const       refGeoPtr_;
        std::vector<PointReferenceT>                   gp_;
    };
//---------------------------------------------------------------------------
    class C2_11_6:public GaussQuadratureRule
    {
    public:
        typedef ReferenceGeometry    ReferenceGeometryT;
        typedef PointReference       PointReferenceT;
    public:
        static C2_11_6& Instance()
        {
            static C2_11_6 theInstance;
            return theInstance;
        }

        virtual std::string getName() const;

        virtual std::size_t            order() const;

        virtual std::size_t            dimension() const;

        virtual std::size_t            nrOfPoints() const;

        virtual double                  weight(std::size_t i) const;

        virtual void                    getPoint(std::size_t i, PointReferenceT& p) const;

        virtual ReferenceGeometryT*     forReferenceGeometry() const;

    private:
        C2_11_6();
        C2_11_6(const C2_11_6&);
        virtual ~C2_11_6();
    private:
        const std::string               name_;
        double                          weight_[36];
        ReferenceGeometryT* const       refGeoPtr_;
        std::vector<PointReferenceT>                   gp_;
    };
//---------------------------------------------------------------------------
} // close namespace QuadratureRules
#endif
