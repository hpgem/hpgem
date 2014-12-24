/* 
 * File:   midPoint.hpp
 * Author: irana
 *
 * Created on December 23, 2014, 1:43 PM
 */

#ifndef MIDPOINT_HPP
#define	MIDPOINT_HPP

#include <vector>
#include "ButcherTableau.hpp"

namespace Base
{

    class MidPoint : public ButcherTableau
    {
    public:

        static MidPoint& instance()
        {
            static MidPoint theInstance;
            return theInstance;
        }

        std::size_t order() const;
        std::size_t numStages() const;
        double a(std::size_t i, std::size_t j) const;
        double b(std::size_t i) const;
        double c(std::size_t i) const;

    private:
        MidPoint();

        virtual ~ MidPoint() { }

        std::size_t order_;
        std::size_t numStages_;
        std::vector<std::vector<double> > a_;
        std::vector<double> b_;
        std::vector<double> c_;
    };
}

#endif	/* MIDPOINT_HPP */

