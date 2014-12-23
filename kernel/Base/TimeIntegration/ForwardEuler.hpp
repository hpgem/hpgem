/* 
 * File:   forwardEuler.hpp
 * Author: irana
 *
 * Created on December 23, 2014, 1:31 PM
 */

#ifndef FORWARDEULER_HPP
#define	FORWARDEULER_HPP

#include <vector>
#include "ButcherTableau.hpp"

namespace Base
{

    class ForwardEuler : public ButcherTableau
    {
    public:

        static ForwardEuler& instance()
        {
            static ForwardEuler theInstance;
            return theInstance;
        }

        std::size_t order() const;
        std::size_t numStages() const;
        double a(std::size_t i, std::size_t j) const;
        double b(std::size_t i) const;
        double c(std::size_t i) const;

    private:
        ForwardEuler();

        virtual ~ ForwardEuler() { }

        std::size_t order_;
        std::size_t numStages_;
        std::vector<std::vector<double> > a_;
        std::vector<double> b_;
        std::vector<double> c_;
    };
}

#endif	/* FORWARDEULER_HPP */