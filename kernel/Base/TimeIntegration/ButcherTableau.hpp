#ifndef BUTCHERTABLEAU_HPP
#define	BUTCHERTABLEAU_HPP


#include <vector>

namespace Base
{
    /** \brief Basic container for Butcher's tableaus.
     *
     *  Butcher tableaus are a tool to write down some numerical methods for solving 
     *  ordinary differential equations compactly.
     *  Every tableau has a matrix a, vector b and vector c containing its coefficients.
     *  Furthermore, the order and the number of stages of the numerical method
     *  corresponding to the Butcher's tableau is stored.
     */
    class ButcherTableau
    {
    public:
        virtual std::size_t order() const = 0;
        virtual std::size_t numStages() const = 0;
        virtual double a(std::size_t i, std::size_t j) const = 0;
        virtual double b(std::size_t i) const = 0;
        virtual double c(std::size_t i) const = 0;
        virtual ~ButcherTableau() {}
    };
}

#endif	/* BUTCHERTABLEAU_HPP */

