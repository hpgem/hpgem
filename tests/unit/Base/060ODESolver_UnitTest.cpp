//test if the ODE solvers give the correct order
//with dy/dx = y, solution y = exp(x)

#include <vector>
#include <iostream>
#include <cmath>
#include "Base/TimeIntegration/AllTimeIntegrators.hpp"
#include "Logger.h"

double executeOneTimeStep(const Base::ButcherTableau *integrator, double u, double dt)
{
    std::vector<double> k;
    //iterate over the stages of the Runge Kutta method, compute temporary solutions
    for (std::size_t level = 0; level < integrator->numStages(); ++level)
    {
        double kNew = u;
        for (std::size_t i = 0; i < level; ++i)
        {
            kNew += integrator->a(level, i) * k[i] * dt;
        }
        k.push_back(kNew);
    }

    //Combine all temporary solutions to the solution for the next time step        
    double newVal = u;
    for (std::size_t i = 0; i < integrator->numStages(); ++i)
    {
        newVal += dt * (integrator->b(i)) * k[i];
    }
    return newVal;
}

int main()
{
    const Base::ButcherTableau *integrator1 = Base::AllTimeIntegrators::Instance().getRule(1, 1);
    const Base::ButcherTableau *integrator2 = Base::AllTimeIntegrators::Instance().getRule(2, 2);
    const Base::ButcherTableau *integrator4 = Base::AllTimeIntegrators::Instance().getRule(4, 4);
    
    
    double dt = 0.04;
    double u1 = 1;
    double u2 = 1;
    double u4 = 1;
    for (double t = 0; t < 1; t += dt)
    {
        u1 = executeOneTimeStep(integrator1, u1, dt);
        u2 = executeOneTimeStep(integrator2, u2, dt);
        u4 = executeOneTimeStep(integrator4, u4, dt);
    }
    
    double error11 = std::abs(std::exp(1) - u1);
    double error21 = std::abs(std::exp(1) - u2);
    double error41 = std::abs(std::exp(1) - u4);
    
    dt /= 2;
    u1 = 1;
    u2 = 1;
    u4 = 1;
    for (double t = 0; t < 1; t += dt)
    {
        u1 = executeOneTimeStep(integrator1, u1, dt);
        u2 = executeOneTimeStep(integrator2, u2, dt);
        u4 = executeOneTimeStep(integrator4, u4, dt);
    }
    
    double error12 = std::abs(std::exp(1) - u1);
    double error22 = std::abs(std::exp(1) - u2);
    double error42 = std::abs(std::exp(1) - u4);
    
    dt /= 2;
    u1 = 1;
    u2 = 1;
    u4 = 1;
    for (double t = 0; t < 1; t += dt)
    {
        u1 = executeOneTimeStep(integrator1, u1, dt);
        u2 = executeOneTimeStep(integrator2, u2, dt);
        u4 = executeOneTimeStep(integrator4, u4, dt);
    }
    
    double error13 = std::abs(std::exp(1) - u1);
    double error23 = std::abs(std::exp(1) - u2);
    double error43 = std::abs(std::exp(1) - u4);
    
    //the first order method should return error1/error2 = 2
    logger.assert_always(error11/error12 > 1.5 && error11/error12 < 2.5, "First order");
    logger.assert_always(error12/error13 > 1.5 && error12/error13 < 2.5, "First order");
    
    //the second order method should return error1/error2 = 4
    logger.assert_always(error21/error22 > 3.5 && error21/error22 < 4.5, "Second order");    
    logger.assert_always(error22/error23 > 3.5 && error22/error23 < 4.5, "Second order");
    
    //the fourth order method should return error1/error2 = 16
    logger.assert_always(error41/error42 > 15.5 && error41/error42 < 16.5, "Fourth order");    
    logger.assert_always(error42/error43 > 15.5 && error42/error43 < 16.5, "Fourth order");
    return 0;
}