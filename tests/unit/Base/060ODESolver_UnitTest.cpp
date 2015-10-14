//test if the ODE solvers give the correct order
//with dy/dx = y, solution y = exp(x)

#include <vector>
#include <iostream>
#include <cmath>
#include "Base/TimeIntegration/AllTimeIntegrators.h"
#include "Logger.h"

double executeOneTimeStep(const TimeIntegration::ButcherTableau *integrator, double u, double dt)
{
    std::vector<double> k;
    //iterate over the stages of the Runge Kutta method, compute temporary solutions
    for (std::size_t level = 0; level < integrator->getNumStages(); ++level)
    {
        double kNew = u;
        for (std::size_t i = 0; i < level; ++i)
        {
            kNew += integrator->getA(level, i) * k[i] * dt;
        }
        k.push_back(kNew);
    }
    
    //Combine all temporary solutions to the solution for the next time step        
    double newVal = u;
    for (std::size_t i = 0; i < integrator->getNumStages(); ++i)
    {
        newVal += dt * (integrator->getB(i)) * k[i];
    }
    return newVal;
}

int main()
{
    const TimeIntegration::ButcherTableau *integrator1 = TimeIntegration::AllTimeIntegrators::Instance().getRule(1, 1);
    const TimeIntegration::ButcherTableau *integrator2 = TimeIntegration::AllTimeIntegrators::Instance().getRule(2, 2);
    const TimeIntegration::ButcherTableau *integrator4 = TimeIntegration::AllTimeIntegrators::Instance().getRule(4, 4);
    const TimeIntegration::ButcherTableau *integrator1Tvd = TimeIntegration::AllTimeIntegrators::Instance().getRule(1, 1, true);
    const TimeIntegration::ButcherTableau *integrator2Tvd = TimeIntegration::AllTimeIntegrators::Instance().getRule(2, 2, true);
    const TimeIntegration::ButcherTableau *integrator3Tvd = TimeIntegration::AllTimeIntegrators::Instance().getRule(3, 3, true);
    
    double dt = 0.04;
    double u1 = 1;
    double u2 = 1;
    double u4 = 1;
    double u1Tvd = 1;
    double u2Tvd = 1;
    double u3Tvd = 1;
    for (double t = 0; t < 1; t += dt)
    {
        u1 = executeOneTimeStep(integrator1, u1, dt);
        u2 = executeOneTimeStep(integrator2, u2, dt);
        u4 = executeOneTimeStep(integrator4, u4, dt);
        u1Tvd = executeOneTimeStep(integrator1Tvd, u1Tvd, dt);
        u2Tvd = executeOneTimeStep(integrator2Tvd, u2Tvd, dt);
        u3Tvd = executeOneTimeStep(integrator3Tvd, u3Tvd, dt);
    }
    
    double error11 = std::abs(std::exp(1) - u1);
    double error21 = std::abs(std::exp(1) - u2);
    double error41 = std::abs(std::exp(1) - u4);
    double error1Tvd1 = std::abs(std::exp(1) - u1Tvd);
    double error2Tvd1 = std::abs(std::exp(1) - u2Tvd);
    double error3Tvd1 = std::abs(std::exp(1) - u3Tvd);
    
    dt /= 2;
    u1 = 1;
    u2 = 1;
    u4 = 1;
    u1Tvd = 1;
    u2Tvd = 1;
    u3Tvd = 1;
    for (double t = 0; t < 1; t += dt)
    {
        u1 = executeOneTimeStep(integrator1, u1, dt);
        u2 = executeOneTimeStep(integrator2, u2, dt);
        u4 = executeOneTimeStep(integrator4, u4, dt);
        u1Tvd = executeOneTimeStep(integrator1, u1Tvd, dt);
        u2Tvd = executeOneTimeStep(integrator2Tvd, u2Tvd, dt);
        u3Tvd = executeOneTimeStep(integrator3Tvd, u3Tvd, dt);
    }
    
    double error12 = std::abs(std::exp(1) - u1);
    double error22 = std::abs(std::exp(1) - u2);
    double error42 = std::abs(std::exp(1) - u4);
    double error1Tvd2 = std::abs(std::exp(1) - u1Tvd);
    double error2Tvd2 = std::abs(std::exp(1) - u2Tvd);
    double error3Tvd2 = std::abs(std::exp(1) - u3Tvd);
    
    dt /= 2;
    u1 = 1;
    u2 = 1;
    u4 = 1;
    u1Tvd = 1;
    u2Tvd = 1;
    u3Tvd = 1;
    for (double t = 0; t < 1; t += dt)
    {
        u1 = executeOneTimeStep(integrator1, u1, dt);
        u2 = executeOneTimeStep(integrator2, u2, dt);
        u4 = executeOneTimeStep(integrator4, u4, dt);
        u1Tvd = executeOneTimeStep(integrator1Tvd, u1Tvd, dt);
        u2Tvd = executeOneTimeStep(integrator2Tvd, u2Tvd, dt);
        u3Tvd = executeOneTimeStep(integrator3Tvd, u3Tvd, dt);
    }
    
    double error13 = std::abs(std::exp(1) - u1);
    double error23 = std::abs(std::exp(1) - u2);
    double error43 = std::abs(std::exp(1) - u4);
    double error1Tvd3 = std::abs(std::exp(1) - u1Tvd);
    double error2Tvd3 = std::abs(std::exp(1) - u2Tvd);
    double error3Tvd3 = std::abs(std::exp(1) - u3Tvd);
    
    //the first order methods should return error1/error2 = 2
    logger.assert_always(error11 / error12 > 1.5 && error11 / error12 < 2.5, "First order");
    logger.assert_always(error12 / error13 > 1.5 && error12 / error13 < 2.5, "First order");
    
    logger.assert_always(error1Tvd1 / error1Tvd2 > 1.5 && error1Tvd1 / error1Tvd2 < 2.5, "First order");
    logger.assert_always(error1Tvd2 / error1Tvd3 > 1.5 && error1Tvd2 / error1Tvd3 < 2.5, "First order");

    //the second order methods should return error1/error2 = 4
    logger.assert_always(error21 / error22 > 3.5 && error21 / error22 < 4.5, "Second order");
    logger.assert_always(error22 / error23 > 3.5 && error22 / error23 < 4.5, "Second order");
    
    logger.assert_always(error2Tvd1 / error2Tvd2 > 3.5 && error2Tvd1 / error2Tvd2 < 4.5, "Second order");
    logger.assert_always(error2Tvd2 / error2Tvd3 > 3.5 && error2Tvd2 / error2Tvd3 < 4.5, "Second order");

    //the third order methods should return error1/error2 = 8
    logger.assert_always(error3Tvd1 / error3Tvd2 > 7.5 && error3Tvd1 / error3Tvd2 < 8.5, "Second order");
    logger.assert_always(error3Tvd2 / error3Tvd3 > 7.5 && error3Tvd2 / error3Tvd3 < 8.5, "Second order");

    //the fourth order methods should return error1/error2 = 16
    logger.assert_always(error41 / error42 > 15.5 && error41 / error42 < 16.5, "Fourth order");
    logger.assert_always(error42 / error43 > 15.5 && error42 / error43 < 16.5, "Fourth order");
    return 0;
}
