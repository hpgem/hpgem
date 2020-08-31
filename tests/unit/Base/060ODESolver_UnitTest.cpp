// test if the ODE solvers give the correct order
// with dy/dx = y, solution y = exp(x)

#include <vector>
#include <iostream>
#include <cmath>
#include <limits>
#include "Base/TimeIntegration/AllTimeIntegrators.h"
#include "Logger.h"

#include "../catch.hpp"

using namespace hpgem;
double executeOneTimeStep(const TimeIntegration::ButcherTableau *integrator,
                          double u, double dt) {
    std::vector<double> k;
    // iterate over the stages of the Runge Kutta method, compute temporary
    // solutions
    for (std::size_t level = 0; level < integrator->getNumberOfStages();
         ++level) {
        double kNew = u;
        for (std::size_t i = 0; i < level; ++i) {
            kNew += integrator->getA(level, i) * k[i] * dt;
        }
        k.push_back(kNew);
    }

    // Combine all temporary solutions to the solution for the next time step
    double newVal = u;
    for (std::size_t i = 0; i < integrator->getNumberOfStages(); ++i) {
        newVal += dt * (integrator->getB(i)) * k[i];
    }
    return newVal;
}

double executeOneTimeStep(const TimeIntegration::ButcherTableau *integrator,
                          double u, double maximumRelativeError, double &dt,
                          double dtMax) {
    // estimate from previous time step
    static double dtEstimate = dt;
    dt = dtEstimate;
    if (dt > dtMax) {
        dt = dtMax;
    }
    INFO(
        "Cannot use dynamic time stepping if the butcher "
        "tableau has no error estimator");
    CHECK(integrator->hasErrorEstimate());

    double currentError = std::numeric_limits<double>::infinity();
    double result = u;
    // cancel first loop
    dt *= 2.;
    while (maximumRelativeError * std::abs(result) < currentError ||
           std::isinf(currentError)) {
        dt /= 2.;
        std::vector<double> k;
        // iterate over the stages of the Runge Kutta method, compute temporary
        // solutions
        for (std::size_t level = 0; level < integrator->getNumberOfStages();
             ++level) {
            double kNew = u;
            for (std::size_t i = 0; i < level; ++i) {
                kNew += integrator->getA(level, i) * k[i] * dt;
            }
            k.push_back(kNew);
        }
        result = u;
        double difference = 0;
        for (std::size_t i = 0; i < integrator->getNumberOfStages(); ++i) {
            result += dt * (integrator->getB(i)) * k[i];
            difference += dt * integrator->getErrorCoefficient(i) * k[i];
        }
        currentError = std::abs(difference);
    }
    dtEstimate =
        0.9 * dt *
        std::pow(maximumRelativeError / currentError * std::abs(result),
                 1. / integrator->getOrder());
    logger(INFO, "dt: %; new estimate: %", dt, dtEstimate);
    return result;
}

TEST_CASE("060ODESolver_UnitTest", "[060ODESolver_UnitTest]") {
    std::vector<const TimeIntegration::ButcherTableau *> integrators{
        TimeIntegration::AllTimeIntegrators::Instance().getRule(1, 1),
        TimeIntegration::AllTimeIntegrators::Instance().getRule(2, 2),
        TimeIntegration::AllTimeIntegrators::Instance().getRule(3, 4),
        TimeIntegration::AllTimeIntegrators::Instance().getRule(4, 4),
        TimeIntegration::AllTimeIntegrators::Instance().getRule(5, 7),
        TimeIntegration::AllTimeIntegrators::Instance().getRule(1, 1, true),
        TimeIntegration::AllTimeIntegrators::Instance().getRule(2, 2, true),
        TimeIntegration::AllTimeIntegrators::Instance().getRule(3, 3, true)};

    for (auto integrator : integrators) {
        double dt = 0.04;
        double u = 1;
        for (double t = 0; t < 1; t += dt) {
            u = executeOneTimeStep(integrator, u, dt);
        }
        double error1 = std::abs(std::exp(1.) - u);
        dt /= 2.;
        u = 1;
        for (double t = 0; t < 1; t += dt) {
            u = executeOneTimeStep(integrator, u, dt);
        }
        double error2 = std::abs(std::exp(1.) - u);
        dt /= 2.;
        u = 1;
        for (double t = 0; t < 1; t += dt) {
            u = executeOneTimeStep(integrator, u, dt);
        }
        double error3 = std::abs(std::exp(1.) - u);
        double expectedRatio = std::pow(2., integrator->getOrder());
        std::string tvd =
            integrator->getTotalVariationDiminishing() ? " tvd" : "";
        INFO("" << integrator->getNumberOfStages() << " stage" << tvd
                << "time integrator of order" << integrator->getOrder()
                << " (with ratio" << error2 / error3 << ")");
        CHECK(error1 / error2 > expectedRatio * 0.875);
        CHECK(error1 / error2 < expectedRatio * 1.25);

        if (integrator->hasErrorEstimate()) {
            double dtMax = 1.;
            double maxError = std::numeric_limits<double>::infinity();
            dt = 0.04;
            u = 1;
            double t = 0;
            while (t < 1 - 1e-14) {
                double uOld = u;
                u = executeOneTimeStep(integrator, u, maxError, dt, dtMax);

                tvd = integrator->getTotalVariationDiminishing() ? " tvd" : "";
                INFO("" << integrator->getNumberOfStages() << "stage" << tvd
                        << " time integrator of order "
                        << integrator->getOrder()
                        << "exceeded error "
                           "tolerance (error was "
                        << std::abs(uOld * std::exp(dt) - u)
                        << " , but should be <" << maxError);
                CHECK(std::abs(uOld * std::exp(dt) - u) <
                      maxError * std::abs(u));

                INFO(""
                     << "stage" << tvd << " time integrator of order "
                     << integrator->getOrder()
                     << " used a larger time "
                        "step than the maximum allowed (time step was "
                     << dt << ",but should be <" << dtMax);
                CHECK(dt < dtMax + 1e-14);

                t += dt;
                dtMax = 1 - t;
            }
            maxError = 1e-10;
            dtMax = 1.;
            dt = 0.04;
            t = 0;
            u = 1;
            while (t < 1 - 1e-14) {
                double uOld = u;
                u = executeOneTimeStep(integrator, u, maxError, dt, dtMax);

                tvd = integrator->getTotalVariationDiminishing() ? " tvd" : "";
                INFO("" << integrator->getNumberOfStages() << "stage" << tvd
                        << " time integrator of order "
                        << integrator->getOrder()
                        << "exceeded error "
                           "tolerance (error was "
                        << std::abs(uOld * std::exp(dt) - u)
                        << " , but should be <" << maxError);
                CHECK(std::abs(uOld * std::exp(dt) - u) <
                      maxError * std::abs(u));

                INFO(""
                     << "stage" << tvd << " time integrator of order "
                     << integrator->getOrder()
                     << " used a larger time "
                        "step than the maximum allowed (time step was "
                     << dt << ",but should be <" << dtMax);
                CHECK(dt < dtMax + 1e-14);
                t += dt;
                dtMax = 1 - t;
            }
        }
    }
}
