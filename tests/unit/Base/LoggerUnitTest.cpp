/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2014, University of Twente
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

 1. Redistributions of source code must retain the above copyright notice, this
 list of conditions and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.

 3. Neither the name of the copyright holder nor the names of its contributors
 may be used to endorse or promote products derived from this software without
 specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <iostream>
#include "Logger.h"

#define CATCH_CONFIG_MAIN
#include "../catch.hpp"

using namespace hpgem;
// --- Declaring a logger.
// --- This allows you to redefine LogLevels based on command line options.
#ifndef LOG_TESTING_LEVEL
#define LOG_TESTING_LEVEL Log::DEBUG
#endif
Logger<LOG_TESTING_LEVEL, false> logger2("Main");

void logMessage(std::string, std::string);

void logTestMessage(std::string, std::string);

TEST_CASE("LoggerUnitTest", "[LoggerUnitTest]") {

    // Basic use cases

    std::size_t x = 3;
    x = 4;  // suppress unused variable
    // the following two lines will abort the program (and fail the test in the
    // process)
    //    logger(ERROR, "Oopsie!");
    //    logger(FATAL, "x is not supposed to be %!!!", x);
    logger(DEBUG, "You won't see me!");
    logger2(DEBUG, "But you will see me!");
    logger(WARN, "Escapes are possible! %\\% sure!", 100.01f);

    // Usage case for redefining with an function
    loggerOutput->onWarn = logMessage;
    logger(WARN, "Custom logger! % + % = %, %!", 3, 5, 3 + 5, "yay");

    // for testing blatant assumptions about the state of the code
    logger.assert_debug(true, "Test %", 3);

    // normal logger will only crash the program in debug mode, use with care
    logger2.assert_debug(false, "Test %", 4);

    // for testing blatant assumptions when it is necessary to check again in
    // release mode (for example, opening files, processing command line input,
    // running unit tests, ...)
    logger.assert_always(true, "Test %", 5);

    // will test even when turned off
    logger2.assert_always(true, "Test %", 6);

    // test if the string parser works correctly (\ tends to accumulate
    // exponentially as the amount of escape layers grows) note that placing \%
    // in a string is technically implementation defined behaviour (to be
    // avoided in portable code) (§2.14.3.3 of the c++11 standard draft)
    loggerOutput->onDebug = logTestMessage;
    logger2(
        DEBUG,
        "If you don't like \\\\\\\\\\%, it is also possible to escape the \\%, "
        "by substituting the % with a % manually",
        '%', "%");

    // Usage case for redefining with a lambda func
    loggerOutput->onFatal = [](std::string module, std::string message) {
        std::cerr << "A fatal error has occurred."
                  << "\n  Module: " << module << "\n  Message: " << message
                  << "\n (This is part of the test.)\n"
                  << std::endl;
    };

    logger2(FATAL, "Null pointer passed!");
    std::cout << "In a normal application you wouldn't see me, but someone "
                 "redefined onFatal for the purpose of this demonstration"
              << std::endl;
}

void logMessage(std::string module, std::string msg) {
    std::cout << "Custom logger. Module " << module << ": " << msg << std::endl;
}

void logTestMessage(std::string module, std::string msg) {
    logger.assert_always(
        msg ==
            "If you don't like \\\\%, it is also possible to escape the %, by "
            "substituting the % with a % manually",
        "% got processed wrong", msg);
    logger(INFO, "% got processed correctly!", msg);
}
