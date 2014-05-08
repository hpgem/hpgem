/*
 This file forms part of hpGEM. This package has been developed over a number of years by various people at the University of Twente and a full list of contributors can be found at
 http://hpgem.org/about-the-code/team
 
 This code is distributed using BSD 3-Clause License. A copy of which can found below.
 
 
 Copyright (c) 2014, Univesity of Twenete
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 
 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <iostream>
#include "Base/Logger.h"

// --- Declaring a logger.
// --- This allows you to redefine LogLevels based on command line options.
#ifndef LOG_MAIN_LEVEL
#define LOG_MAIN_LEVEL Log::FATAL
#endif
Logger<LOG_MAIN_LEVEL> log("Main");

void logMessage(std::string, std::string);

int main(int argc, char** argv) {
    
    int x = 3;
    //Basic use cases
    log.log(Log::ERROR, "Oopsie!");
    log.log(Log::FATAL, "Mweh. x = %", x);
    log.log(Log::DEBUG, "You won't see me!");
    log.log(Log::WARN, "Escapes are possible! %\% sure!", 100.01f);
    
    //Usage case for redefining with an function
    loggerOutput->onWarn = logMessage;
    log.log(Log::WARN, "Custom logger! % + % = %, %!",
            3, 5, 3+5, "yay");
    
    //Usage case for redefining with a lambda func
    loggerOutput->onFatal = [](std::string module, std::string message) {
        std::cerr << "A fatal error has occurred."
        << "\n  Module: " << module
        << "\n  Message: " << message << std::endl;
//        std::exit(-1);
    };
    
    log.log(Log::FATAL, "Null pointer passed!");
    std::cout << "You shouldn't see me." << std::endl;
    return 0;
}

void logMessage(std::string module, std::string msg) {
    std::cout << "Custom logger. Module " << module << ": " << msg << std::endl;
}
