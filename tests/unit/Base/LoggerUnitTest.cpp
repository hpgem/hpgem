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
