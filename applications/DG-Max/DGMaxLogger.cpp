/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2018, University of Twente
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
#include "DGMaxLogger.h"
#include "Base/CommandLineOptions.h"
#include <iostream>
#include <iomanip>
#include <algorithm>

#ifdef HPGEM_USE_MPI
#include <mpi.h>

using namespace hpgem;

#endif

Logger<HPGEM_LOGLEVEL> DGMaxLogger("DGMax");

static int commRank = 0;
// Base 10 logarithm of the commsize, rounded up. This is the number of digits
// needed.
static int log10CommSize = 1;
static bool logAllEnabled = false;
static time_t initTime;

static void printMessage(std::string module, std::string msg) {
    time_t currentTime;
    time(&currentTime);
    time_t relTime = currentTime - initTime;
    std::cout << '[' << std::setw(4) << relTime << ']';

    std::cout << "@" << std::setw(log10CommSize) << commRank;
    std::cout << " " << module << ": " << msg << std::endl;
}

static void printMessage0(std::string module, std::string msg) {
    if (commRank == 0 || logAllEnabled) {
        printMessage(std::move(module), std::move(msg));
    }
}

void logAll(std::function<void()> log) {
    logAllEnabled = true;
    log();
    logAllEnabled = false;
}

bool loggingSuppressed() { return commRank != 0 && !logAllEnabled; }

static Base::CommandLineOption<std::string>* restrictionOption = nullptr;

void registerLogLevelCommandLineFlag() {
    restrictionOption = &Base::register_argument<std::string>(
        '\0', "mpiloglevel",
        "Level from which to log messages from all MPI processes, possible "
        "values error, warn, info, verbose, debug",
        false, "");
}

namespace {

using LogFunction = std::function<void(std::string, std::string)>;

/**
 * Select the function to log a message based on the restriction level.
 *
 * @param level Level of the messages to log
 * @param restriction Level below which not to log from all processes
 * @return
 */
LogFunction restrictedLogFunction(Log level, Log restriction) {
    if (level <= restriction) {
        // Lower log levels are more important
        return printMessage;
    } else {
        return printMessage0;
    }
}

}  // namespace

void initDGMaxLogging(Log restrict, bool useCommandLine) {
    time(&initTime);
#ifdef HPGEM_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &commRank);
    int commSize;
    MPI_Comm_size(MPI_COMM_WORLD, &commSize);
    while (commSize > 1) {
        log10CommSize++;
        commSize /= 10;
    }
    if (useCommandLine && restrictionOption != nullptr &&
        restrictionOption->isUsed()) {
        if (!restrictionOption->hasArgument()) {
            logger(WARN, "Command line log level requires an argument");
        }
        // Convert to lower case for user friendlyness
        std::string value = restrictionOption->getValue();
        std::transform(value.begin(), value.end(), value.begin(),
                       [](unsigned char c) { return std::tolower(c); });
        if (value == "error") {
            restrict = Log::ERROR;
        } else if (value == "warn") {
            restrict = Log::WARN;
        } else if (value == "info") {
            restrict = Log::INFO;
        } else if (value == "verbose") {
            restrict = Log::VERBOSE;
        } else if (value == "debug") {
            restrict = Log::DEBUG;
        } else {
            logger(WARN, "Unknown value for logger restriction %",
                   restrictionOption->getValue());
        }
    }

#endif
    static LoggerOutput output = {
        // Keep the fatal and error as original, it prints stacktraces.
        loggerOutput->onFatal, loggerOutput->onError,
        restrictedLogFunction(Log::WARN, restrict),
        restrictedLogFunction(Log::INFO, restrict),
        restrictedLogFunction(Log::VERBOSE, restrict),
        restrictedLogFunction(Log::DEBUG, restrict),
        // Don't override the failure handler
        loggerOutput->onFail};
    loggerOutput = &output;
}
