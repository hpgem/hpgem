
#include "DGMaxLogger.h"
#include <iostream>
#include <iomanip>

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

void initDGMaxLogging() {
    time(&initTime);
#ifdef HPGEM_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &commRank);
    int commSize;
    MPI_Comm_size(MPI_COMM_WORLD, &commSize);
    while (commSize > 1) {
        log10CommSize++;
        commSize /= 10;
    }
#endif
    static LoggerOutput output = {
        // Keep the error as original, it prints stacktraces.
        loggerOutput->onError,
        printMessage,   // Warnings are important for each processor
        printMessage0,  // Info, Verbose and Debug only on the first processor
        printMessage0, printMessage0};
    loggerOutput = &output;
}
