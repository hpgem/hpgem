
#include "DGMaxLogger.h"

#ifdef HPGEM_USE_MPI
#include <mpi.h>
#endif

Logger<HPGEM_LOGLEVEL> DGMaxLogger("DGMax");

static int rank = 0;
static bool logAllEnabled = false;

static void printMessage(std::string module, std::string msg)
{
    std::cout << "Proc " << rank << " " << module << ": " << msg << std::endl;
}

static void printMessage0(std::string module, std::string msg)
{
    if (rank == 0 || logAllEnabled)
    {
        printMessage(std::move(module), std::move(msg));
    }
}

void logAll(std::function<void()> log)
{
    logAllEnabled = true;
    log();
    logAllEnabled = false;
}

bool loggingSuppressed()
{
    return rank != 0 && !logAllEnabled;
}

void initDGMaxLogging()
{
#ifdef HPGEM_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
    static LoggerOutput output = {
            // Keep the error as original, it prints stacktraces.
            loggerOutput->onError,
            printMessage,  // Warnings are important for each processor
            printMessage0, // Info, Verbose and Debug only on the first processor
            printMessage0,
            printMessage0
    };
    loggerOutput = &output;
}
