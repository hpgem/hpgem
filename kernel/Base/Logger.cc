#include <iostream>
#include "Logger.h"

LL<Log::FATAL>     FATAL;
LL<Log::ERROR>     ERROR;
LL<Log::WARN>      WARN;
LL<Log::INFO>      INFO;
LL<Log::DEFAULT>   DEFAULT;
LL<Log::VERBOSE>   VERBOSE;
LL<Log::DEBUG>     DEBUG;

static void printMessage(std::string module, std::string msg) {
  std::cout << "Module " << module << ": " << msg << std::endl;
}
static void printError(std::string module, std::string msg) {
  std::cerr << "Module " << module << ": " << msg << std::endl;
  std::exit(2);
}

LoggerOutput loggerOutputDefaultImpl = {
  printError, //onFatal
  printError, //onError
  printMessage, //onWarn
  printMessage, //onInfo
  printMessage, //onVerbose
  printMessage  //onDebug
};

LoggerOutput* loggerOutput = &loggerOutputDefaultImpl;
