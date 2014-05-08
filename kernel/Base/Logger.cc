#include <iostream>
#include "Logger.h"

static void printMessage(std::string module, std::string msg) {
  std::cout << "Module " << module << ": " << msg << std::endl;
}
static void printError(std::string module, std::string msg) {
  std::cerr << "Module " << module << ": " << msg << std::endl;
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
