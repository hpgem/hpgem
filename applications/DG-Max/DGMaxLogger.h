#ifndef HPGEM_APP_DGMAXLOGGER_H
#define HPGEM_APP_DGMAXLOGGER_H

#include <Logger.h>

using namespace hpgem;


extern Logger<HPGEM_LOGLEVEL> DGMaxLogger;

void initDGMaxLogging();
void logAll(std::function<void()> log);
bool loggingSuppressed();

#endif  // HPGEM_APP_DGMAXLOGGER_H
