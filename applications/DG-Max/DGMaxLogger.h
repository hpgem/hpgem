#ifndef HPGEM_DGMAXLOGGER_H
#define HPGEM_DGMAXLOGGER_H

#include <Logger.h>

extern Logger<HPGEM_LOGLEVEL> DGMaxLogger;

void initDGMaxLogging();
void logAll(std::function<void()> log);
bool loggingSuppressed();

#endif //HPGEM_DGMAXLOGGER_H
