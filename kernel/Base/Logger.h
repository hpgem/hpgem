#ifndef LOGGER_H
#define LOGGER_H

#include <string>
#include <sstream>
#include <functional>
#include <type_traits>

#ifndef LOG_LEVEL
#define LOG_LEVEL Log::DEFAULT
#endif

/*!
 * \brief
 */
enum class Log : signed char {
  FATAL =   -20,
  ERROR =   -15,
  WARN  =   -10,
  INFO  =   -5,
  DEFAULT = 0,
  VERBOSE = 5,
  DEBUG =   10,
};

/*!
 * \brief 
 */
constexpr bool operator<=(const Log rhs, const Log lhs) {
  return ((static_cast<signed char>(rhs)) <= (static_cast<signed char>(lhs)));
}
 
/*!
 * \class LoggerOutput
 * \brief
 */
class LoggerOutput {
  public:
    std::function<void(std::string,std::string)> onFatal;
    std::function<void(std::string,std::string)> onError;
    std::function<void(std::string,std::string)> onWarn;
    std::function<void(std::string,std::string)> onInfo;
    std::function<void(std::string,std::string)> onVerbose;
    std::function<void(std::string,std::string)> onDebug;
};

/*!
* \brief 
*/
extern LoggerOutput* loggerOutput;

template<Log L = Log::DEFAULT> class Logger;

/*!
 * \class Logger
 * \brief 
 */

template<Log Level>
class LL {
public:
  
};

extern LL<Log::FATAL>     FATAL;
extern LL<Log::ERROR>     ERROR;
extern LL<Log::WARN>      WARN;
extern LL<Log::INFO>      INFO;
extern LL<Log::DEFAULT>   DEFAULT;
extern LL<Log::VERBOSE>   VERBOSE;
extern LL<Log::DEBUG>     DEBUG;

template<Log L>
//template<typename LL>
class Logger {
  private:
      
    /*!
     * \brief 
     */
    const std::string module;
    
  public:

    /*!
     * \brief constructor
     */
     Logger(const std::string name) : module(name) { }
    /*!
     * \brief destructor
     */
     ~Logger() {}

     template<Log LOGLEVEL, typename... Args>
     typename std::enable_if<!(L < LOGLEVEL), void>::type
     operator()(const LL<LOGLEVEL> log, const std::string& format, Args&&... arg) {
         std::stringstream msgstream;
         createMessage(msgstream, format.c_str(), arg...);
         if        (LOGLEVEL <= Log::FATAL) {
           loggerOutput->onFatal(module, msgstream.str());
         } else if (LOGLEVEL <= Log::ERROR) {
           loggerOutput->onError(module, msgstream.str());
         } else if (LOGLEVEL <= Log::WARN) {
           loggerOutput->onWarn(module, msgstream.str());
         } else if (LOGLEVEL <= Log::INFO) {
           loggerOutput->onInfo(module, msgstream.str());
         } else if (LOGLEVEL <= Log::VERBOSE) {
           loggerOutput->onVerbose(module, msgstream.str());
         } else {
           loggerOutput->onDebug(module, msgstream.str());
         }
     }
     
     template<Log LOGLEVEL, typename... Args>
     typename std::enable_if<L < LOGLEVEL, void>::type
     operator()(const LL<LOGLEVEL> log, const std::string& format, Args&&... arg) {
         
     }
     
    
    /*!
     * \brief 
     */     
     template<typename... Args>
     void log(const Log loglevel, const std::string& format, Args&&... arg) {
       if (loglevel <= L || loglevel <= LOG_LEVEL) {
         std::stringstream msgstream;
         createMessage(msgstream, format.c_str(), arg...);
         if        (loglevel <= Log::FATAL) {
           loggerOutput->onFatal(module, msgstream.str());
         } else if (loglevel <= Log::ERROR) {
           loggerOutput->onError(module, msgstream.str());
         } else if (loglevel <= Log::WARN) {
           loggerOutput->onWarn(module, msgstream.str());
         } else if (loglevel <= Log::INFO) {
           loggerOutput->onInfo(module, msgstream.str());
         } else if (loglevel <= Log::VERBOSE) {
           loggerOutput->onVerbose(module, msgstream.str());
         } else {
           loggerOutput->onDebug(module, msgstream.str());
         }
       }
     }
  private:
    /*!
     * \brief 
     */
     //Put all the messages 
     template<typename Arg1, typename... Args>
     void createMessage(std::stringstream& msg, const char* fmt, 
                                                  Arg1&& arg, Args&&... args) {
       bool doSkipNext = false;
       while (*fmt != '%' && !doSkipNext) {
         doSkipNext = false;
         //Make sure we're not running past the end of our formatting string.
         if (*fmt == '\0')
           return;
         if (*fmt == '\\') { //Escape for the %sign
           doSkipNext = true;
         } else {
           msg << *fmt;
           fmt++;
         }
       }

       fmt++; //Consume the % sign
       msg << arg;
       createMessage(msg, fmt, args...);
     }
     
    /*!
     * \brief 
     */
     template<typename Arg1>
     void createMessage(std::stringstream& msg, const char* fmt, Arg1&& arg) {
       bool doSkipNext = false;
       while (*fmt != '%' && !doSkipNext) {
         doSkipNext = false;
         if (*fmt == '\0')
           return;
         if (*fmt == '\\') { //Escape for the %sign
           doSkipNext = true;
         } else {
           msg << *fmt;
           fmt++;
         }
       }
       fmt++; //Consume the % sign
       msg << arg;
       while (*fmt != '\0') {
         msg << *fmt;
         fmt++;
       }
     }
     

    /*!
     * \brief 
     */
     //No arg case.
     void createMessage(std::stringstream& msg, const char* message) {
       msg << message;
     }
};


#endif
