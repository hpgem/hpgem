/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2014, University of Twente
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

#ifndef HPGEM_KERNEL_LOGGER_H
#define HPGEM_KERNEL_LOGGER_H

#include <string>
#include <sstream>
#include <functional>
#include <type_traits>

#ifndef HPGEM_LOGLEVEL
#define HPGEM_LOGLEVEL Log::DEFAULT
#endif

#ifdef HPGEM_FORCE_ASSERTS
#define HPGEM_ASSERTS true
#else
#ifdef HPGEM_NO_ASSERTS
#define HPGEM_ASSERTS false
#else
#ifdef NDEBUG
#define HPGEM_ASSERTS false
#else
#define HPGEM_ASSERTS true
#endif
#endif
#endif

/* IMPLEMENTATION DETAIL
 *    - by dducks
 *
 * A brief explanation how this class works. Beware, there is black magic
 * going on here.
 *
 * So, the previous version of the Logger used the syntax
 * Logger<Log::LEVEL> logger;
 * logger.log(Log::LEVEL, "Message", args...);
 * However, people didn't like the slightly large amount of log's in a single
 * statement. Therefore, the operator() has been chosen. The next problem I
 * personally had with the logger is the fact that it relied on a compile-time
 * resolvable if-statement which could often (always) be optimised out. The
 * problem however is that you would need optimisation enabled to have no
 * speed penalty for loglevels below visibility.
 *
 * The problem however, is that we want to generate different code based on
 * the first parameter; the loglevel. Because we actually want to generate
 * code and people don't like the preprocessor, we have to use the template
 * system to do this.
 *
 * One of the caveats of the template system is that it is only able to
 * resolve templates based on parameters, not on values. It also allows you
 * to ommit template parameters in case where the template parameters can be
 * deduced.
 *
 * Based on these two rules, we now need not a value per loglevel, but a type.
 * Therefore we use tagging, in which we use the class LL with a template
 * argument to create the different tags. Our loglevel, the enum class Log, is
 * considered to be a valid non-type template parameter - it must be either
 * integral, enum (which is integral) or a pointer - so we can use it to tag a
 * class to create diffent types.
 *
 * Now, our Logger instance has a template argument which is the loglevel
 * associated with our Logger; anything of lesser priority is ignored; It
 * also honours the symbol MERCURY_LOGLEVEL, and compares it to the lesser of
 * the template argument and the preprocessor define.
 *
 * The operator() function now should be based on the template parameter of the
 * logger class (or MERCURY_LOGLEVEL), and the loglevel of the current message.
 * Because we can't depend on the value but just on the type, we have to do some
 * magic. Now, we want the function to resolve to something that produces output
 * if the level is high enough, or nothing at all if the level is of a too low
 * priority for the current loglevel. For this we utilize std::enable_if<>
 * to select the right function.
 *
 * Because in C++ a templated function which leads to ill-formed code upon
 * instantiation is rejected and does NOT result in a compile error,
 * std::enable_if allows us to make either the version WITH output or
 * WITHOUT output, based on template parameters.
 *
 * As you can see below, the class LL is completely empty; However,
 * there are a few instances; one for every loglevel. Now, remember,
 * when using the logger, Log::DEFAULT resolves to the enum class value
 * while DEFAULT resolves to the instance of class LL with template
 * argument Log::DEFAULT.
 *
 * However, LL<Log::DEFAULT> differs in type from LL<Log::ERROR>, as the
 * template argument is different. Because of deduction, the compiler can
 * figure out the actual values for the template arguments. In case of
 * operator(), the first argument is Log loglevel, which is deduced from
 * the declaration;
 *
 * template<Log LogLevel, typename... Args>
 * void operator(const LL<LogLevel> ll, std::string text, Args... args);
 *
 * Now, the template parameter LogLevel gets filled with Log::ERROR in case
 * we give ERROR as the first argument, or Log::DEBUG with DEBUG as the first
 * argument. Now, we want to resolve to two different functions, so as the
 * return type we use std::enable_if which would lead to ill-formed code
 * in case the predicate of the std::enable_if is false. So, based on the
 * tag ll we now can select the different implementations.
 *
 * It then buffers everything into a stringstream. so as long as operator<<
 * is defined properly for your object,. This then gets redirected towards the
 * correct output channel.
 *
 * Please note that operator() is an inline, templated function. In case the
 * function body is empty, this code is very, very likely to not emit any
 * instructions at all. If you don't give arguments which have only non-const
 * functions, the function call can be considered invariant which means it
 * can completely be taken out, in the case it's a seperate function it just
 * lowers the cost.
 *
 * As said, black magic below.
 * END OF IMPLEMENTATION DETAIL
 */

/*!
 * \brief The different loglevels.
 *
 * The different loglevels, represented as signed characters,
 * in descending order of severeness. Worst is FATAL, best is DEBUG.
 *
 * Please, use the tags FATAL/ERROR/etc without class/enum/namespace instead.
 */
namespace hpgem {
enum class Log : signed char {
    FATAL = -20,
    ERROR = -15,
    WARN = -10,
    INFO = -5,
    DEFAULT = 0,
    VERBOSE = 5,
    DEBUG = 10
};

/*!
 * \brief Internally used to filter on loglevel.
 * Do not edit, as this is required for an optimised logger.
 */
constexpr bool operator<=(const Log lhs, const Log rhs) {
    return ((static_cast<signed char>(lhs)) <= (static_cast<signed char>(rhs)));
}

/*!
 * \class LoggerOutput
 * \brief Default functions for output generation
 *
 * These handlers will be called on generation of the message.
 * The functions are of signature
 *    void (std::string moduleName, std::string message);
 *
 * These functions may not return but call std::exit() instead.
 * They may also throw any exception to allow code to gracefully
 * recover.
 */
class LoggerOutput {
   public:
    std::function<void(std::string, std::string)> onFatal;
    std::function<void(std::string, std::string)> onError;
    std::function<void(std::string, std::string)> onWarn;
    std::function<void(std::string, std::string)> onInfo;
    std::function<void(std::string, std::string)> onVerbose;
    std::function<void(std::string, std::string)> onDebug;
    // Failure handler. This MUST not return, but instead either stop the
    // program or throw an appropriate exception
    std::function<void(std::string, std::string)> onFail;
};

/*!
 * \brief Declaration of the output functions.
 * If the output needs to be redirected, please
 * swap the loggerOutput pointer to your preferred
 * LoggerOutput instance, and make sure this exists
 * until _AFTER_ an std::exit() invocation.
 * (e.g. clean up with std::atexit())
 */
extern LoggerOutput* loggerOutput;

// Forward declaration..
template <Log L = Log::DEFAULT, bool ASSERTS = HPGEM_ASSERTS>
class Logger;

/*!
 * \brief Tag for template metaprogramming
 *
 * This tag class serves as a way to correctly
 * resolve the implementation (if any) of the
 * operator() or .log() method of the Logger.
 * Please, don't change it at all nor give it
 * any members.
 */
template <Log Level>
class LL {
   public:
};

/*!
 * These are the loglevels which should be used to
 * control the logger. They are declared here
 * but defined in Logger.cpp
 */
/*!
 * \brief Fatal log level
 *
 * Fatal, as in, the program has suffered from the worst possible failure and
 * there is no way it can gracefully recover.
 *
 * Example: No memory allocations possible
 *
 * Default behaviour: log to std::cerr, followed by std::exit().
 */
extern LL<Log::FATAL> FATAL;
/*!
 * \brief Error log level
 *
 * Error, as in, the program has found a severe problem which it cannot resolve
 * any further. It does not know how to recover in a sane way.
 *
 * Example: Negative timestep, Infinite end time and no override of the
 * continuation function.
 *
 * Default behaviour: log to std::cerr, followed by std::exit().
 */
extern LL<Log::ERROR> ERROR;
/*!
 * \brief Warning log level
 *
 * Warning, as in, the program has detected a problem but does know a solution.
 * The simulation can continue with this fix, but the user should look at
 * fixing his / her simulation so this won't occur in the future.
 *
 * Example: Setting a smaller Xmax than Xmin.
 *
 * Default behaviour: log to std::cerr, returns afterwards.
 */
extern LL<Log::WARN> WARN;
/*!
 * \brief Info log level
 *
 * Useful information, small oddities found which should be of no real effect
 * to the user. Also information about the current state and progress of the
 * program.
 *
 * Example: Finished inserting particles.
 *
 * Default behaviour: log to std::cout, returns afterwards.
 */
extern LL<Log::INFO> INFO;
/*!
 * \brief Default log level
 *
 * Only useful for defining the loglevel of the logger itself. Should not
 * actually be used.
 */
extern LL<Log::DEFAULT> DEFAULT;
/*!
 * \brief Verbose information
 *
 * Information which is not useful to anybody except those looking for weird
 * behaviour and statistics. These should however still be clear in meaning.
 *
 * Example: Inserted 381 particles in Insertion Boundary #1.
 *
 * Default behaviour: ignore.
 */
extern LL<Log::VERBOSE> VERBOSE;
/*!
 * \brief Debug information
 *
 * Only used for internal development. Can be very cryptic, as it is only meant
 * for finding bugs / oddities by the internal development team.
 *
 * Example: Collission found between Particle #38201 and Wall #5
 *
 * Default behaviour: ignore.
 */
extern LL<Log::DEBUG> DEBUG;

/*!
 * \brief Logger
 *
 * \arg L The log level. Messages of higher level are ignored
 *
 * Usage: logger(FATAL, "Error in (here) because % < %!\n", var1, var2);
 *
 * Define custom loggers by:
 * #ifndef HG_LOGLEVEL_CUSTOMMOD
 * #define HG_LOGLEVEL_CUSTOMMOD Log::Debug
 * #endif
 * Logger<HG_LOGLEVEL_CUSTOMMOD> customLogger;
 */
template <Log L, bool ASSERTS>
class Logger {
   private:
    /*!
     * \brief The module name of this actual logger
     */
    const std::string module;

   public:
    /*!
     * \brief constructor
     * \arg name The name in this module used in output messages.
     */
    Logger(const std::string name) : module(name) {}
    /*!
     * \brief destructor
     */
    ~Logger() = default;

    /*
     *
     * \brief Log implementation of this function
     *
     * Actual implementation of the log function.
     * At compile time evaluates to this implementation,
     * or an empty body implementation, depending on the
     * loglevel parameter of the Logger itself.
     *
     * \arg log Loglevel, either FATAL, ERROR, WARN, INFO, VERBOSE, DEBUG
     * \arg format Message format, where % can be used as a placeholder for
     * arguments. \arg arg... Any arguments which needs to be replaced.
     */
    template <Log LOGLEVEL, typename... Args>
    typename std::enable_if<!((L < LOGLEVEL) && (HPGEM_LOGLEVEL < LOGLEVEL)),
                            void>::type
        operator()(const LL<LOGLEVEL> log, const char* format, Args&&... arg) {
        std::stringstream msgstream;
        createMessage(msgstream, format, arg...);
        if (LOGLEVEL <= Log::FATAL) {
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

    template <Log LOGLEVEL, typename... Args>
    typename std::enable_if<(L < LOGLEVEL && HPGEM_LOGLEVEL < LOGLEVEL),
                            void>::type
        operator()(const LL<LOGLEVEL> log, const char* format, Args&&... arg) {}

    // std::string is sometimes convenient, but always slow, so where possible,
    // don't convert the const char* to a string before converting it back
    template <Log LOGLEVEL, typename... Args>
    void operator()(const LL<LOGLEVEL> log, std::string& format,
                    Args&&... arg) {
        (*this)(log, format.c_str(), std::forward<Args>(arg)...);
    }

    /*
     *
     * \brief Asserts on this logger
     *
     * Evaluates an assertion, prints an error message and aborts
     * in case of a failure. This message can be redirected,
     * and will be send to loggerOuptput->onFatal.
     *
     * \arg assertion An assertion, which must be true
     * \arg format Message format, where % can be used as a placeholder for
     * arguments. \arg arg... Any arguments which needs to be replaced.
     */

    // the conversion from "" to a std::sting is so slow, it takes 50% of the
    // total run time for a release build...
    template <typename... Args>
    typename std::enable_if<(ASSERTS) && (sizeof...(Args) >= 0), void>::type
        assert_debug(bool assertion, const char* format, Args&&... arg) {
        assert_always(assertion, format, std::forward<Args>(arg)...);
    }

    template <typename... Args>
    typename std::enable_if<!((ASSERTS) && sizeof...(Args) >= 0), void>::type
        assert_debug(bool assertion, const char* format, Args&&... arg) {}

    template <typename... Args>
    void assert_debug(bool assertion, const std::string format, Args&&... arg) {
        assert_debug(assertion, format.c_str(), std::forward<Args>(arg)...);
    }

    template <typename... Args>
    void assert_always(bool assertion, const char* format, Args&&... arg) {
        if (!assertion) {
            std::stringstream msgstream;
            createMessage(msgstream, format, std::forward<Args>(arg)...);
            loggerOutput->onFatal(module, msgstream.str());
        }
    }

    template <typename... Args>
    void assert_always(bool assertion, const std::string format,
                       Args&&... arg) {
        assert_always(assertion, format.c_str(), std::forward<Args>(arg)...);
    }

    /// Trigger handling for an unrecoverable failure.
    template <typename... Args>
    void fail [[noreturn]] (const char* format, Args&&... arg) {
        std::stringstream msgstream;
        createMessage(msgstream, format, std::forward<Args>(arg)...);
        loggerOutput->onFail(module, msgstream.str());
        // Fall back, should never be triggered for a correct onFail.
        exit(1);
    }

    template <typename... Args>
    void fail [[noreturn]] (const std::string format, Args&&... arg) {
        fail(format.c_str(), std::forward<Args>(arg)...);
    }

    /*!
     * \brief Oldskool log method.
     * \deprecated Use operator() instead.
     */
    template <typename... Args>
    void log(const Log loglevel, const std::string& format, Args&&... arg) {
        if (loglevel <= L || loglevel <= HPGEM_LOGLEVEL) {
            std::stringstream msgstream;
            createMessage(msgstream, format.c_str(),
                          std::forward<Args>(arg)...);
            if (loglevel <= Log::FATAL) {
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

    template <typename T>
    void suppressWarnings(T body) {
        auto oldWarn = loggerOutput->onWarn;
        loggerOutput->onWarn = [](std::string, std::string) {};
        body();
        loggerOutput->onWarn = oldWarn;
    }

   private:
    /*!
     * \brief Actual implementation to recursively replace all the '%' signs by
     * actual values.
     */
    template <typename Arg1, typename... Args>
    void createMessage(std::stringstream& msg, const char* fmt, Arg1&& arg,
                       Args&&... args) {
        bool doSkipNext = false;
        while (*fmt != '%' || doSkipNext) {
            // Make sure we're not running past the end of our formatting
            // string.
            if (*fmt == '\0') return;

            if (*fmt == '\\' && !doSkipNext) {  // Escape for the %sign
                doSkipNext = true;
                fmt++;
            } else {
                msg << *fmt;
                fmt++;
                doSkipNext = false;
            }
        }

        fmt++;  // Consume the % sign
        msg << arg;
        createMessage(
            msg, fmt,
            std::forward<Args>(args)...);  // and recursively call ourselve /
                                           // the method below.
    }

    /*!
     * \brief Terminating case / argument call
     */
    template <typename Arg1>
    void createMessage(std::stringstream& msg, const char* fmt, Arg1&& arg) {
        bool doSkipNext = false;
        while (*fmt != '%' || doSkipNext) {
            if (*fmt == '\0')  // End of string
                return;

            if (*fmt == '\\' && !doSkipNext) {  // Escape for the %sign and the
                                                // \sign
                doSkipNext = true;
                fmt++;
            } else {  // invoke the replacement
                msg << *fmt;
                fmt++;
                doSkipNext = false;
            }
        }
        fmt++;  // Consume the % sign
        msg << arg << fmt;
    }

    /*!
     * \brief Terminating case / no argument call
     */
    void createMessage(std::stringstream& msg, const char* message) {
        msg << message;
    }
};

/*! Default logger.
 * Use this for general logging.
 *
 * For very specific modules, define your own logger.
 * If you want to make extensive use of Debug messages,
 * please use a custom logger as well, to prevent polluting
 * the output.
 */
extern Logger<HPGEM_LOGLEVEL> logger;

// just emptying the functions is not sufficiently aggressive in disabling the
// actual (costly) comparison
#if !HPGEM_ASSERTS
#define assert_debug(e, ...) assert_debug(true, "")
#endif

}  // namespace hpgem

#endif  // HPGEM_KERNEL_LOGGER_H
