/*
 * CommandLineOptions.hpp
 *
 *  Created on: Jun 16, 2014
 *      Author: dducks
 */

#ifndef COMMANDLINEOPTIONS_H_
#define COMMANDLINEOPTIONS_H_

#include <functional>
#include <map>
#include <vector>
#include <iostream>
#include <sstream>
#include <string>

namespace Base
{

    /**
     * Actually parses the arguments and uses the linked in CLO's
     * @param argc argc as specified by a main() call
     * @param argv argv as specified by a main() call
     * @return anything else than 0 means something went wrong.
     *         preferably crash the program and return the same value
     */
    int parse_options(int argc, char** argv);

    /**
     * Checks whether or not the arguments have been parsed.
     *
     */
    bool parse_isDone();

    /*
     * Don't go here. This is nasty and is hidden on purpose.
     * Black magic happens here, so don't.
     *
     * @dducks
     */
    namespace Detail
    {

        class CLOParser
        {
            std::size_t currCount;
            std::size_t count;
            char** pData;
            friend int Base::parse_options(int argc, char** argv);

            CLOParser(const CLOParser& other) = delete;
            CLOParser(CLOParser&& other) = delete;
            CLOParser& operator=(const CLOParser& other) = delete;
            CLOParser& operator=(CLOParser&& other) = delete;

            inline CLOParser(int argc, char** argv) : currCount(), count(argc), pData(argv) { }

            int go();
        public:

            inline std::string operator*() const
            {
                return std::string(pData[currCount]);
            }

            inline CLOParser& operator++()
            {
                currCount ++;
                return *this;
            }

            inline std::size_t remaining() const
            {
                return count - currCount;
            }
        };

        class CommandLineOptionBase
        {
        protected:
            char tag;
            std::string long_tag;
            std::string description;
            bool used;
            bool required;

            CommandLineOptionBase(char tag, std::string long_tag, std::string description, bool required)
            : tag(tag), long_tag(long_tag), description(description), used(false), required(required) { }
        public:

            inline char getTag()
            {
                return tag;
            }

            inline std::string getLongTag()
            {
                return long_tag;
            }
            virtual void parse(CLOParser& p) = 0;
            virtual bool hasArgument() = 0;

            bool isUsed()
            {
                return used;
            }

            bool isRequired()
            {
                return required;
            }

            std::string getDescription()
            {
                return description;
            }

            virtual ~ CommandLineOptionBase() { }

        };
        std::map<char, CommandLineOptionBase*>& getCLOMapping_short();
        std::map<std::string, CommandLineOptionBase*>& getCLOMapping_long();
        std::vector<CommandLineOptionBase*>& getCLOList();

        template<typename T>
        typename std::enable_if<std::is_same<T, bool>::value, T>::type
        parse_argument(CLOParser& p)
        {
            return true;
        }

        template<typename T>
        typename std::enable_if<std::is_arithmetic<T>::value &&
        ! std::is_same<T, bool>::value, T>::type
        parse_argument(CLOParser& p)
        {
            ++ p;
            std::istringstream arg(*p);
            T retVal = 0;
            arg >> retVal;
            if (! arg.eof())
            {
                std::string error = "Argument could not be converted to ";
                if (std::is_unsigned<T>::value)
                    error += "unsigned ";

                if (std::is_integral<T>::value)
                {
                    error += "integer.";
                }
                else if (std::is_floating_point<T>::value)
                {
                    error += "float.";
                }
                else
                {
                    error += typeid (T).name();
                }
                throw error;

            }

            return retVal;
        }

        /*
        template<typename T>
        typename std::enable_if<std::is_same<T, ::Vec3D>::value, T>::type
        parse_argument(CLOParser& p) {
            Vec3D retVal;
            std::size_t idx;
            if (p.remaining() < 4) {
                std::string error = "Not enough parameters.";
                throw error;
            }
            for (idx = 0; idx < 3; idx++) {
                ++p;
        
                std::istringstream coord(*p);
                double tmp;
                coord >> tmp;
                if (!coord.eof()) {
                    std::string error = "Argument could not be converted to float.";
                    throw error;
                }
                retVal.setComponent(idx, tmp);
        
            }
            return retVal;
        } */

        template<typename T>
        typename std::enable_if<std::is_same<T, std::string>::value, T>::type
        parse_argument(CLOParser& p)
        {
            ++ p;
            std::string copy = * p;
            return copy;
        }
    } // namespace Detail


    template<typename T>
    class CommandLineOption;

    /**
     * Returns a REFERENCE, yes really REFERENCE which you need to bind by
     * guess what, REFERENCE, to the actual CLO. Bind this to some kind of
     * static variable, and use it within the code to set and get properties.
     * make sure it is initialized before actually calling the parse_arguments()
     * function as documented above the magic part.
     *
     * @argument tag is the short version of the argument, or '\0'/0 if no such tag
     * is available. It has to be prefixed by '-'.
     * @argument long_tag is the long version of the argument, which when entered needs to
     * be prefixed with '--'.
     * @argument Description description of the argument, will show up in the --help command
     * @argument required Crash if it isn't supplied.
     */
    template<typename T>
    CommandLineOption<T>& register_argument(char tag, std::string long_tag,
            std::string description, bool required = false, T defaultValue = T());

    template<typename T>
    class CommandLineOption : public Detail::CommandLineOptionBase
    {
        /*template<typename T2>
        friend CommandLineOption<T2>& register_argument(char tag, std::string long_tag,
                std::string description,bool required,T2 defaultValue);
         */


        T value;

        CommandLineOption(CommandLineOption&& other) = delete; //yeah you're doing something wrong.
        CommandLineOption(const CommandLineOption& other) = delete; //lets prevent that mistake...
        CommandLineOption& operator=(CommandLineOption&& other) = delete;
        CommandLineOption& operator=(const CommandLineOption& other) = delete;
    public:

        CommandLineOption(char tag, std::string long_tag, std::string description, bool required, T defaultValue)
        : CommandLineOptionBase(tag, long_tag, description, required), value(defaultValue) { }

        void parse(Detail::CLOParser& p) override
        {
            value = Detail::parse_argument<T>(p);
            used = true;
        }

        T& getValue()
        {
            return value;
        }

        bool hasArgument() override
        {
            return std::is_same<T, bool>::value ? false : true;
        }


    };

    template<typename T>
    CommandLineOption<T>& register_argument(char tag, std::string long_tag,
            std::string description, bool required, T defaultValue)
    {
        CommandLineOption<T> * t = new CommandLineOption<T>{tag, long_tag, description, required, defaultValue};
        std::map<std::string, Detail::CommandLineOptionBase*> & lmapping = Detail::getCLOMapping_long();
        std::map<char, Detail::CommandLineOptionBase*> & smapping = Detail::getCLOMapping_short();
        std::vector<Detail::CommandLineOptionBase*> & req_list = Detail::getCLOList();
        if (lmapping.find(long_tag) == lmapping.end())
        {
            lmapping[long_tag] = t;
        }
        if (tag != '\0' && smapping.find(tag) == smapping.end())
        {
            smapping[tag] = t;
        }
        //   if (required) {
        req_list.push_back(t);
        //  }

        return *t;
    }

}

#endif /* COMMANDLINEOPTIONS_H_ */
