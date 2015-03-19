/*
 This file forms part of hpGEM. This package has been developed over a number of years by various people at the University of Twente and a full list of contributors can be found at
 http://hpgem.org/about-the-code/team
 
 This code is distributed using BSD 3-Clause License. A copy of which can found below.
 
 
 Copyright (c) 2014, University of Twente
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 
 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef COMMANDLINEOPTIONS_H_
#define COMMANDLINEOPTIONS_H_

#include <functional>
#include <map>
#include <vector>
#include <iostream>
#include <sstream>
#include <string>
#include "Logger.h"

namespace Base
{
    
    /**
     * Actually parses the arguments and uses the linked in CLO's
     * @param argc argc as specified by a main() call
     * @param argv argv as specified by a main() call
     */
    void parse_options(int argc, char** argv);
    
    /**
     * Checks whether or not the arguments have been parsed.
     */
    bool parse_isDone();
    
    /*
     * The code below this comment is difficult to understand, so don't change
     * anything unless you know what you're doing and you test your changes 
     * thoroughly before committing them.
     */
    namespace Detail
    {
        
        class CLOParser
        {
            std::size_t currCount;
            std::size_t count;
            char** pData;
            friend void Base::parse_options(int argc, char** argv);

            CLOParser(const CLOParser& other) = delete;
            CLOParser(CLOParser&& other) = delete;
            CLOParser& operator=(const CLOParser& other) = delete;
            CLOParser& operator=(CLOParser&& other) = delete;
            
            inline CLOParser(int argc, char** argv)
                    : currCount(), count(argc), pData(argv)
            {
            }
            
            //does the actual parsing
            //returns the number of arguments parsed
            int go();
        public:
            
            inline char* operator*() const
            {
                return pData[currCount];
            }
            
            inline CLOParser& operator++()
            {
                currCount++;
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
                    : tag(tag), long_tag(long_tag), description(description), used(false), required(required)
            {
            }
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
            
            virtual ~ CommandLineOptionBase()
            {
            }
            
        };
        std::map<char, CommandLineOptionBase*>& getCLOMapping_short();
        std::map<std::string, CommandLineOptionBase*>& getCLOMapping_long();
        std::vector<CommandLineOptionBase*>& getCLOList();
        
        template<typename T>
        typename std::enable_if<std::is_same<T, bool>::value, T>::type parse_argument(CLOParser& p)
        {
            return true;
        }
        
        template<typename T>
        typename std::enable_if<std::is_arithmetic<T>::value && !std::is_same<T, bool>::value, T>::type parse_argument(CLOParser& p)
        {
            ++p;
            std::istringstream arg(*p);
            T retVal = 0;
            arg >> retVal;
            if (!arg.eof())
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
                    error += typeid(T).name();
                }
                logger(ERROR, error);
                
            }
            
            return retVal;
        }
        

        template<typename T>
        typename std::enable_if<std::is_same<T, std::string>::value, T>::type parse_argument(CLOParser& p)
        {
            ++p;
            std::string copy = *p;
            return copy;
        }
    } // namespace Detail
    
    template<typename T>
    class CommandLineOption;
    
    /**
     * Returns a REFERENCE which you need to bind by REFERENCE to the actual CLO.
     * Bind this to some kind of static variable, and use it within the code to
     * set and get properties. Make sure it is initialized before actually calling
     * the function parse_options(int, char**).
     *
     * @argument tag is the short version of the argument, or '\0'/0 if no such tag
     * is available. It has to be prefixed by '-'.
     * @argument long_tag is the long version of the argument, which when entered needs to
     * be prefixed with '--'.
     * @argument Description description of the argument, will show up in the --help command
     * @argument required Crash if it isn't supplied.
     */
    template<typename T>
    CommandLineOption<T>& register_argument(char tag, std::string long_tag, std::string description, bool required = false, T defaultValue = T());
    
    template<typename T>
    class CommandLineOption : public Detail::CommandLineOptionBase
    {

        T value;

        CommandLineOption(CommandLineOption&& other) = delete; //yeah you're doing something wrong.
        CommandLineOption(const CommandLineOption& other) = delete; //lets prevent that mistake...
        CommandLineOption& operator=(CommandLineOption&& other) = delete;
        CommandLineOption& operator=(const CommandLineOption& other) = delete;
    public:
        
        CommandLineOption(char tag, std::string long_tag, std::string description, bool required, T defaultValue)
                : CommandLineOptionBase(tag, long_tag, description, required), value(defaultValue)
        {
        }
        
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
    CommandLineOption<T>& register_argument(char tag, std::string long_tag, std::string description, bool required, T defaultValue)
    {
        CommandLineOption<T> * t = new CommandLineOption<T> {tag, long_tag, description, required, defaultValue};
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
        req_list.push_back(t);
        
        return *t;
    }

}

#endif /* COMMANDLINEOPTIONS_H_ */
