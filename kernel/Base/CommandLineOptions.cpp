/*
 * CommandLineOptions.cpp
 *
 *  Created on: Jun 16, 2014
 *      Author: dducks
 */

#include "CommandLineOptions.hpp"
#include <cstring>

std::map<std::string, Base::Detail::CommandLineOptionBase*>&
Base::Detail::getCLOMapping_long() {
    static std::map<std::string, CommandLineOptionBase*> mapping;
    return mapping;
}

std::map<char, Base::Detail::CommandLineOptionBase*>&
Base::Detail::getCLOMapping_short() {
    static std::map<char, CommandLineOptionBase*> mapping;
    return mapping;
}

std::vector<Base::Detail::CommandLineOptionBase*> &
Base::Detail::getCLOList() {
    static std::vector<Base::Detail::CommandLineOptionBase*> list;
    return list;
}

static auto& printHelp = Base::register_argument<bool>('?', "help", "Prints this help message", false, false);

static bool hasParsed = false;

bool Base::parse_isDone() {
    return hasParsed;
}

int Base::parse_options(int argc, char** argv) {
    if (hasParsed)
        throw ("Arguments have already been parsed");
    Base::Detail::CLOParser parser(argc, argv);
    hasParsed = true;
    return parser.go();
}

int Base::Detail::CLOParser::go()
{
    auto& lmapping = getCLOMapping_long();
    auto& smapping = getCLOMapping_short();

    ++(*this);
    try
    {
        while (remaining())
        {

            std::string tag = *(*this);

            Base::Detail::CommandLineOptionBase* pBase;


            if (tag.size() >= 2 && tag[0] == '-') {
                if (tag[1] == '-') {
                    std::string realtag = tag.substr(2);
                    auto it = lmapping.find(realtag);
                    if (it == lmapping.end()) {
                        std::cerr << "Argument " << tag << " not found.\n";
                        std::exit(1);
                    }

                    if (it->second->hasArgument() && remaining() <= 1) {
                        std::cerr << "Argument " << tag << "has an argument, but was not provided one.\n";
                        std::exit(1);
                    }
                        
                    pBase = it->second;

                    pBase->parse(*this);
                } else {
                    for (std::string::size_type i = 1; i < tag.size(); i++) {
                        auto it = smapping.find(tag[i]);
                        if (it == smapping.end()) {
                            std::cerr << "Argument " << tag << " not found.\n";
                            std::exit(1);
                        }

                        pBase = it->second;
                        if (pBase->hasArgument() && (i < (tag.size() - 1) || remaining() <= 1)) {
                            std::cerr << "Argument " << tag << "has an argument, but was not provided one.\n";
                            std::exit(1); // Tag with arg, with more flags following..
                        }
                        pBase->parse(*this);
                    }
                }
            } else {
                std::cerr << "Unknown argument '" << tag << "' found.\n";
                std::exit(1);
            }

            ++(*this);
            //if ()
        }
    } catch (const std::string& explain) {
        std::cerr << explain << std::endl;
        for (int i = 0; i < count; i++) {
            std::cerr << pData[i] << ' ';
        }
        std::cerr << std::endl;
        for (int i = 0; i < currCount; i++) {
            for (int j = 0; j < std::strlen(pData[i]); j++) {
                std::cerr << ' ';
            }
            std::cerr << ' ';
        }
        for (int j = 0; j < std::strlen(pData[currCount]); j++) {
            std::cerr << '^';
        }
        std::cerr << std::endl;
        std::exit(1);
    }

    bool kill = false;
    for (Base::Detail::CommandLineOptionBase* base : Base::Detail::getCLOList()) {
        if (base->isRequired() && !base->isUsed()) {
            std::cerr << "Argument \'" << base->getLongTag() << "\' required.\n";
            kill = true;
        }
    }
    if (printHelp.getValue()) {
        for (Base::Detail::CommandLineOptionBase* base : Base::Detail::getCLOList()) {
            std::cout << '\t';
            if (base->getTag() != '\0')
                std::cout << "'-" << base->getTag() << "', ";
            std::cout << "'--" << base->getLongTag() << "' ";
            if (base->hasArgument())
                std::cout << "[args...]";
            if (base->isRequired())
                std::cout << " - REQUIRED";
            std::cout << "\n\t  " << base->getDescription() << "\n\n";
        }
        kill = true;
    }
    
    if (kill)
        std::exit(1);
    return 0;
}
