#ifndef PRIVATE_LOAD_SIGNAL_HPP
#define PRIVATE_LOAD_SIGNAL_HPP
#include <filesystem>
#include <string>
#include <fstream>
#include <uSignal/vector.hpp>

namespace
{

template<typename T>
USignal::Vector<T> loadSignal(const std::filesystem::path &fileName)
{
    if (!std::filesystem::exists(fileName))
    {   
        throw std::invalid_argument(std::string {fileName} + " does not exist");
    }
    std::ifstream signalFile;
    signalFile.open(fileName);
    USignal::Vector<T> signal;
    signal.reserve(12000);
    std::string line;
    while (std::getline(signalFile, line))
    {   
        double yi; 
        std::sscanf(line.c_str(), "%lf\n", &yi);
        signal.push_back(static_cast<T> (yi));
    }   
    signalFile.close();
    return signal;
}

}
#endif
