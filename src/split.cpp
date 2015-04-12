#include "split.h"

#include <sstream>

std::vector<std::string>& split(const std::string& string, char delimeter, std::vector<std::string>& elements)
{
    std::stringstream stringStream(string);
    std::string item;

    while (std::getline(stringStream, item, delimeter))
    {
        elements.push_back(item);
    }

    return elements;
}

std::vector<std::string> split(const std::string& string, char delimeter)
{
    std::vector<std::string> elements;

    split(string, delimeter, elements);

    return elements;
}

std::vector<std::string> split(const char*& string, char delimeter)
{
    std::string tempString(string);

    return split(tempString, delimeter);
}
