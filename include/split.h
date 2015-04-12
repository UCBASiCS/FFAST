#ifndef HELPERS_H
#define HELPERS_H

#include <string>
#include <vector>

std::vector<std::string>& split(const std::string& string, char delimeter, std::vector<std::string>& elements);
std::vector<std::string> split(const std::string& string, char delimeter);
std::vector<std::string> split(const char*& string, char delimeter);

#endif // HELPERS_H
