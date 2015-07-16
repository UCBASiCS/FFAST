#include "lcm.h"
#include <iostream>

int gcdPair(int a, int b)
{
    if (b == 0)
    {
        return a;
    }
    else
    {
        return gcdPair(b, a%b);
    }
}

int lcmPair(int a, int b)
{
    return a*b / gcdPair(a, b);
}

int lcmArray(int* array, int length)
{
    int lcm = 1;
    
    for (int i=0; i<length; ++i)
    {
        lcm = lcmPair(lcm,array[i]);
    }
    
    return lcm;
}

int lcmVector(std::vector<int> array) {
    int lcm = 1;
    for(auto it=array.cbegin(); it != array.cend(); ++it)
        lcm = lcmPair(lcm,*it);
    return lcm;
}

bool areCoprime(std::vector<int> array) {
    int d = array.size();
    bool result = true;
    int i=0;

    while (result && i < d) {
        int j=i+1;
        while (result && j < d) {
            result = (gcdPair(array[i],array[j]) == 1);
            ++j;
        }
        ++i;
    }
    return result;
}
