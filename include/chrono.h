#ifndef CHRONO_H
#define CHRONO_H

#include <map>
#include <string>
#include <sys/time.h>
#include <time.h>
#include <vector>

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif // __MACH__

class Chrono
{
private:
    timespec time;
    std::map< std::string, std::vector<double> > sectionsStartTimes;
    std::map< std::string, std::vector<double> > sectionsStopTimes;
    std::map< std::string, std::vector<double> > sectionsIterationTimes;

public:
    Chrono() {}
    void start(const std::string& section);
    void stop(const std::string& section);
    double iteration(const std::string& section);
    double average(const std::string& section);

private:
    double getTime();
    double getSectionTime(const std::string& section, int index);
    bool notEmpty(const std::string& section) const;
    std::vector<double>& getSectionStartTimes(const std::string& section);
    std::vector<double>& getSectionStopTimes(const std::string& section);
};

#endif // CHRONO_H
