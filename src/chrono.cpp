#include "chrono.h"

#include <iostream>

void Chrono::start(const std::string& section)
{
    getSectionStartTimes(section).push_back(getTime());
}

void Chrono::stop(const std::string& section)
{
    getSectionStopTimes(section).push_back(getTime());
}

double Chrono::iteration(const std::string& section)
{
    if (notEmpty(section))
    {
        return getSectionTime(section, getSectionStopTimes(section).size()-1);
    }

    return 0;
}

double Chrono::average(const std::string& section)
{
    if (notEmpty(section))
    {
        double sum = 0;
        unsigned int timesNb = getSectionStopTimes(section).size();

        for (unsigned int i=0; i<timesNb; i++)
        {

	    sum += getSectionTime(section, i);
        }

        return sum/timesNb;
    }

    return 0;
}

double Chrono::getTime()
{
#ifdef __MACH__
    clock_serv_t cclock;
    mach_timespec_t mts;
    host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
    clock_get_time(cclock, &mts);
    mach_port_deallocate(mach_task_self(), cclock);
    time.tv_sec = mts.tv_sec;
    time.tv_nsec = mts.tv_nsec;
#else
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time);
#endif

    return time.tv_sec + time.tv_nsec * 1e-9;
}

double Chrono::getSectionTime(const std::string& section, int index)
{
    return getSectionStopTimes(section).at(index) - getSectionStartTimes(section).at(index);
}

bool Chrono::notEmpty(const std::string& section) const
{
    if (sectionsStopTimes.count(section) > 0)
    {
        return true;
    }

    return false;
}

std::vector<double>& Chrono::getSectionStartTimes(const std::string& section)
{
    return sectionsStartTimes[section];
}

std::vector<double>& Chrono::getSectionStopTimes(const std::string& section)
{
    return sectionsStopTimes[section];
}
