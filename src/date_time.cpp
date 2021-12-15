#include <date_time.hpp>
#include <cassert>
#include <casa/Quanta/MVEpoch.h>
#include <casa/Quanta/MVTime.h>
#include <sstream>

using namespace std;
int month2num (const std::string &month_name)
{
    const static std::string month_names[] = { "Jan", "Feb", "Mar", "Apr", "May", "Jun",
                                               "Jul", "Aug", "Sep", "Oct", "Nov", "Dec" };
    for (int i = 0; i < 12; ++i)
        {
            if (month_name == month_names[i])
                {
                    return i + 1;
                }
        }
    assert (0);
    return 0;
}

double parse_21cma_date (const std::string &date_string)
{
    using namespace casacore;
    std::string wd, mn, day, time, year;
    istringstream iss (date_string);
    iss >> wd >> mn >> day >> time >> year;
    std::string date_string2 = year + "-" + to_string (month2num (mn)) + "-" + day + "T" + time + "+8";
    Quantity qn;
    std::cerr << date_string2 << std::endl;
    MVTime::read (qn, date_string2);
    return qn.getValue ("s");
}
