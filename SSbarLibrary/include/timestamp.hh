#ifndef timestamp_h
#define timestamp_h

//2015-02-01 - taken from 2013 Lorentz Angle work. Functions that spit out time stamps formatted for file names.

#include <sstream>
#include <string>
#include <vector>
#include <TDatime.h>

#include <ctime>
#include <chrono>
#include <iomanip>

// Spits out date and time, by default in YYYY-MM-DD HH:MM:SS
static std::string timeStamp(char d = '-', char b = ' ', char t = ':');
// Spits out date and time, by default in YYYY-MM-DD_HHMMSS, for use in filenames
static std::string fileTimeStamp();
// Calculates time duration
static std::string timeDuration( std::string Tstart, std::string Tend );


static std::string timeStamp( char d, char b, char t ) {
  // print_date_time - takes characters to separate date & time.
  // if b and/or t are equal to '0', only prints date.
    std::stringstream output;
    TDatime now;
    int date = now.GetDate();
    int time = now.GetTime();
    if( time < 0 )
    {
        time-= 240000;
        date+= 1;
    }

    output << date/10000 << d << (date%10000)/1000 << (date%1000)/100 << d << (date%100)/10 << date%10;
    if( b!='0' && t!='0' ) output << b << time/10000 << t << (time%10000)/1000 << (time%1000)/100 << t << (time%100)/10 << time%10;
return output.str();
}


// Spits out date and time, by default in YYYY-MM-DD_HHMMSS, for use in filenames
static std::string fileTimeStamp() {
    std::stringstream output;
    TDatime now;
    int date = now.GetDate();
    int time = now.GetTime();
    if( time < 0 )
    {
        time-= 240000;
        date+= 1;
    }

    output << date/100000 << (date%100000)/10000 << '-' << (date%10000)/1000 << (date%1000)/100 << '-' << (date%100)/10 << date%10;
    output << '_' << time/100000 << (time%100000)/10000 << (time%10000)/1000 << (time%1000)/100 << (time%100)/10 << time%10;
return output.str();
}


std::string timeDuration(std::string Tstart, std::string Tend) {
    std::tm tm_start = {}, tm_end = {};
    std::istringstream ss_start(Tstart), ss_end(Tend);
    ss_start >> std::get_time(&tm_start, "%Y-%m-%d %H:%M:%S");
    ss_end >> std::get_time(&tm_end, "%Y-%m-%d %H:%M:%S");
    auto time_start = std::chrono::system_clock::from_time_t(std::mktime(&tm_start));
    auto time_end = std::chrono::system_clock::from_time_t(std::mktime(&tm_end));
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(time_end - time_start);
    int hours = duration.count() / 3600;
    int minutes = (duration.count() % 3600) / 60;
    int seconds = duration.count() % 60;
    std::ostringstream output;
    output << std::setw(2) << std::setfill('0') << hours << ":"
           << std::setw(2) << std::setfill('0') << minutes << ":"
           << std::setw(2) << std::setfill('0') << seconds;
    return output.str();
}

#endif