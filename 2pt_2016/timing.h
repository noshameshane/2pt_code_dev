#ifndef TIMING_H
#define TIMING_H
#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>

using namespace std;

class TIME
{
public:
  TIME () {
   timezero=timelast=timecurrent=0;
   clockzero=clocklast=clockcurrent=0;
   day=hr=min=sec=0;
   cday=86400;
   chr=3600;
   cmin=60;
   csec=CLOCKS_PER_SEC;
   settimezero();
  }
  ~TIME() {
  }
  time_t timezero;
  time_t timelast;
  time_t timecurrent;
  clock_t clockzero;
  clock_t clocklast;
  clock_t clockcurrent;
  time_t day;
  time_t hr;
  time_t min;
  time_t sec;
  double  dsec;
  time_t cday;  //# seconds per day
  time_t chr;   //# seconds per hour
  time_t cmin;  //# seconds per min
  time_t csec;  //# seconds per sec

  int myclock();
  int settimezero();
  int gettimetotal(ostream *);
  int gettimeincrement(ostream *);
  int getclockincrement(ostream *);
  int outtime(time_t,time_t,clock_t,clock_t,char *,ostream *);
  int outclock(clock_t,clock_t,char *,ostream *);
  
};

#endif
