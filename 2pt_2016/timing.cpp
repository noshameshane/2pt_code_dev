#include "timing.h"

int TIME::myclock()
{
    clockcurrent=clock();
    timecurrent=time(NULL);
    return 0;
}

int TIME::settimezero()
{
    myclock();
    timelast=timezero=timecurrent;
    clocklast=clockzero=clockcurrent;
    return 0;
}

int TIME::gettimetotal(ostream *outf)
{
    myclock();
    return outtime(timezero,timecurrent,clockzero,clockcurrent,"Total",outf);
}

int TIME::gettimeincrement(ostream *outf)
{
    timelast=timecurrent;
    clocklast=clockcurrent;
    myclock();
    return outtime(timelast,timecurrent,clocklast,clockcurrent,"Increment",outf);
}

int TIME::getclockincrement(ostream *outf)
{
    timelast=timecurrent;
    timecurrent=myclock();
    return outclock(timelast,timecurrent,"Increment",outf);
}

int TIME::outtime(time_t itime,time_t ftime,clock_t iclock,clock_t fclock,char *str,ostream *outf)
{
    time_t difft;
    clock_t diffc;

    difft=ftime-itime;
    diffc=fclock-iclock;

    day=  difft/cday; 
    hr = (difft%cday)/chr;
    min=((difft%cday)%chr)/cmin;
    sec=(((difft%cday)%chr)%cmin);
    dsec=sec+(diffc%csec)/(double)csec;
    char null[1024];
    sprintf(null,"%s time used %d day %d hrs %d min %.2f sec.\n",str,day,hr,min,dsec);
    *outf<<null<<endl;

    return 0;
}

int TIME::outclock(clock_t iclock,clock_t fclock,char *str,ostream *outf)
{
    clock_t diff;

    diff=(fclock-iclock);

    char null[1024];
    sprintf(null,"%s clock ticks is %d, %f sec.\n",str,diff,double(diff)/csec);
    *outf<<null<<endl;

    return 0;
}

