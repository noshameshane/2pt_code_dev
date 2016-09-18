#include "statistics.h"

int STAT::init()
{
    max=-1e99; min=-max;
    sum=sum2=avg=std=0;
    npt=0;
    if(binflag) init_bin();
    if(dataflag) init_XYdata();
    if(lsfitflag) init_XYdata();
    return 0;
}

int STAT::add_data(double newdata)
{
    if(newdata>max) max=newdata;
    if(newdata<min) min=newdata;
    sum+=newdata;
    sum2+=newdata*newdata;
    npt++;
    if(binflag) add_bin(newdata);
    if(dataflag) add_Ydata(newdata);
    return 0;
}

int STAT::add_data(double newx,double newy)
{
    if(newy>max) max=newy;
    if(newy<min) min=newy;
    sum+=newy;
    sum2+=newy*newy;
    npt++;
    if(binflag) add_bin(newy);
    if(dataflag) add_Ydata(newy);
    if(lsfitflag) add_XYdata(newx,newy);
    return 0;
}

int STAT::normalize(double factor)
{
    sum/=factor;
    sum2/=factor*factor;
    return 0;
}

int STAT::cal_avg()
{
    if(npt==0) return 1;
    else if(npt==1) {
       avg=sum;
       std=-999.99;
    } else {
       avg=sum/npt;
       if(max==min) std=0;
       else std=sqrt((sum2-avg*avg*npt)/(npt-1));
       //else std=sqrt((sum2/(npt-1)-sum*sum/npt/(npt-1)));
    }
    return 0;
}

int STAT::add_Ydata(double y)
{
    if(npt>nXY) return 1;
    XYdata[npt-1][1]=y;
    return 0;
}

int STAT::init_bin()
{
    int i;
    nbinpt=nbinexceed=0;
    nbin=(int)((bmax-bmin)/bsize)+1;   
    if(nbin<=0) { 
        binflag=0;
        return 1;
    }
    if(bin!=NULL) delete [] bin;
    bin=new double [nbin];
    for(i=0;i<nbin;i++) bin[i]=0;
    return 0;
}

int STAT::add_bin(double newdata)
{
    int index;
    index=(int) ((newdata-bmin)/bsize);
    if(index>=nbin || index<0) { 
      nbinexceed++;
    } else {
      bin[index]++;
      nbinpt++;
    }
    return 0;
}

int STAT::normalize_bin()
{
    if(binflag==0) return 1;
    int i;

    for(i=0;i<nbin;i++) bin[i]/=nbinpt;

    return 0;
}

int STAT::init_XYdata()
{
    if(nXY==0) return 1;
    if(XYdata!=NULL) delete [] XYdata;
    XYdata=new double [nXY][2];
    Ax=Ay=Sxx=Sxy=Syy=0;
    sigma2x=sigma2y=covxy=R2=s2=s=0; //variance,correlation coefficient, variance
    a=b=0;   //Y= a + b X
    SEa=SEb=0;   //standard errors in aXY and bXY
    
    return 0;
}

int STAT::add_XYdata(double x,double y)
{
    if(npt>nXY) return 1;
    XYdata[npt-1][0]=x;
    XYdata[npt-1][1]=y;
    Ax+=x;
    Ay+=y;
    Sxx+=x*x;
    Syy+=y*y;
    Sxy+=x*y;
    return 0;
}

int STAT::lsfit_XYdata()
{
    Ax /=npt;
    Ay /=npt;
    Sxx =(Sxx-npt*Ax*Ax);
    Syy =(Syy-npt*Ay*Ay);
    Sxy =(Sxy-npt*Ax*Ay);
    sigma2x=Sxx/npt;
    sigma2y=Syy/npt;
    covxy=Sxy/npt;
    b  = Sxy/Sxx;  //=covxy/sigma2x
    a  = Ay-(b*Ax);
    R2 = (Sxy*Sxy)/(Sxx*Syy);
    s  = sqrt( ( Syy-b*Sxy )/(npt-2) );
    SEa = s *sqrt(1.0/npt+Ax*Ax/Sxx);
    SEb = s /sqrt(Sxx);
    //printf("npt %d Ax %f Ay %f Sxx %f Syy %f Sxy %f b %f a %f R2 %f SEa %f SEb %f\n",npt,Ax,Ay,Sxx,Syy,Sxy,b,a,R2,SEa,SEb);
    return 0;
}

