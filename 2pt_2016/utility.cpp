/* Routines from Numerical Recipes (slightly modified). */
#include "utility.h"

int ROTATE(double **a,int i,int j,int k,int l,double *tau,double *s)
{
    double g,h;
    g=a[i][j];
    h=a[k][l];
    a[i][j]=g-(*s)*(h+g*(*tau));
    a[k][l]=h+(*s)*(g-h*(*tau));
    return 0;
}

int jacobi(double **a, int n, double *d, double **v)
{
     int j,iq,ip,i;
     double tresh,theta,tau,t,sm,s,h,g,c,*b,*z;

     b=new double [n];
     z=new double [n];
     for (ip=0;ip<n;ip++) {
          v[ip][ip]=1.0;
          for (iq=ip+1;iq<n;iq++) v[ip][iq]=v[iq][ip]=0.0;
     }
     for (ip=0;ip<n;ip++) {
          b[ip]=d[ip]=a[ip][ip];
          z[ip]=0.0;
     }
     int nrot=0;
     for (i=0;i<50;i++) {
          sm=0.0;
          for (ip=0;ip<n-1;ip++) {
               for (iq=ip+1;iq<n;iq++)
                    sm += fabs(a[ip][iq]);
          }
          if (sm == 0.0) {
               delete [] z;
               delete [] b;
               return 0;
          }
          if (i < 3) tresh=0.2*sm/(n*n);
          else tresh=0.0;
          for (ip=0;ip<n-1;ip++) {
               for (iq=ip+1;iq<n;iq++) {
                    g=100.0*fabs(a[ip][iq]);
                    if (i > 3 && (fabs(d[ip])+g) == fabs(d[ip])
                         && (fabs(d[iq])+g) == fabs(d[iq]))
                         a[ip][iq]=0.0;
                    else if (fabs(a[ip][iq]) > tresh) {
                         h=d[iq]-d[ip];
                         if ((fabs(h)+g) == fabs(h)) t=(a[ip][iq])/h;
                         else {
                              theta=0.5*h/(a[ip][iq]);
                              t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
                              if (theta < 0.0) t = -t;
                         }
                         c=1.0/sqrt(1+t*t);
                         s=t*c;
                         tau=s/(1.0+c);
                         h=t*a[ip][iq];
                         z[ip] -= h;
                         z[iq] += h;
                         d[ip] -= h;
                         d[iq] += h;
                         a[ip][iq]=0.0;
                         for (j=0;j<ip;j++) {
                              ROTATE(a,j,ip,j,iq,&tau,&s);
                         }
                         for (j=ip+1;j<iq;j++) {
                              ROTATE(a,ip,j,j,iq,&tau,&s);
                         }
                         for (j=iq+1;j<n;j++) {
                              ROTATE(a,ip,j,iq,j,&tau,&s);
                         }
                         for (j=0;j<n;j++) {
                              ROTATE(v,j,ip,j,iq,&tau,&s);
                         }
                         ++(nrot);
                    }
               }
          }
          for (ip=0;ip<n;ip++) {
               b[ip] += z[ip];
               d[ip]=b[ip];
               z[ip]=0.0;
          }
     }
     printf("Too many iterations in routine jacobi nrot %d sum %f tresh %f\n",nrot,sm,tresh);
     if(1) {
       printf("original matrix\n");
       for(i=0;i<n;i++) {
          for(j=0;j<n;j++) printf("%8.4e ",a[i][j]);
          printf("\n");
       }
     }
     delete [] z;
     delete [] b;
     return 1;
}

int eigsrt(double *d, double **v, int n)
{
    int k,j,i;
    double p;

    for (i=0;i<n-1;i++) {
      for (j=i+1;j<n;j++) {
         if (d[j] > d[i]) {
            SWAP(&d[j],&d[i]);
            for(k=0;k<n;k++) SWAP(&v[k][i],&v[k][j]);
         }
      }
    }
    return 0;
}

int eigsrt0(double *d, double **v, int n)
{
    int k,j,i;
    double p;

    for (i=0;i<n-1;i++) {
     p=d[k=i];
     for (j=i+1;j<n;j++)
          if (d[j] >= p) p=d[k=j];
     if (k != i) {
          d[k]=d[i];
          d[i]=p;
          for (j=0;j<n;j++) {
               p=v[j][i];
               v[j][i]=v[j][k];
               v[j][k]=p;
          }
     }
    }
    return 0;
}

int solve_cubic_eq(double a,double b,double c,double *c3rts)
{
   double Q,R,D;
   Q=(3*b-a*a)/9.0;
   R=(9*a*b-27*c-2*a*a*a)/54;
   D=Q*Q*Q+R*R;
   if(D>0) { //one real root
     double T=R-sqrt(D);
     double U=R+sqrt(D);
     if(U<0) c3rts[0]=-pow(-U,1.0/3.0)-pow(-T,1.0/3.0)-a/3.0;
     else if(T>0) c3rts[0]=pow(R+sqrt(D),1.0/3.0)+pow(T,1.0/3.0)-a/3.0;
     else c3rts[0]=pow(R+sqrt(D),1.0/3.0)-pow(-T,1.0/3.0)-a/3.0;
     c3rts[1]=c3rts[2]=c3rts[0];
   } else if (D<0) { //three unequal real roots (Q<0)
     double theta;
     double PI=acos(-1.0);
     theta=acos(R/sqrt(-1.0*Q*Q*Q));
     c3rts[0]=2*sqrt(-Q)*cos(theta/3.0)-a/3.0;
     c3rts[1]=2*sqrt(-Q)*cos((theta+2*PI)/3.0)-a/3.0;
     c3rts[2]=2*sqrt(-Q)*cos((theta+4*PI)/3.0)-a/3.0;
     if(c3rts[0]<c3rts[1]) SWAP(&c3rts[0],&c3rts[1]);
     if(c3rts[0]<c3rts[2]) SWAP(&c3rts[0],&c3rts[2]);
     if(c3rts[1]<c3rts[2]) SWAP(&c3rts[1],&c3rts[2]);
   } else { //D==0, at least two equal real roots
     if(R>0) {
        c3rts[0]=2*pow(R,1.0/3.0)-a/3.0;
        c3rts[1]=-pow(R,1.0/3.0)-a/3.0;
     } else {
        c3rts[0]=-2*pow(-R,1.0/3.0)-a/3.0;
        c3rts[1]=+pow(-R,1.0/3.0)-a/3.0;
     }
     if(c3rts[0]<c3rts[1]) {
         c3rts[2]=c3rts[0];
         c3rts[0]=c3rts[1];
     } else {
         c3rts[2]=c3rts[1];
     }
   }
   return 0;
}

/* The following subroutines calculates the integration of weighting    */
/* function for entropy from 0+ to 1/2 fmin. The algorithm is based on  */
/* mathematica result:                                                  */
/* w[u_] := (u/ (Exp[u] - 1) - Log[1 - Exp[-u]] )                       */
/* Integrate[  w[u], {u, 0, x}] =                                       */
/* (Pi^2)/2 -x^2 -x Log[1-Exp[-x]] + 2 PolyLog[2,Exp[-x]])              */

double scweighting(double upper)
{
  if(upper==0) return 0;
  return (2*upper-upper*log(upper))/upper;
}

double sqweighting(double upper)
{
  if(upper==0) return 0;
  double pi;
  double wsq;
  pi=3.1415926535897932385;

  wsq= pi*pi/3.0 -upper*upper +upper*log(-1+exp(upper))
     - 2*polylog(2.0,exp(-upper));

  return wsq/upper;
}

double polylog(double n, double z)
{
 int k;
 double sum,sum_old;
 k=1;
 sum=0.0;
 sum_old=0.0;

 do {
      sum_old=sum;
      sum+= pow(z,k)/pow(k,n);
      k++;
/*    printf("\n n=%f z=%f sum=%f sum_old=%f k=%d",n,z,sum,sum_old,k); */
   } while (fabs(sum-sum_old)!=0);
 return sum;
}

double HSDF(double *pwr,int lenth,double fmin,int nmol,double fract_f,double fmf)
{
     double hsdf,tmpg,tmps;
     int j;
     hsdf=tmpg=tmps=0;
     for(j=0;j<lenth;j++) {
         twoPT(&tmpg,&tmps,pwr[0],pwr[j],fmin*j,nmol,fract_f,fmf);
         if(j==0 || j==lenth-1) { hsdf += tmpg*fmin*0.5;}
         else { hsdf += tmpg*fmin; }
     }
     return hsdf;
}

double TTDF(double *pwr,int lenth,double fmin)
{
     double ttdf;
     int j;
     ttdf=0;
     for(j=0;j<lenth;j++) {
         if(j==0 || j==lenth-1) { ttdf += pwr[j]*fmin*0.5;}
         else { ttdf += pwr[j]*fmin; }
     }
     return ttdf;
}

int HSweighting(double *wep,double *wap,double *wcvp,double *wsehs,double *wsp,double *wspd,double y,double mass,double nmol,double transT,double rotT,double volume,double *wsr,double *war,double *rT,double rs)
{
    *wep=*wap=*wcvp=*wsehs=*wsp=*wspd=*wsr=*war=0;
    //Ideal gas
    if(y>0.74) return 0; //packing fraction too large, ignore hard sphere correction
    else {
       *wsp= 5.0/2.0 + log(pow(2.0*PI*mass*1e-3/Na*kb*transT/h/h,1.5)*volume*1e-30/(nmol)); // indistinguishable particles 
       *wspd= 3.0/2.0 + log(pow(2.0*PI*mass*1e-3/Na*kb*transT/h/h,1.5)*volume*1e-30); // distinguishable particles 

       //Carnahan-Starling
       *wsehs=log((1.0+y+y*y-y*y*y)/pow(1.0-y,3.0))+y*(3.0*y-4.0)/pow(1.0-y,2.0);

       //Hard sphere gas
       *wsp  = (*wsp  + *wsehs)/3.0;
       *wspd = (*wspd + *wsehs)/3.0;
       *wep  = 0.5;
       *wap  = *wep-*wsp; /*wa=we-ws, no need of factor T in front of ws*/
       *wcvp = 0.5;

       //rigid rotor
       if(rT[0]>0) { //for none-atomic molecules
         if(rT[2]<0) { //linear
           *wsr  = (1.0 + log(rotT/sqrt(rT[0]*rT[1])/rs))/2.0;
         } else { //non-linear
           *wsr  = (3.0/2.0 + log(sqrt(PI*rotT*rotT*rotT/(rT[0]*rT[1]*rT[2]))/rs))/3.0;
         }
         
         *war  = *wep-*wsr;
       } else *wsr = *war = 0;
       cout<<"\nstlin T "<<rotT<<" rT "<<rT[0]<<" "<<rT[1]<<" "<<rT[2]<<" rs "<<rs<<" ws "<<*wsr<<endl;
    }
    return 0;
}

int twoPTmf(double *pwr,double s0,int n,double pwrfreq,int nmol,double f2pt,double &fmf)
{   //by Ming-Hsien 2015
    int i;
    double A,B ;
    fmf=0.9*f2pt; //initial guess for fmf  

    cal_AB(A,B,s0,f2pt,fmf);
    // printf("A %f B %f f2pt %f fmf %f s0 %f nmol %d\n",A,B,f2pt,fmf,s0,nmol);
    if(0) {
      for(i=0;i<n;i++) cout<<i*pwrfreq<<" "<<Sgmf(i*pwrfreq,s0,f2pt,nmol,fmf)<<endl;
    }

    double nxtfmf=0;
    double dev=fabs(nxtfmf/fmf-1);
    double sum,min,tmp;
    int nint=n-1,cc=0;
    tmp=s0;
    while (tmp>pwr[nint])
    {
        nint--;
        tmp=s0/(1.0+pow(PI*s0*pwrfreq*nint/(6.0*nmol*f2pt),2.0));
    }

    while(dev>1e-8) {
       sum=0;min=s0;
       for(i=nint;i<n;i++) { 
          tmp=Sgmf(i*pwrfreq,s0,f2pt,nmol,fmf)-pwr[i];
          if(tmp>0) sum+=tmp*pwrfreq;//sum=tmp/Sgmf(i*pwrfreq,s0,f2pt,nmol,fmf);//tmp;
          else if(min>fabs(tmp)) min=fabs(tmp); 
       }
       if(sum>0) min=0; //test
       nxtfmf=fmf*(1.0-sum/(3.0*nmol)+min/pwr[nint]);
       dev=fabs(nxtfmf/fmf-1);
       //printf("loop=%d fmf=%e nxtfmf=%e dev=%e",cc,fmf,nxtfmf,dev); cout<<endl;
       if (nxtfmf>1) fmf=f2pt;
       else fmf=nxtfmf;
       cc++;
    }
    cout<<endl<<"fmf loop="<<cc<<" fhs="<<f2pt<<" fmf="<<fmf<<endl;
    if(0) {
      for(i=0;i<n;i++) { 
        tmp=s0/(1.0+pow(PI*s0*pwrfreq*i/(6.0*nmol*f2pt),2.0));
        // cout<<i*pwrfreq<<" "<<Sgmf(i*pwrfreq,s0,f2pt,nmol,fmf)<<" "<<tmp<<" "<<pwr[i]<<endl;
      }
    }

    return 0;
}

double cal_AB(double &A, double &B, double s0, double f2pt, double fmf)
{

    double alpha,a,b,c;
    alpha=f2pt/s0;
    a=16.0*pow(PI*(s0*s0/(fmf*fmf)-1.0/(alpha*alpha)),2.0);
    b=-8.0*PI*((4.0+PI)*(s0*s0/(fmf*fmf))+(4.0-PI)/(alpha*alpha));
    c=(PI-4.0)*(PI-4.0);
    B=(-b+pow(pow(b,2.0)-4.0*a*c,0.5))/2/a;       //in unit of
    A=pow(4*B*pow(fmf,2.0)/PI/pow(s0,2.0),0.5);   //in unit of

    return 0;
}

double Sgmf(double v,double s0, double f2pt,int nmol,double fmf)
{
   //gas component DoS in unit of cm
   double A,B;
   cal_AB(A,B,s0,f2pt,fmf);
   double vv,y;
   vv=PI*v/6.0/nmol;
   y=vv/pow(4.0*B,0.5);
   double D,c1,c2;
   c1=214.0/155.0;
   c2=503.0/754.0;
   if(y==0) {D=0;}
   else{
     D=(1.0-exp(-c1*y*y))/2.0/y+(1.0-(1.0+y*y)*exp(-y*y))/4.0/pow(y,3.0)+(0.875-0.50*c1)*y*exp(-c2*y*y);
     //Approximate solution to Dawson Integral : Reference:
   }
   if(fmf==f2pt)
   {return s0/(1.0+pow(PI*s0*v/(6.0*nmol*f2pt),2.0));}
   else{
   return (A*pow(PI/4.0/B,0.5)*exp(-pow(y,2.0)))*fmf/(pow(A,2.0)*PI/4.0/B*exp(-2*pow(y,2.0))+pow(A*D,2.0)/B-4.0*A*y*D+pow(vv,2.0));}
}

int twoPT(double *tmpg,double *tmps,double s0,double sv,double v,int nmol,double fract_f,double fmf)
{
    if(fract_f==0) { *tmps=sv; *tmpg=0; }
    else if(v==0) { *tmpg=s0; *tmps=0; }
    else {
       if(fract_f==fmf) *tmpg= s0/(1.0+pow(PI*s0*v/(6.0*nmol*fract_f),2.0)); //2PT
       else *tmpg=Sgmf(v,s0,fract_f, nmol, fmf);//s0/(1.0+pow(PI*s0*v/(6.0*nmol*fract_f),2.0)); //2PT-MF
       if(*tmpg>sv) *tmpg=sv;
       *tmps=sv-*tmpg;
    }
    return 0;
}

/*determine fraction factor f from constant K*/
double search2PT(double K)
{
        double P,fold,fnew,dPdf,tol;
        int count;
        fold=0.0;
        fnew=0.7293*pow(K,0.5727); /*initial guess*/
        if(fnew>0.5) fnew=0.5;
        tol=1e-10;
        count=0;

        while( fabs(fnew-fold)>tol && count<999) {
                dPdf=0.0;
                fold=fnew;
                P=2.0*pow(K,-4.5)*pow(fnew,7.5)-6.0*pow(K,-3.0)*pow(fnew,5.0)-pow(K,-1.5)*pow(fnew,3.5)
                 +6.0*pow(K,-1.5)*pow(fnew,2.5)+2.0*fnew-2;
                dPdf = 15.0*pow(K,-4.5)*pow(fnew,6.5)-30.0*pow(K,-3.0)*pow(fnew,4.0)-3.5*pow(K,-1.5)*pow(fnew,2.5)
                 +15.0*pow(K,-1.5)*pow(fnew,1.5)+2.0;
                fnew=fold-(P)/dPdf;
                count++;
                /*printf("\n%d P=%lf f=%lf dPdf=%lf",count,P,fnew,dPdf);*/
        }
        return fnew;
}

double dist(double a1[3],double a2[3])
{
       double v[3]; //v = R(C1-C2)
       int k;
       double ip=0;
       for(k=0;k<3;k++) {
           v[k]=(a1[k]-a2[k]);
           ip+=(v[k]*v[k]);
       }
       return sqrt(ip);
}

double angl(double a1[3],double a2[3],double a3[3])
{
       double v[2][3];  //vectors of the three bonds
       int k;
       double ip00,ip11,ip01;

       ip00=ip11=ip01=0.0;

       for(k=0;k<3;k++) {
         v[0][k]=(a1[k]-a2[k]);
         v[1][k]=(a3[k]-a2[k]);
         ip00+=(v[0][k]*v[0][k]);
         ip11+=(v[1][k]*v[1][k]);
         ip01+=(v[0][k]*v[1][k]);
       }
       return acos(ip01/sqrt(ip00*ip11))*180/PI;
}

double dihe(double a1[3],double a2[3],double a3[3],double a4[3])
{      
       double v[3][3];  //vectors of the three bonds
       double n[3];
       int i,j,k;
       double ip00,ip11,ip22,ip01,ip02,ip12,ip01x2,theta;
       double x,y;

       ip00=ip11=ip22=ip01=ip02=ip12=ip01x2=0.0;

       for(k=0;k<3;k++) {
         v[0][k]=(a2[k]-a1[k]);
         v[1][k]=(a3[k]-a2[k]);
         v[2][k]=(a4[k]-a3[k]);
         ip00+=(v[0][k]*v[0][k]);
         ip11+=(v[1][k]*v[1][k]);
         ip22+=(v[2][k]*v[2][k]);
         ip01+=(v[0][k]*v[1][k]);
         ip02+=(v[0][k]*v[2][k]);
         ip12+=(v[1][k]*v[2][k]);
       }

       i=1; j=2;
       n[0]=v[i][1]*v[j][2]-v[i][2]*v[j][1];
       n[1]=v[i][2]*v[j][0]-v[i][0]*v[j][2];
       n[2]=v[i][0]*v[j][1]-v[i][1]*v[j][0];

       for(k=0;k<3;k++) ip01x2+=(v[0][k]*n[k]);

       x=-ip02/sqrt(ip00*ip22)+(ip01*ip12)/(ip11*sqrt(ip00*ip22));
       y=-ip01x2/sqrt(ip00*ip11*ip22);
       theta=acos(x/sqrt(x*x+y*y));
       if( y<0 ) theta=2*PI-theta;

       return theta*180/PI;
}

int ludcmp(double **a, int n, int *indx, double *d)
/* a[nxn]: matrix and also output of L and U matrix
   indx[1xn]: output vector that records the row permuatation effected by the partial pivoting
   d: output +1 or -1 depending on whether the number of row interchanges was even or odd
   will be needed in lubksb for matrix inversion
   */
{
    int i,imax,j,k;
    double big,dum,sum,temp;
    double *vv=new double[n];

    *d=1.0;
     for(i=0;i<n;i++) {
	big=0.0;
	for(j=0;j<n;j++) {
		if( ( temp = fabs(a[i][j]) ) > big ) big=temp;
	}
	if(big == 0.0) cerr<<"Singular matrix in routine ludcmp";
	vv[i]=1.0/big;
     }
     for(j=0;j<n;j++) {
	for(i=0;i<j;i++) {
		sum=a[i][j];
		for(k=0;k<i;k++) sum -= a[i][k]*a[k][j];
		a[i][j]=sum;
	}
	big=0.0;
	for(i=j;i<n;i++) {
		sum=a[i][j];
		for(k=0;k<j;k++) sum -= a[i][k]*a[k][j];
		a[i][j]=sum;
		if( (dum=vv[i]*fabs(sum)) >= big) {
			big=dum;
			imax=i;
		}
	}
	if(j != imax) {
		for(k=0;k<n;k++) {
			dum=a[imax][k];
			a[imax][k]=a[j][k];
			a[j][k]=dum;
		}
		*d = -(*d);
		vv[imax]=vv[j];
	}
	indx[j]=imax;
	if (a[j][j] == 0.0) a[j][j]=ludTINY;
	if (j != n) {
		dum=1.0/(a[j][j]);
		for (i=j+1;i<n;i++) a[i][j] *= dum;
	}
     }
     delete [] vv;
     return 0;
}

void lubksb(double **a, int n, int *indx, double *b)
/* LU back substitution
 *    a[nxn]: input LU decomposition of the original matrix
 *       b[1xn]: input vector
 *       */
{
     int i,ii=0,ip,j;
     double sum;

     for (i=0;i<n;i++) {
	ip=indx[i];
	sum=b[ip];
	b[ip]=b[i];
	if (ii+1)
		for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
	else if (sum) ii=i;
	b[i]=sum;
     }
     for (i=n-1;i>=0;i--) {
	sum=b[i];
	for (j=i+1;j<n;j++) sum -= a[i][j]*b[j];
	b[i]=sum/a[i][i];
     }
}

int matrix_inverse(double **A,int rank, double **invA, double &detA)
{
    int i,j,*indx;
    double d,*col;

    indx=new int[rank];
    col=new double[rank];    

    double **Ac; //copy matrix A to Ac
    new_2d_matrix(Ac,rank,rank);
    for(i=0;i<rank;i++)
	for(j=0;j<rank;j++) Ac[i][j]=A[i][j];

    //LU decomposition of Ac
    ludcmp(Ac,rank,indx,&d);
    //calculate determinant
    detA=d;
    for(i=0;i<rank;i++) detA *= Ac[i][i];

    //determine the inverse
    for(j=0;j<rank;j++) {
	for(i=0;i<rank;i++) col[i]=0.0;
	col[j]=1.0;
	lubksb(Ac,rank,indx,col);
	for(i=0;i<rank;i++) invA[i][j]=col[i];
    }

    delete_2d_matrix(Ac,rank,rank);

    delete [] indx;
    delete [] col;

    return 0;
}

int matrix_product(double **A,int ra,int ca,double **B,int rb,int cb,double **C)
{
    int i,j,k;
    if(ca!=rb) return 1;
    for(i=0;i<ra;i++) {
	for(j=0;j<cb;j++) {
		C[i][j]=0;
		for(k=0;k<ca;k++) C[i][j]+= A[i][k]*B[k][j];      
	}
    }
    return 0;
}

int matrix_product(double **A,int ra,int ca,double *B,int rb,double *C)
{
    int i,j,k;
    if(ca!=rb) return 1;
    for(i=0;i<ra;i++) {
        C[i]=0;
        for(k=0;k<ca;k++) C[i]+= A[i][k]*B[k];
    }
    return 0;
}

int matrix_BMBt(double **B,int rank,double **Mi,double **G)
{

    double **temp;
    new_2d_matrix(temp,rank,rank);
    int i,j,k;

    for(i=0;i<rank;i++) {
        for(j=0;j<rank;j++) {
                temp[i][j]=0;
                for(k=0;k<rank;k++) temp[i][j]+= B[i][k]*Mi[k][j];
        }
    }

    for(i=0;i<rank;i++) {
        for(j=0;j<rank;j++) {
                G[i][j]=0;
                for(k=0;k<rank;k++) G[i][j]+= temp[i][k]*B[j][k];
        }
    }

    delete_2d_matrix(temp,rank,rank);
    return 0;
}


int matrix_sum(double **A,double **B,int row,int col,double **C)
{
    int i,j;
    for(i=0;i<row;i++) {
	for(j=0;j<col;j++) C[i][j]=A[i][j]+B[i][j];
    }    
    return 0;    
}

int matrix_sum(double **A,double **B,int row,int col,double **C,double fac)
{
    int i,j;
    for(i=0;i<row;i++) {
        for(j=0;j<col;j++) C[i][j]=fac*(A[i][j]+B[i][j]);
    }
    return 0;
}

int matrix_sub(double **A,double **B,int row,int col,double **C)
{
    int i,j;
    for(i=0;i<row;i++) {
	for(j=0;j<col;j++) C[i][j]=A[i][j]-B[i][j];
    }    
    return 0;    
}

int matrix_scale(double **A,int row,int col,double scale)
{
    int i,j;
    for(i=0;i<row;i++) {
        for(j=0;j<col;j++) A[i][j]*=scale;
    }
    return 0;
}

int matrix_sqrt(double **A,int rank,double **Y)
{   //By Denman–Beavers iteration wikipedia
   
    double **Z,**Zi,**Yi,det;
    int i,j,k;

    new_2d_matrix(Z,rank,rank);
    new_2d_matrix(Zi,rank,rank);
    new_2d_matrix(Yi,rank,rank);

    //clean 
    for(i=0;i<rank;i++) {
       Y[i][i]=A[i][i];
       Z[i][i]=Zi[i][i]=1;
       for(j=i+1;j<rank;j++) {
           Z[i][j]=Z[j][i]=Zi[i][j]=Zi[j][i]=0;
           Y[i][j]=A[i][j];
           Y[j][i]=A[j][i];
       }
    }

    do {
      matrix_inverse(Y,rank,Yi,det);
      matrix_inverse(Z,rank,Zi,det);
      //show_2d_matrix(Y,rank,rank);
      //show_2d_matrix(Yi,rank,rank);
      //show_2d_matrix(Z,rank,rank);
      //show_2d_matrix(Zi,rank,rank);
      matrix_sum(Y,Zi,rank,rank,Y,0.5);
      matrix_sum(Z,Yi,rank,rank,Z,0.5);
      det=0;
      for(i=0;i<rank;i++) for(j=0;j<rank;j++) {
        det+=fabs(Y[i][j]-Zi[i][j]);
      }
      //cout<<det<<endl;
    } while (det>tolerance);

    delete_2d_matrix(Z,rank,rank);
    delete_2d_matrix(Zi,rank,rank);
    delete_2d_matrix(Yi,rank,rank);

    return 0;
}


int new_2d_matrix(double **&M,int row,int col)
{   
	int i;
	M= new double *[row]; 
	for(i=0;i<row;i++) M[i]=new double [col];   
	return 0;
}

int delete_2d_matrix(double **&M,int row,int col)
{   
	int i;
	for(i=0;i<row;i++) delete [] M[i];   
	delete [] M;
	return 0;
}

int show_2d_matrix(double **&M,int row,int col)
{   
	int i,j;
	for(i=0;i<row;i++) {
		for(j=0;j<col;j++) {cout.width(8); cout.precision(3); cout.flags(ios::fixed); cout<<M[i][j];}
		cout<<endl;
	}
	cout<<endl;
	return 0;
}
