#include "forcefield.h"

int FFITEM::init_param()
{
    if(nmax==0) return 1;
    if(param!=NULL) delete [] param;
    param=new double [nmax];
    return 0;
}

int FORCEFIELD::init_atmt()
{
    if(natmt==0) return 1;
    if(atmt!=NULL) delete [] atmt;
    atmt=new ATOMTYPE [natmt];
    return 0;
}

int FORCEFIELD::init_ff(FFITEM **item, int number, int nbody)
{
    if(number==0) return 1;
    if(*item!=NULL) delete [] item;
    *item=new FFITEM [number];
    int i;
    for(i=0;i<number;i++) (*item)[i].nbody=nbody;
    return 0;
}

int FORCEFIELD::init_vdwij()
{
    int i,j,k,l;
    if(vdwij!=NULL) {
       for(i=0;i<nvdw;i++) delete [] vdwij[i];
       delete [] vdwij; vdwij=NULL;
    }
    if(nvdw==0) return 1;
    vdwij=new FFITEM *[nvdw];
    for(i=0;i<nvdw;i++) vdwij[i]=new FFITEM [nvdw];

    //on diagonal vdw
    for(i=0;i<nvdw;i++) {
       vdwij[i][i].nbody=2;
       vdwij[i][i].fatmtid[0]=vdwij[i][i].fatmtid[1]=vdw[i].fatmtid[0];
       strcmp(vdwij[i][i].eetype,vdw[i].eetype);
       vdwij[i][i].nparam=vdw[i].nparam;
       vdwij[i][i].nmax=vdw[i].nmax;
       vdwij[i][i].ntorf=vdw[i].ntorf;
       vdwij[i][i].init_param();
       for(k=0;k<vdwij[i][i].nmax;k++) vdwij[i][i].param[k]=vdw[i].param[k];
    }
    //establish default off diagonal vdw
    for(i=0;i<nvdw;i++) {
       for(j=i+1;j<nvdw;j++) {
           vdwij[i][j].nbody=vdwij[j][i].nbody=vdw[i].nbody;
           vdwij[i][j].fatmtid[0]=i; vdwij[i][j].fatmtid[1]=j;
           vdwij[j][i].fatmtid[0]=j; vdwij[j][i].fatmtid[1]=i;
           strcmp(vdwij[i][j].eetype,vdw[i].eetype);
           strcmp(vdwij[j][i].eetype,vdw[i].eetype);
           vdwij[i][j].nparam=vdwij[j][i].nparam=vdw[i].nparam;
           vdwij[i][j].nmax=vdwij[j][i].nmax=vdw[i].nmax;
           vdwij[i][j].ntorf=vdwij[j][i].ntorf=vdw[i].ntorf;
           vdwij[i][j].init_param(); vdwij[j][i].init_param();
           if(vdw_combination==GEOMETRIC) {
             vdwij[i][j].param[0]=sqrt(vdw[i].param[0]*vdw[j].param[0]); ////geometric mean for Ro
             vdwij[i][j].param[1]=sqrt(vdw[i].param[1]*vdw[j].param[1]); //geometric mean for Do
           } else if (vdw_combination==ARITHMETIC) {
             vdwij[i][j].param[0]=0.5*(vdw[i].param[0]+vdw[j].param[0]); //arithmetic mean for Ro
             vdwij[i][j].param[1]=sqrt(vdw[i].param[1]*vdw[j].param[1]); //geometric mean for Do
           } else if (vdw_combination==SIXTHPOWER) {
             vdwij[i][j].param[0]=pow(0.5*(pow(vdw[i].param[0],6.0)+pow(vdw[j].param[0],6.0)),1.0/6.0);
             vdwij[i][j].param[1]=(2.0*sqrt(vdw[i].param[1]*vdw[j].param[1])*pow(vdw[i].param[0],3.0)*pow(vdw[i].param[0],3.0))/(pow(vdw[i].param[0],6.0)+pow(vdw[j].param[0],6.0)); 
           }
           if(vdw[i].nparam==3&&vdw[j].nparam==3) { //exponential-6
              vdwij[i][j].param[2]=sqrt(vdw[i].param[2]*vdw[j].param[2]); //geometric mean for Do
           }
           for(k=0;k<vdwij[i][j].nmax;k++) vdwij[j][i].param[k]=vdwij[i][j].param[k];
       }
    }
    //fill in off diagonal vdw if available
    for(k=0;k<noffvdw;k++) {
       i=offvdw[k].fatmtid[0];
       j=offvdw[k].fatmtid[1];
       vdwij[i][j].nbody=vdwij[j][i].nbody=offvdw[k].nbody;
       vdwij[i][j].fatmtid[0]=i; vdwij[i][j].fatmtid[1]=j;
       vdwij[j][i].fatmtid[0]=j; vdwij[j][i].fatmtid[1]=i;
       strcmp(vdwij[i][j].eetype,offvdw[k].eetype);
       strcmp(vdwij[j][i].eetype,offvdw[k].eetype);
       vdwij[i][j].nparam=vdwij[j][i].nparam=offvdw[k].nparam;
       vdwij[i][j].nmax=vdwij[j][i].nmax=offvdw[k].nmax;
       vdwij[i][j].ntorf=vdwij[j][i].ntorf=offvdw[k].ntorf;
       vdwij[i][j].init_param(); vdwij[j][i].init_param();
       for(l=0;l<offvdw[k].nmax;l++) vdwij[i][j].param[l]=vdwij[j][i].param[l]=offvdw[k].param[l];
    }

    if(0) {
       for(i=0;i<nvdw;i++) {
          for(j=0;j<nvdw;j++) {
             printf("%d %d %8f/%8f ",i,j,vdwij[i][j].param[0],vdwij[i][j].param[1]);
          }
          printf("\n");
       }
    }
    //param[0] Ro in A
    //param[1] Do in kcal/mol

    //change param[0] to 12*Do*Ro^12, unit kcal/mol-A12
    //change param[1] to 12*Do*Ro^6, unit kcal/mol-A6
    double v0,v1;
    for(i=0;i<nvdw;i++) {
          v0=vdwij[i][i].param[0]; v1=vdwij[i][i].param[1];
          vdwij[i][i].param[0]=12*v1*pow(v0,12.0);
          vdwij[i][i].param[1]=12*v1*pow(v0,6.0);
       for(j=i+1;j<nvdw;j++) {
          v0=vdwij[i][j].param[0]; v1=vdwij[i][j].param[1];
          vdwij[i][j].param[0]=vdwij[j][i].param[0]=12*v1*pow(v0,12.0);
          vdwij[i][j].param[1]=vdwij[j][i].param[1]=12*v1*pow(v0,6.0);
       }
    }

    if(0) { //print the vdwij matrix
       for(i=0;i<nvdw;i++) {
          for(j=0;j<nvdw;j++) {
             printf("%d %d %8f/%8f ",i,j,vdwij[i][j].param[0],vdwij[i][j].param[1]);
          }
          printf("\n");
       }
    }

    return 0;
}

int FORCEFIELD::ftype2id(char *type)
{
    int id;
    id=natmt-1;
    while(strcmp(type,atmt[id].fftype)!=0) {
       id--;
       if(id<0) break;
    }
    
    return id;
}

