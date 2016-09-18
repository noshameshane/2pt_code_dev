#ifndef FORCEFIELD_H
#define FORCEFIELD_H
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include "constant.h"

using namespace std;

class FFPREFERENCE
{
public:
   FFPREFERENCE () {
     exocyclictor=0.4;
     dielectric=1.0;
     coul_cut=8.5;
     coul_spline_off=8.5;
     coul_spline_on=8.0;
     coul_ewald_acc=0.001;
     coul_ewald_eta=2.5;
     coul_ewald_kcut=0.5;
     coul_ewald_rcut=6.0;
     coul_14_factor=1.0;  
     coul_method=EWALD;
     vdw_cut=8.5;
     vdw_spline_off=8.5;
     vdw_spline_on=8.0;
     vdw_ewald_acc=0.001;
     vdw_ewald_eta=2.5;
     vdw_ewald_kcut=0.5;
     vdw_ewald_rcut=6.0;
     vdw_14_factor=1.0;
     vdw_method=EWALD;
     vdw_combination=GEOMETRIC;
   }
   ~FFPREFERENCE () {
   }
   double exocyclictor;
   double dielectric;
   double coul_cut;
   double coul_spline_on;
   double coul_spline_off;
   double coul_ewald_acc;
   double coul_ewald_eta;
   double coul_ewald_kcut;
   double coul_ewald_rcut;
   double coul_14_factor;
   int    coul_method; 
   double vdw_cut;
   double vdw_spline_on;
   double vdw_spline_off;
   double vdw_ewald_acc;
   double vdw_ewald_eta;
   double vdw_ewald_kcut;
   double vdw_ewald_rcut;
   double vdw_14_factor;
   int    vdw_method;
   int    vdw_combination;

};

class ATOMTYPE
{
public:
   ATOMTYPE () {
      mass=chg=0;
      strcpy(name,"X");
      strcpy(fftype,"X_");
      spn=imh=lp=0;
   }
   ~ATOMTYPE () {
   }
   double mass;
   double chg;
   int    spn;  //hybridization
   int    imh;  //implicit hydrogens
   int    lp;   //electron lone pairs
   char   fftype[6];
   char   name[3];
};

class FFITEM
{
public:
   FFITEM () {
      nmax=5;
      param=NULL;
      nparam=0;
      ntorf=1;
      init_param();
      fatmtid[0]=fatmtid[1]=fatmtid[2]=fatmtid[3]=-1;
   }
   ~FFITEM () {
      if(param!=NULL) delete [] param;
   }
   int    nbody;
   int    fatmtid[4]; //ATOMTYPE ID
   char   eetype[25]; //energy expression type
   int    nparam;     //number of parameters
   int    nmax;       //number of entries in item list (nmax=nparam*ntorf)
   int    ntorf;      //number of fourier terms for torsion
   double *param;     //ff parameter
   int    init_param();
};

class FORCEFIELD : public FFPREFERENCE
{
public:
   FORCEFIELD () {
     version=0;
     natmt=nvdw=noffvdw=nbond=nangle=ntor=ninv=nimp=0;
     atmt=NULL;  
     vdw=NULL;  
     offvdw=NULL;  
     bond=NULL;  
     angle=NULL;  
     tor=NULL;  
     inv=NULL;  
     imp=NULL;  
     vdwij=NULL;  
   }
   ~FORCEFIELD () {
     if(atmt   !=NULL) delete [] atmt; 
     if(vdw    !=NULL) delete [] vdw; 
     if(offvdw !=NULL) delete [] offvdw; 
     if(bond   !=NULL) delete [] bond; 
     if(angle  !=NULL) delete [] angle; 
     if(tor    !=NULL) delete [] tor; 
     if(inv    !=NULL) delete [] inv; 
     if(imp    !=NULL) delete [] imp; 
     if(vdwij  !=NULL) { int i;
        for(i=0;i<nvdw;i++) delete [] vdwij[i];
        delete [] vdwij; 
     }
   }

   ATOMTYPE *atmt;  //atom types
   FFITEM *vdw;     //van der Waals
   FFITEM *offvdw;  //user defined off-diagonal van der Waals
   FFITEM **vdwij;  //van der Waals interaction matrix
   FFITEM *bond;    //bond
   FFITEM *angle;   //angle
   FFITEM *tor;     //torsion
   FFITEM *inv;     //inversion
   FFITEM *imp;     //improper torsion
   int natmt;
   int nvdw;
   int noffvdw;
   int nbond;
   int nangle;
   int ntor;
   int ninv;
   int nimp;
   int init_ff(FFITEM **,int,int);
   int init_atmt();
   int init_vdwij();  //initialize the vdW interaction matrix
   int ftype2id(char *);

   char fffile[1024];
   int  version;

};

#endif
