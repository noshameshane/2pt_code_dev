#ifndef STRUCTURE_H
#define STRUCTURE_H
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <iostream>
#include "constant.h"

using namespace std;

class CELL
{
public:
   CELL () {
     H[0][0]=H[0][1]=H[0][2]=0;
     H[1][0]=H[1][1]=H[1][2]=0;
     H[2][0]=H[2][1]=H[2][2]=0;
     la=lb=lc=alpha=beta=gamma=volume=0;
     a[0]=a[1]=a[2]=0;
     b[0]=b[1]=b[2]=0;
     c[0]=c[1]=c[2]=0;
     o[0]=o[1]=o[2]=0;
     Hinv[0][0]=Hinv[0][1]=Hinv[0][2]=0;
     Hinv[1][0]=Hinv[1][1]=Hinv[1][2]=0;
     Hinv[2][0]=Hinv[2][1]=Hinv[2][2]=0;
     cb[0]=cb[1]=cb[2]=0;
     sgrid=0;
     ngrid[0]=ngrid[1]=ngrid[2]=0;
     vgrid=NULL;
   }
   ~CELL () {
      free_grid();
   }
   double H[3][3];
   double Hinv[3][3];
   double a[3];  //vector of a
   double b[3];  //vector of b
   double c[3];  //vector of c
   double la;    //length of a
   double lb;    //length of b
   double lc;    //length of c
   double alpha;
   double beta;
   double gamma;
   double volume;
   double o[3];  //origin
   double cb[3]; //center of box
   int s2r2H(double *);
   int H2others();
   int labc2H();
   int cal_Hinv();
   int cal_volume();
   int Map2UnitCell(double*ori,double*shift);
   int mindst(double*ori,double*pv);
   int out_cell(ostream *);
   double sgrid;    //grid size
   int    ngrid[3]; //number of grids in each dimension
   char   ***vgrid; //grid for void analysis
   int    set_grid(double);
   int    free_grid();
};

class AGRID
{
public:
   AGRID() { 
      sgrid=0;
      ngrid=0;
      vgrid=NULL;
   }
   ~AGRID() { 
      free_grid();
   }
   double radius;   //atom radius
   double sgrid;    //grid size
   int    ngrid;    //number of grids in each dimension
   char   ***vgrid; //grid for void analysis
   int    set_grid(double);
   int    free_grid();
};

class ATOM
{
public:
   ATOM () {
      mass=chg=eng=0;
      ffid=mol=id=-1;
      ncnt=0;
      int i;
      for(i=0;i<maxcnt;i++) { bod[i]=1; cnt[i]=-1; }
      faratm=-1;
      fardst=0;
      strcpy(name,"X");
      strcpy(molname,"X");
      strcpy(fftype,"X_");
      strcpy(label,"X");
      gd=NULL;
      color=-1;
      atomsize=bondsize=labelsize=1;          //display size
      radius=-1;
      crad=-1;
      shift[0]=shift[1]=shift[2]=0;
      nxcl13=nxcl14=nlist=0;
      xcl13=NULL;
      xcl14=NULL;
      list=NULL;
   }
   ~ATOM () {
      gd=NULL;
      if(xcl13!=NULL) delete [] xcl13;
      if(xcl14!=NULL) delete [] xcl14;
      if(list !=NULL) delete [] list;
   }
   double f[FORCETYPE][3];  //net forces (bond,angle,torsion,inversion,vdw,vdw corr,elec,ewaldsum,net) in kcal/mol-A
   double pv[6];    //position vector in angstroms, pv[4-6] stores pv of previous step, used for Verlet
   double vel[3];   //velocity in A/ps
   double vt[3];    //translational velocity
   double vr[3];    //rotational velocity
   double vv[3];    //vibrational velocity
   double shift[3]; //shift to base unit cell, used in atom remapping
   double mass;     //mass in g/mol
   double chg;      //charge in electrons
   double eng;      //atomic energy in kcal/mol
   int    ffid;     //ATOMTYPE
   char   fftype[6];//forcefild type
   int    mol;      //molecule id
   int    id;       //atom id
   int    ncnt;     //number of atoms connected to it
   int    cnt[maxcnt];   //ids of atoms connected to it
   int    bod[maxcnt];   //bond order
   int    nxcl13;   //number of 13 exclusion
   int    *xcl13;   //list of 13 exclusion
   int    nxcl14;   //number of 14 exclusion
   int    *xcl14;   //list of 14 exclusion
   int    nlist;    //number of items in the nonbond list
   int    *list;    //nonbond list, created in neighbor.cpp
   char   name[3];  //atom name
   char   molname[6];//molecule name
   int    faratm;   //id of the farthermost connect atom
   int    fardst;   //number of atoms away from faratm
   int    color;    //atom color for visualization
   char   label[50];//atom label for visualization
   double atomsize; //display atom size
   double bondsize; //display bond size
   double labelsize; //display label size
   double radius;   //atomic radius
   double crad;     //covalent radius for covalent bond calculation
   AGRID  *gd;      //for void analysis
private:
};

class MOLECULE
{
public:
   MOLECULE () {
      atm=NULL;
      natom=0;
      mass=chg=eng=mu[0]=mu[1]=mu[2]=0;
      head=tail=-1;
      linear=0;
      atom=NULL;
      inertia=NULL;
      inertia=new double *[3];
      int i;
      for(i=0;i<3;i++) {
             inertia[i]=new double [3];
      }
      iccnt=NULL;
      mqvib=NULL;
      moi=NULL;
   }
   ~MOLECULE () {
      if(atm!=NULL) delete [] atm;
      if(atom!=NULL) {
         delete [] atom;
         atom=NULL;
      }
      if(inertia!=NULL) {
         delete [] inertia[0];
         delete [] inertia[1];
         delete [] inertia[2];
         delete [] inertia;
      }
      if(iccnt!=NULL) delete [] iccnt;
      if(mqvib!=NULL) delete [] mqvib;
      if(moi!=NULL) delete [] moi;
   }
   int natom;
   int *atm;
   int init_atom();
   double mass;     //mass
   double chg;      //charge
   double eng;      //molecular energy in kcal/mol
   int head;
   int tail;
   int linear;      //linear molecule
   double pv[3];    //position vector (center of mass) in angstroms
   double vel[3];   //center of mass velocity (A/ps)
   double mu[3];    //dipole
   double **inertia;  //moment of inertia tensor
   double omega[3]; //angular velocity
   double angmom[3];//angular momentum
   double anguv[3]; //angular velocity
   double pI[3];    //principle moment of inertia (kg*m2)
   int    (*iccnt)[5]; //connectivity for internal coordinates
   double *mqvib;    //mass weighted velocity of each internal degree of freedom (3*natom) 
   double *moi;      //moment of inertia of each internal degree of freedom (3*natom) 
   ATOM   **atom;
   int cal_mass();   //calculate center of mass 
   int cal_chgmu();  //calculate charge and dipole
   int cal_vc();     //calculate velocity components
   int cal_vc1();    //calculate velocity components, first version
   int cal_Zmat(double(*)[3],double **);   //calculate Zmat from iccnt
   int find_cm();    //calculate center of mass pv[3]
   int ini_iccnt();  //initialization for iccnet[natom][5]
};

class BOND
{
public:
   BOND () {
      ffid=id=-1;
      lbd=-1;
      len=-1;
   }
   ~BOND () {
   }
   int atm[2];    //atom pairs
   double len;      //current bond length angstroms
   int    id;       //bond id
   int    ffid;     //id in forcefield
   int    lbd;      //bond when folded into one cell (for proton transfer)
   ATOM   *atom[2];
   int    cal_len();
private:
};

class ANGLE
{
public:
   ANGLE () {
      ffid=id=-1;
      theta=theta0=costh=0;
   }
   ~ANGLE () {
   }
   int atm[3];      //angle atoms
   double theta0;   //equilibrium angle in radians
   double theta;    //current angle in radians
   double costh;    //cos(theta)
   int    id;       //angle id
   int    ffid;
   ATOM   *atom[3];
   int    cal_theta();
private:
};

class TORSION
{
public:
   TORSION () {
      ffid=id=-1;
      deg=1;
      theta=theta0=costh=0;
      exo=1;
   }
   ~TORSION () {
   }
   int atm[4];      //torsion atoms
   double theta0;   //equilibrium torsional angle in radians
   double theta;    //current torsional angle in radians
   double costh;    //cos(theta)
   double exo;      //exocyclic factor
   int    id;       //torsion id
   int    ffid;
   int    deg;      //degeneracy 
   ATOM   *atom[4];
   int    cal_theta();
   int    cal_theta1();
private:
};

class INVERSION
{
public:
   INVERSION () {
      ffid=id=-1;
      theta=theta0=costh=0;
   }
   ~INVERSION () {
   }
   int atm[4];      //inversion atoms
   double theta0;   //equilibrium inversion angle in radians
   double theta;    //current inversion angle in radians
   double costh;    //cos(theta)
   int    id;       //inversion id
   int    ffid;
   ATOM   *atom[4];
   int    cal_theta();
private:
};

class IMPROPER
{
public:
   IMPROPER () {
      ffid=id=-1;
   }
   ~IMPROPER () {
   }
   int atm[4];      //improper torsion atoms
   double theta0;   //equilibrium improper torsional angle in radians
   double theta;    //current improper torsional angle in radians
   int    id;       //improper id
   int    ffid;
   ATOM   *atom[4];
private:
};

class GROUP
{
public:
  GROUP() {
     natom=0;
     nmol=0;
     constraint=0;  //fixed bonds/angles
     rotsym=1;      //rotational symmetry
     linear=0;      //linear molecule
     pmvr=1;        //partial molar volume ratio
     chg=mass=eng=0;
     atm=NULL;
     mol=NULL;
     color=-1;
     atomsize=bondsize=labelsize=1;
     mu[0]=mu[1]=mu[2]=0;
  }
  ~GROUP() {
     natom=0;
     if(atm!=NULL) delete [] atm;
     if(mol!=NULL) delete [] mol;
  }
  int natom;
  int nmol;
  int constraint;  //number of constraints
  int rotsym;      //rotational symmetry
  int linear;      //linear molecule
  int color;
  double pmvr;     //partial molar volume ratio
  double atomsize; 
  double bondsize;
  double labelsize;
  double chg;  //electrons
  double mu[3];
  double mass; //g/mol
  double eng;  //kcal/mol
  int *atm;
  int *mol;
  int init_atom();
};

class SBOND
{ //special bonds, currently for displaying purpose only
public:
   SBOND () {
   }
   ~SBOND () {
   }
   ATOM   atom[2];
   double mpv[2][3]; //mid point position, for visualization
private:
};

class HBOND
{ //Hbond
public:
   HBOND () {
      oodist=coshoo=0;
      hodist=cosoho=0;
      hov[0]=hov[1]=hov[2]=0;
      cosdih=0;
      bc[0]=bc[1]=bc[2]=0;
   }
   ~HBOND () {
   }
   int atm[3];      //h-bond atoms O[2]...H[0]-O[1]
   double hov[3];   //H[0]-O[2] vector
   double oodist;
   double hodist;
   double coshoo;
   double cosoho;
   double cosdih;   //dihedral angle between the two h-bonded water molecules
   int    bc[3];    //recording the crossing of unit cell boundary
   double mpv[2][3]; //mid point vector, for visualization
private:
};

class RING
{
public:
   RING() {
       size=0;
       atm=NULL;
   }
   ~RING() {
       if(atm!=NULL) delete [] atm;
   }
   int size;
   int *atm;
   int init_atm();
private:
};


#endif
