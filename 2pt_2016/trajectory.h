#ifndef TRAJECTORY_H
#define TRAJECTORY_H
#include <iostream>
#include <fstream>
#include <cstring>
#include <iomanip>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include "constant.h"
#include "property.h"
#include "structure.h"

using namespace std;

const int MAXFILE=24;
typedef int   logical;

class TRJheader
{
public:
    TRJheader() {
	mvatmofst=NULL;
    }
    ~TRJheader() {
        if(mvatmofst!=NULL) delete [] mvatmofst;
    }
    char hdr[5];	//Header
    int icntrl[20];	//Control (icntrl[0]=trj file version, default=2010)
    int ntrjti;		//number of comments
    char trjtic[10][80];	//comments
    int neexti;		//number of EEX comments
    char eextic[10][80];	//EEX comments
    logical period;		//Periodicity (1,2,3 =1D,2D,3D)
    logical molxtl;		//Molecular Crystal(unknown usage)
    logical lcanon;		//Canonical (NVT)
    logical defcel;		//Cell definition (make it true to be safe)
    logical prtthrm;	//Perturbation Theory
    logical lnose;		//NoseorHoover
    logical lnpecan;	//NPT Canonical
    logical ltmpdamp;	//Temperature damping
    int  nflusd;		//number of files unsed (make it 1)
    int  mvatmpfu[MAXFILE];	//Moveable atoms (must give mvatmpfu[0])
    int  natmpfu[MAXFILE];	//Total atoms (must give natmpfu[0])
    char  decusd[MAXFILE][8];//Descriptor
    int  totmov;		//Total moveable atoms (totmov=mvatmpfu[0] when nflusd=1)
    int  *mvatmofst; 	//movable atom ids (starts from 1)
    int leexti;		//length of EEX title
    char eextit[80];	//EEX title
    int lparti;		//length of parameter file
    char partit[80];	//parameter file
    int natom;		//total number of atoms		
    int nmovatm;
    int movatm1;
    int movatmn;
    int version;
    
    int xyzfreq,nxyz;  //coordinate dump frequency, skips when reading traj
    int velfreq,nvel;  //velocity dump frequency
    int engfreq,neng;  //energy dump frequency
    int strfreq,nstr;  //stress dump frequency
    double timestep; //simulation time step in ps
    int totaccustep;   //total accumulated steps (needed for CPMD trj)

    int init_header();
    int CleanTRJheader();
    int ReadBinHeader(ifstream *intrj);
    int ReadAscHeader(ifstream *intrj);
    int ReadCPMDHeader(ifstream *intrj,ifstream *ineng,ifstream *instr);
    int ReadLMPHeader(ifstream *intrj,ifstream *ineng);
    int ReadUSRHeader(ifstream *intrj);
    int WriteASCIIheader(ostream *, int oformat, int ndel, int iddel[]);
    int WriteBinaryheader(ostream *, int oformat);

private:
};

class TRJ //double precision for version 2010,2012
{
public:
    TRJ() {
           itrjformat=0;
           otrjformat=0;
           ascfid=datasize2fid=0;
           location0e=datasize2fide=0;
           location0s=datasize2fids=0;
    }
    ~TRJ() {
           if(intrj.is_open()) intrj.close();
    }

    CELL *cell;
    PROPERTY *prp;
    ATOM *atom;

    TRJheader header;
        
    //double curtim;	//current time (ps)
    //int itstep;	//step
    //double system;	//? system temp?
    double avetem;	//
    double dstep;	//time step
    double firstt;	//initial set temperature K
    double finalt;	//final set temperature
    //double  e;		//potential energy kcal/mol
    //double  eb;		//bond energy
    //double  et;		//angle
    //double  ep;		//torsion
    //double  ei;		//inversion
    //double  enb;	//van der Waals
    //double  eel;	//Columb 
    //double  ehb;	//Hydrogen bond
    double  ec;		//constraint
    double  eu;		//user energy
    double  tint;	//total interanl
    double  tnb;	//total nonbond
    double  ea;		//trj average properties
    double  eba;
    double  eta;
    double  epa;
    double  eia;
    double  enba;
    double  eela;
    double  ehba;
    double  eca;
    double  eua;
    double  tinta;
    double  tnba;	//trj average properties
    //double  tote;	//total energy
    //double  totke;	//kinetic energy
    double  totea;
    double  tkea;
    double	dum[12];	//unknown meanings
    double	duma[12];	//nuknown averaged properties
    int  iconmp;
    int imstep;
    logical lvelwr;	//velocity output
    logical lfrcwr;	//force output
    logical lengwr;	//atomic energy output
    int iconfs;
    int icstep;

    double  pressur; //, pressura;	//pressure in GPa
    double  vol; //, vola;		//volume A3
    double  pvtot, pvtota;
    double  pvkin, pvkina;
    double  pvpot, pvpota;
    double  radgyr, radgyra;

    double  signose;
    double  zfrict;
    double	zprfrict;
    double  snose,snoseh,ssdot;
    double	qcanon;
    double	sigdyn[2];
    double	gamtmp;

    double	tcel,ucel,tcela,ucela;
    double	s2ra[6];//,s2r[6];
    double	s2rdot[6];
    
    int natmcel;
    double  strsa[6];//,strs[6]; xx,yy,zz,xy,xz,yz
    double  extstrs,extstrsa; //! added after 220 but version called 210

    double  eabtota;
    double  eabvala;
    double  eabelha;
    double  eabnba;
    double  eabmisa;
    double  dltfaba;
    double  expprta;
    //double *x,*y,*z; //coordinates in Angstroms
    //double *velx,*vely,*velz; //velocities in Angstroms/ps

    int itrjformat;     //in trj format 0:binary 1:ascii
    int otrjformat;     //out trj format 0:binary 1:ascii
    ifstream intrj;     //input trj stream
    ifstream ineng;     //input energy stream
    ifstream instr;     //input stress stream
    unsigned long datasize,fcount,totframe;
    unsigned long location0,location1,location2;//after header,after 1st frame,end
    unsigned long ascfid,datasize2fid; //frame id, datasize upto frame fid
    unsigned long location0e,datasize2fide; //after header, datasize upto frame lf
    unsigned long location0s,datasize2fids; //after header, datasize upto frame lf

    int init_trj(char *);
    int init_c2bintrj(char *);
    int init_c2asctrj(char *);
    int init_cpmdtrj(char *); //cpmd
    int init_lmptrj(char *);  //lammps
    int init_usrtrj(char *);  //lammps
    int rd_frame(int ); //read frame
    int wt_frame(ostream *); //write current content to frame
    int rd_header();    //read header
    int wt_header(ostream *); //write header
    int init_contentDouble();
    int CleanTRJcontentDouble();
    int ReadBinEng();
    int ReadBinFrame();
    int ReadAscFrame();
    int ReadCPMDFrame();
    int ReadLMPFrame();
    int ReadUSRFrame();
    int WriteASCIIcontent(ostream *, int oformat, int ndel, int iddel[]);
    int WriteBinarycontent(ostream *, int oformat);

private:
};

class TRAJECTORY //multiple cerius2 trj
{
public:
      TRAJECTORY() {
         strj=NULL;
         iframe=NULL;
         fframe=NULL;
      }
      ~TRAJECTORY() {
          if(strj!=NULL) delete [] strj;
          if(iframe!=NULL) delete [] iframe;
          if(fframe!=NULL) delete [] fframe;
      }
       TRJ *strj;
       int itrjf; //in trj format
       int otrjf; //out trj format
       int ntrjf; //number of trj files
       int ctrjf;  //current trj file
       int cframe; //current frame id
       int *iframe;      
       int *fframe;      
       int totframe;
       int init_trj(int, char (*)[1024],int,int,CELL*,ATOM*,PROPERTY*);
       int rd_frame(int );
       int wt_frame(ostream *);
       int wt_header(ostream *);
private:
};

#endif
