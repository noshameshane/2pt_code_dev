#ifndef MODEL_H
#define MODEL_H
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector>
#include "timing.h"
#include "constant.h"
#include "structure.h"
#include "property.h"
#include "forcefield.h"
#include "trajectory.h"
#include "control.h"
#include "statistics.h"

using namespace std;

class MODEL
{
public:
   MODEL () {
      int i;
      nmol=nbond=nangle=ntor=ninv=nimp=ngrp=0;
      atom=NULL;
      bond=NULL;
      angle=NULL;
      tor=NULL;
      inv=NULL;
      imp=NULL;
      mol=NULL;
      grp=NULL;
      periodic=0;
      xyzfreq=velfreq=engfreq=strfreq=0;
      timestep=0;
      clrmax=clrmin=0;
   }
   ~MODEL () {
      int i;
      if(atom!=NULL)  delete [] atom;
      if(bond!=NULL)  delete [] bond;
      if(angle!=NULL) delete [] angle;
      if(tor!=NULL)   delete [] tor;
      if(inv!=NULL)   delete [] inv;
      if(imp!=NULL)   delete [] imp;
      if(mol!=NULL)   delete [] mol;
      if(grp!=NULL)   delete [] grp;
   }
   int nmol;    //number of molecules
   int natom;   //number of atoms
   int nbond;   //number of bonds
   int nangle;  //number of angles
   int ntor;    //number of torsions
   int ninv;    //number of inversions
   int nimp;    //number of improper torsions
   int ngrp;    //number of groups

   int periodic;  //periodicity
   CONTROL ctl;       //control file
   CELL cell;         //cell properties
   TRAJECTORY trj;    //trajectory
   PROPERTY prp;      //model properties
   FORCEFIELD ff;     //force field settings
   TIME timing;       //timing

   char   name[1024]; //name of the model

   ATOM *atom;
   BOND *bond;
   ANGLE *angle;
   TORSION *tor;
   INVERSION *inv;
   IMPROPER *imp;
   MOLECULE *mol;
   GROUP *grp;

   int xyzfreq;  //coordinate dump frequency
   int velfreq;  //velocity dump frequency
   int engfreq;  //energy dump frequency
   int strfreq;  //stress dump frequency
   double timestep; //simulation time step in ps
   double clrmax;  //max scale in coloring
   double clrmin;  //max scale in coloring

   int  init_atom();
   int  init_bond();
   int  init_angle();
   int  init_tor();
   int  init_inv();
   int  init_imp();
   int  init_mol();
   int  init_grp();

   int rd_strt();
   int rd_ff();
   int rd_lmpdata_old(char *); //read lammps data file (fortran version)
   int rd_lmpdata(char *); //read lammps data file (c++ version)
   int rd_gro(char *);     //read gromacs data file (c++ version)
   int rd_g96(char *);     //read gromacs g96 data file (c++ version)
   int rd_bgf(char *);     //read bgf file
   int rd_cpmdin(char *);  //read cpmd in file
   int rd_grp(char *);

   int out_grp(ostream *);
   int out_bgf(ostream *); //output bgf file
   int out_g96(ostream *); //output gro file
   int out_lmpdata(ostream *); //output LAMMPS data file (c++ version)
   int out_cpmdin(ostream *); //output CPMD in file
   int out_structure(ostream *); //output molecular structure

   int ffsetup();     //find ffid 
   int bond2mol();    //find molecules from connected bonds
   int bond2cnt();    //build connectivity from bonding info
   int cnt2valence(); //build bond, angle, torsion, inversion from connectivity
   int cnt2bod();     //build bond order from connectivity info
   int dst2bond();    //find connected bonds from atomic distances
   int mol2bond();    //find connected bonds from atomic distances within a molecule
   int farneigh();    //find farthermost connected atom
   int headtail();    //find head and tail atoms for each molecule in the model
   int fixlongbond(); //fix long bonds (by translation of atoms) caused by chemical rxn
   int findlongbond(); //find long bonds
   int setgrids4atm(); //set grids for atom, for void analysis
   int find_void();    //find voids, for void analysis
   int ck_exocyclic(int,int); //check for exocyclic torsion

   int cvtmol2grp();  //create group file for each molecule
   int cvteqmass2grp();  //create group file for atoms equivalent masses
   int cvteqmolmass2grp();  //create group file for molecules of equivalent masses
   int mass2types(); //check existance of atom name and types, if nonexistant fill in with default values
   int element2prp(); //element to properties
   int ff2mass(); //force field type to mass
   int cal_mass(); //calc total mass of the model
   int cal_grp_mass(); //calc mass of each group
   int ck_range(); //check atom range for nonperiodic system;
   int cal_multipole(); //calc multipole moments of the model
   int cal_vcomp(); //calc velocity components
   int cal_T();  //calculate temperature of model using atomic velocities
   int cal_P();  //calculate pressure of model from the virial tensor strs[]
   int find_molcm(); //find the center of mass position of each molecule
   int find_cm(); //find the center of mass position of the model
   int find_molingrp(); //find molecules in each group
   int setdspcolor(int); //set atom color
   int setdsplabel(int); //set atom label
   int setdspatomsize(double); //set atom display size
   int setdspbondsize(double); //set bond display size
   int setdsplabelsize(double); //set label display size
   int setgrpdisplay(); //set group display

   //model manipulation
   int mapatom(int); //fold atoms to unit cell
   int setcenter(double *); //set the center of box to 
   int translate(double *); //translate model position
   int superlattice(); //translate model position
   int init_trj();

};

#endif
