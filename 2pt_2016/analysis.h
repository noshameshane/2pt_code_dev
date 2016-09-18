#ifndef ANALYSIS_H
#define ANALYSIS_H
#include <iostream>
#include <fstream>
#include "/ul/stlin/local/include/fftw3.h"
#include <cstdlib>
#include <cmath>
#include "constant.h"
#include "utility.h"
#include "statistics.h"
#include "model.h"

using namespace std;

class ANALYSIS
{
public:
   ANALYSIS () {  
      nframe=0;
      Pvac=NULL;
      PvacInt=NULL;
      vacatom=NULL;
      vacvv=NULL;
      vacT=NULL;
      vacDF=NULL;
      vacTInt=NULL;
      vacDFInt=NULL;
      pI=NULL;
      vacread=NULL;
   }
   ~ANALYSIS () {
      int i,j,k;
      if(Pvac!=NULL) {
        for(i=0;i<model->ngrp;i++) {
           for(j=0;j<vactype;j++) delete [] Pvac[i][j];
           delete [] Pvac[i];
         }
         delete [] Pvac;
      }
      if(PvacInt!=NULL) {
        for(i=0;i<model->ngrp;i++) {
           int df=3*model->mol[(model->grp[i].mol[0])].natom-6;
           for(j=0;j<df;j++) delete [] PvacInt[i][j];
           delete [] PvacInt[i];
         }
         delete [] PvacInt;
      }
      if(vacvv!=NULL) {
         for(k=0;k<vactype;k++) {
            for(i=0;i<vacnsteps;i++) {
               for(j=0;j<vacatmperpass;j++) delete [] vacvv[k][i][j];
               delete [] vacvv[k][i];
            }
            delete [] vacvv[k];
         }
         delete [] vacvv;
      }
      if(vacatom!=NULL) delete [] vacatom;
      if(vacT!=NULL) {
         for(i=0;i<model->ngrp;i++) delete [] vacT[i];
         delete [] vacT;
      }
      if(vacDF!=NULL) {
         for(i=0;i<model->ngrp;i++) delete [] vacDF[i];
         delete [] vacDF;
      }
      if(vacTInt!=NULL) {
         for(i=0;i<model->ngrp;i++) delete [] vacT[i];
         delete [] vacT;
      }
      if(vacDFInt!=NULL) {
         for(i=0;i<model->ngrp;i++) delete [] vacDF[i];
         delete [] vacDF;
      }
      if(pI!=NULL) delete [] pI;
      if(vacread!=NULL) delete [] vacread;
   }
   MODEL *model;

   int nframe;

   int ana_init();
   int ana_ana();
   int ana_report(ostream *); //output analysis results

   ofstream trjf;      //output trj 
   int ini_trj();
   int out_trj();
   int ini_fix_gro_trj();  //initialization for fix_gro_trj()
   int fix_gro_trj();  //fix the "no jump" issue in gromacs trj
   int *gtrjchk;
   double (*gtrjfix)[3];

   char outname[1024];
   int out_bgf();
   int out_lmp();

   STAT ***Pvac;   //Pvac[ngrp][vtype][nstep]
   STAT ***PvacInt;   //Pvac[ngrp][3*natom-6][nstep]
   ATOM *vacatom,*tatom;
   CELL vaccell;
   PROPERTY vacprp;
   int  vacmaxf,vacnsteps; //correlation steps, number of origins
   int  c2pt; //2pt correction type
   int  cmol; //consider molecules
   int  c2ptmf; //consider menory function correction
   int  vactype; //types of velocities
   int  *vacread; //atom read list
   double **vacT;  //average temperature  //vacT[ngrp][vtype]
   double **vacDF; //degrees of freedom   //vacDF[ngrp][vtype]
   double **vacTInt;  //average temperature  //vacT[ngrp][3*natom-6]
   double **vacDFInt; //degrees of freedom   //vacDF[ngrp][3*natom-6]
   double vacE;  //average energy
   double vacV;  //average volume
   double vacP;  //average pressure
   double trjT;  //average temperature from trj file
   double vacdtime;
   double (*pI)[3];//principle moments of inertia
   float ****vacvv; //vacvv[vtype][step][atom][xyz]
   int  vacatmperpass,vacpassatms,vacatmend,vacatmpassed;
   int  vacnpass,vacipass; //number of passes for reading the whole trj
   int  vactatmpassed,vactmolpassed,vacatmpass,vacmolpass,vaclastmol; //number of atoms read each pass
   int ini_vac();  //mean square displacement analysis
   int ana_vac();  //mean square displacement analysis
   int out_vac(ostream *);  //same method as in vac
  
   int doit(MODEL *);

};

#endif
