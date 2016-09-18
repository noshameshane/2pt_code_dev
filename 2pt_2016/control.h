#ifndef CONTROL_H
#define CONTROL_H
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include "constant.h"

using namespace std;

class CONTROL
{
public:
   CONTROL () {
      in_strt_flag=0;
      in_str_flag=0;
      in_ff_flag=0;
      in_trj_flag=0;
      in_grp_flag=0;
      out_grp_flag=0;
      out_trj_flag=0;
      out_bgf_flag=0;
      out_lmp_flag=0;
      out_img_flag=0;
      out_anim_flag=0;
      ana_flag=0;
      ana_stat_flag=0;
      ana_slab_flag=0;
      ana_rdf_flag=0;
      ana_msd_flag=0;
      ana_vac_flag=0;
      ana_mss_flag=0;
      ana_ptt_flag=0;
      ana_pto_flag=0;
      ana_ptv_flag=0;
      ana_tor_flag=0;
      ana_str_flag=0;
      ana_eps_flag=0;
      ana_prt_flag=0;
      ana_prdf_flag=0;
      ana_vod_flag=0;
      ana_view_flag=0;
      ana_iframe=1;
      ana_fframe=0;
      ana_sframe=1;
      ana_stat_type=0;
      ana_stat_blk=0;
      ana_slab_hbond=0; //hbond analysis
      ana_slab_direction=2; //z-direction
      ana_slab_thickness=0.1; //in angstroms
      ana_slab_nref=0;
      ana_slab_ref=NULL;
      strcpy(in_grp,"none");
      strcpy(in_cam,"none");
      in_trj=NULL;
      in_trj_n=0;
      ana_rdf_rcut=9;
      ana_rdf_rincrmnt=0.02;
      ana_rdf_grp[0]=ana_rdf_grp[1]=0;
      ana_prdf_rcut=9;
      ana_prdf_rincrmnt=0.02;
      ana_prdf_grp[0]=-69;
      ana_prdf_grp[1]=0;
      ana_msd_corlen=0.5;
      ana_msd_mem=500; //memory allocation in Mega Bytes
      ana_msd_cmc=0;   //no correction for center of mass motion
      ana_vac_corlen=0.5;
      ana_vac_mem=500; //memory allocation in Mega Bytes
      ana_vac_2pt=0; //2pt calculations 0:no, 1:yes, 2:2pt for molecules
      strcpy(ana_vac_const,"g");  //vac constraints
      strcpy(ana_vac_rotsym,"g"); //vac rotational symmetry
      strcpy(ana_vac_linear,"g"); //vac linear molecule
      strcpy(ana_vac_cnt,"\0"); //vac connectivity
      ana_tor_id=TORTYPE;
      ana_tor_binsize=30;
      ana_vod_blk=0;
      ana_vod_binsize=0;
      ana_str_direction=2;
      ana_str_avgblock=0;
      ana_eps_ewald=0.0001;
      ana_ptv_ang=(acos(0.98)*180.0/PI);
      ana_ptv_len=sqrt(80.0);
      ana_pto_delt=15;
      visu_flag=0;
      visu_iframe=0;
      visu_fframe=0;
      visu_sframe=0;
      visu_windx=640;
      visu_windy=480;
      visu_bkclr=0;
      strcpy(img_fmt,"gif");
      strcpy(out_img,"image");
      strcpy(out_anim,"anim");
      md_flag=0;       //flag for md simulation
      md_visu=1;       //flag for visualization
      md_runstep=100;  //number of steps to be run
      md_Tseed=123;    //random seed for initial temperature
      md_msgfreq=100;  //message output frequency
      md_dt=0.001;     //integration time step (in ps)
      md_Tset=-1;      //initial temperature
      md_nlskin=2;     //skin thickness for neighbor list (in A)

   }
   ~CONTROL () {
      if(in_trj!=NULL) delete [] in_trj;
      if(ana_slab_ref!=NULL) delete [] ana_slab_ref;
   }
   //flags
   int in_strt_flag;  //structure flag
   int in_str_flag;   //stress file flag
   int in_ff_flag;    //forcefield flag
   int in_trj_flag;   //input trj flag: c2trj or asctrj
   int in_grp_flag;   //group file flag
   int out_grp_flag;
   int out_trj_flag;
   int out_bgf_flag;
   int out_lmp_flag;
   int out_img_flag;
   int out_anim_flag;
   int ana_flag;       //analysis flag
   int ana_stat_flag;  //statistical analysis flag
   int ana_slab_flag;  //slab analysis flag
   int ana_rdf_flag;   //rdf analysis flag
   int ana_msd_flag;   //mean squre displacement analysis flag
   int ana_vac_flag;   //velocity autocorrelation analysis flag
   int ana_mss_flag;   //radius of gyration, moment of inertia, end2end dist analysis flag
   int ana_ptt_flag;   //polymer PTT order parameter calculation
   int ana_pto_flag;   //polymer PTT torsion calculation
   int ana_ptv_flag;   //polymer PTT vector (orientation) calculation
   int ana_tor_flag;   //torsion angle distribution analysis
   int ana_str_flag;   //Young's Modulus calculation
   int ana_eps_flag;   //dielectric constant
   int ana_prt_flag;   //proton analysis
   int ana_prdf_flag;   //proton rdf analysis
   int ana_vod_flag;   //void analysis
   int ana_view_flag;  //view model
   int visu_flag;  //graphics from opengl
   //values
   char in_strt[1024]; //structure file
   char in_stress[1024]; //stress file
   char in_ff[1024];   //forcefield file
   char (*in_trj)[1024];
   char in_cam[1024];  //camera parameter file
   char in_grp[1024];
   char out_grp[1024];
   char out_trj[1024];
   char out_bgf[1024];
   char out_lmp[1024];
   char out_img[1024];
   char img_fmt[4];   //image format
   char out_anim[1024];
   int in_trj_n;       //number of input trj file
   int ana_iframe;
   int ana_fframe;
   int ana_sframe;
   int ana_stat_type;
   int ana_stat_blk;
   int ana_slab_direction;
   double ana_slab_thickness; 
   int ana_slab_hbond; //hbond analysis
   int *ana_slab_ref;  //list of reference atoms
   int ana_slab_nref;  //number of reference atoms
   double ana_rdf_rcut;      //maximum r in rdf
   double ana_rdf_rincrmnt;  //increment in r 
   int ana_rdf_grp[2];
   double ana_msd_corlen;  //maximum correlation length (percentage of the total frame number) in MSD calc
   double ana_msd_mem; //memory allocation in Mega Bytes
   int    ana_msd_cmc; //correction for center of mass motion
   double ana_vac_corlen;  //maximum correlation length (percentage of the total frame number) in VAC calc
   double ana_vac_mem; //memory allocation in Mega Bytes
   int    ana_vac_2pt; //flag for 2pt calculations, 0:no, 1:yes, 2:2pt for molecules
   char   ana_vac_const[512]; //constraint degree of freedom
   char   ana_vac_rotsym[512]; //rotational symmetry number
   char   ana_vac_linear[512]; //flag for linear molecule
   char   ana_vac_cnt[1024]; //filename for connectivity
   int    ana_ptv[3];
   int    ana_tor_id;
   double ana_tor_binsize;
   double ana_prdf_rcut;      //maximum r in prdf
   double ana_prdf_rincrmnt;  //increment in r
   int ana_prdf_grp[2];       //grp1 ZEHO and grp2 id
   int    ana_vod_blk;      //size for block average
   double ana_vod_binsize;  //binsize for void analysis
   int ana_str_direction;
   int ana_str_avgblock;
   double ana_pto_delt;   //polymer PTT torsion calculation
   double ana_ptv_ang;    //angle used to define crystalline PTT
   double ana_ptv_len;    //length used to define crystalline PTT
   double ana_eps_ewald; //ewald accuracy in kcal/mol
   int visu_iframe;
   int visu_fframe;
   int visu_sframe;
   int visu_windx;
   int visu_windy;
   int visu_bkclr;  //background color
   int    md_flag;       //flag for md simulation
   int    md_visu;       //flag for visualization
   int    md_runstep;    //number of steps to be run
   int    md_Tseed;      //random seed for initial temperature
   int    md_msgfreq;    //message output frequency
   double md_dt;         //integration time step (in ps)
   double md_Tset;       //initial temperature
   double md_nlskin;     //skin thickness for neighbor list (in A)

   int rd_ctl(char*);
   int out_ctl(ostream *);
   int printsameline(char*,int,char*);
   int printsameline(char*,int,char*,int);
};

#endif
