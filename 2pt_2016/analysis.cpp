#include <iostream>
#include <fstream>
#include <cmath>
#include "constant.h"
#include "utility.h"
#include "analysis.h"

int ANALYSIS::doit(MODEL *in_model)
{
    model=in_model;
    if(model->ctl.ana_flag==0) return 0;

    if(model->ctl.in_trj_flag==filenotopen && model->ctl.ana_str_flag==0) return 0;

    int i,j;
    printf("#Analyzing Trj %s from frame %d to frame %d step %d\n",model->ctl.in_trj,model->ctl.ana_iframe,model->ctl.ana_fframe,model->ctl.ana_sframe);
    
    nframe=0;
    model->rd_grp(model->ctl.in_grp);

    if(model->ctl.ana_vac_flag)  ini_vac();
    for(i=model->ctl.ana_iframe;i<=model->ctl.ana_fframe;i+=model->ctl.ana_sframe) {
      model->trj.rd_frame(i);
      if(model->ctl.ana_vac_flag)  ana_vac();
      nframe++;
    }
    if(model->ctl.ana_vac_flag)  out_vac(&cout);

    return 0;

}

int ANALYSIS::ini_vac()
{
    int i,j,k;
    char null[1024];
    vacnsteps=(int)((model->ctl.ana_fframe-model->ctl.ana_iframe+1.0)/(model->ctl.ana_sframe));
    vacmaxf=(int)((model->ctl.ana_vac_corlen)*(vacnsteps-1))+1;
    int tot_N=vacnsteps+vacmaxf; //tot_N=nsteps+cor_length=vacnsteps+fraction*vacnsteps=vacnsteps+ana_vac_corlen*vacnsteps=vacnsteps+vacmaxf
    c2pt=model->ctl.ana_vac_2pt; 
    c2ptmf=0;

    switch(c2pt) {
      case 0: cout<<" Mode analysis will be performed"<<endl; c2pt=cmol=0; break;
      case 1: cout<<" Mode analysis with 2PT corrections will be performed"<<endl; c2pt=1; cmol=0; break;
      case 2: cout<<" Mode analysis with molecular corrections will be performed"<<endl; c2pt=0; cmol=1; break;
      case 3: cout<<" Mode analysis with molecular and 2PT corrections will be performed"<<endl; c2pt=cmol=1; break;
      case 4: cout<<" Mode analysis with molecular and 2PT-MF corrections will be performed"<<endl; c2pt=cmol=c2ptmf=1; break;
      case 5: cout<<" Mode analysis with molecular internal coordinates and 2PT corrections will be performed"<<endl; c2pt=cmol=1; break;
    }

    if(cmol==0) vactype=1; 
    else vactype=VELTYPE+1; 

    //specify if molecule in each group is linear
    if(model->ctl.ana_vac_linear[0]=='G'||model->ctl.ana_vac_linear[0]=='g') ;
    else {
      for(i=0;i<model->ngrp;i++) model->grp[i].linear = atoi(model->ctl.ana_vac_linear);
    }

    model->find_molingrp();

    //read connecitvity
    if(model->ctl.ana_vac_2pt==5) {
      for(i=0;i<model->nmol;i++) {
          model->mol[i].ini_iccnt();
          ifstream icnf(model->ctl.ana_vac_cnt,fstream::in);
          icnf.getline(null,1024);
          for(j=0;j<model->mol[i].natom;j++) {
              for(k=0;k<5;k++) icnf>>model->mol[i].iccnt[j][k];
              for(k=0;k<3;k++) model->mol[i].iccnt[j][k]--;
          }
      }
    }

    if(0) {
      for(i=0;i<model->nmol;i++)  {  
         for(j=0;j<model->mol[i].natom;j++)
           printf("m%d a%d b%d a%d d%d f%d s%d\n",i+1,j+1,model->mol[i].iccnt[j][0],model->mol[i].iccnt[j][1],model->mol[i].iccnt[j][2],model->mol[i].iccnt[j][3],model->mol[i].iccnt[j][4]);
      }
    }

    if(1) {
     for(i=0;i<model->ngrp;i++) printf("Group %d mass %f nmol %d\n",i+1,model->grp[i].mass,model->grp[i].nmol);
    }

    //memory usage statistics
    //Pvac:    ngrp*vactype*tot_N*sizeof(STAT)
    //vacT:    ngrp*vactype*sizeof(double)
    //vacDF:   ngrp*vactype*sizeof(double)
    //vacatom: natom*sizeof(ATOM)
    //pI:      nmol*3*sizeof(double)
    //in:      tot_N*sizeof(fftw_complex)
    //out:     tot_N*sizeof(fftw_complex)
    //pwr:     2*ngrp*tot_N/2*sizeof(double);
    double memusage= model->ngrp*vactype*( tot_N*sizeof(STAT) +sizeof(double)*2 )
                   + model->natom*sizeof(ATOM)
                   + model->nmol*3*sizeof(double)
                   + 2*tot_N*sizeof(fftw_complex)
                   + model->ngrp*( tot_N*sizeof(double) );
    memusage/=MEGABYTE;
    //vacread: vacatmperpass*sizeof(int)
    //vacvv:   vactype*vacnsteps*vacatmperpass*3*sizeof(float)

    //calc the memory needed to read in coordinates of an atom from the whole trj
    double memperatm=vactype*vacnsteps*sizeof(float)*3.0/MEGABYTE;
    vacatmperpass=(int)( (model->ctl.ana_vac_mem-memusage)/memperatm);
  
    if(vacatmperpass<1) {
       cout<<" Allocated memory ("<<model->ctl.ana_vac_mem<<" MB) insufficient for the calculation"<<endl;
       cout<<" Memory needed for reading one atom: "<<memperatm<<" MB"<<endl;
       cout<<" Memory needed for additional arrays: "<<memusage<<" MB"<<endl;
       cout<<" Minimum memory needed is "<<memusage+memperatm<<" MB"<<endl;
       cout<<" Please increase ANALYSIS_VAC_MEMORYMB"<<endl;
       exit(0);
    }
    cout<<" VAC extra memory needed: "<<memusage<<" MB"<<endl;
    cout<<"     Pvac[ngrp][vactype][tot_N]: "<<model->ngrp*vactype*tot_N*sizeof(STAT)/MEGABYTE<<endl;
    cout<<"     vacT[ngrp][vactype]       : "<<model->ngrp*vactype*sizeof(double)/MEGABYTE<<endl;
    cout<<"     vacDF[ngrp][vactype]      : "<<model->ngrp*vactype*sizeof(double)/MEGABYTE<<endl;
    cout<<"     vacatom[natom]            : "<<model->natom*sizeof(ATOM)/MEGABYTE<<endl;
    cout<<"     pI[nmol][3]               : "<<model->nmol*3*sizeof(double)/MEGABYTE<<endl;
    cout<<"     fftw in out               : "<<2*tot_N*sizeof(fftw_complex)/MEGABYTE<<endl;
    cout<<"     pwr[2][ngrp][tot_N/2]     : "<<model->ngrp*tot_N*sizeof(double)/MEGABYTE<<endl;
    cout<<"     ngroup "<<model->ngrp<<" vactype "<<vactype<<" tot_N "<<tot_N<<" STAT "<<sizeof(STAT)/MEGABYTE<<endl;
    cout<<" VAC memory needed for one atom: "<<memperatm<<" MB"<<endl;

    vacnpass=(int)(model->natom/vacatmperpass+1.0);
    vacatmperpass=(vacatmperpass>model->natom?model->natom:vacatmperpass);
    vacnpass=(int)(model->natom/vacatmperpass)+(model->natom%vacatmperpass==0?0:1);
    printf(" ALLOC MEM %5.1f MB, Max Atom per Pass %d, nPass %d\n",model->ctl.ana_vac_mem,vacatmperpass,vacnpass);

    for(i=0;i<model->nmol;i++) {
        if(model->mol[i].natom>vacatmperpass) {
           printf(" Error: Max allowed atoms per Pass (%d) is less than the number of atoms (%d) contained in molecule %d\n",vacatmperpass,model->mol[i].natom,i+1);
           printf(" Please increase ANALYSIS_VAC_MEMORYMB\n");
           exit(0);
        }
    }

    Pvac=new STAT **[model->ngrp];
    for(i=0;i<model->ngrp;i++) {
        Pvac[i]=new STAT *[vactype];
        for(j=0;j<vactype;j++) Pvac[i][j]=new STAT [tot_N];
    }

    vacatom=new ATOM [model->natom];
    for(i=0;i<model->natom;i++) {
        vacatom[i].mass=model->atom[i].mass;
        vacatom[i].mol=model->atom[i].mol;
        vacatom[i].id=model->atom[i].id;
    }

    tatom=model->atom;
    model->atom=vacatom;
    for(i=0;i<model->trj.ntrjf;i++) {
        model->trj.strj[i].atom=model->atom;
        model->trj.strj[i].cell=&vaccell;
        model->trj.strj[i].prp =&vacprp;
    }

    model->trj.rd_frame(model->ctl.ana_iframe);
    vacdtime=vacprp.time;
    model->trj.rd_frame(model->ctl.ana_iframe+model->ctl.ana_sframe);
    vacdtime=vacprp.time-vacdtime;

    model->atom=tatom;
    for(i=0;i<model->trj.ntrjf;i++) {
        model->trj.strj[i].atom=model->atom;
        model->trj.strj[i].cell=&(model->cell);
        model->trj.strj[i].prp =&(model->prp);
    }


    vacread=new int [vacatmperpass]; //id of atoms read in each pass

    vacvv=new float ***[vactype];
    for(k=0;k<vactype;k++) {
        vacvv[k]=new float **[vacnsteps];
        for(i=0;i<vacnsteps;i++) {
           vacvv[k][i]=new float *[vacatmperpass];
           for(j=0;j<vacatmperpass;j++) vacvv[k][i][j]=new float [3];
        }
    }
    vacipass=0;
    vactatmpassed=0;  //total atoms passed
    vacatmpass=0;     //# atom in current pass
    vacmolpass=0;     //# molecule in current pass
    vactmolpassed=0;  //total molecules passed
    vaclastmol=0;     //last molecule id from previous pass

    if(cmol==0) {
       if(vacatmperpass>model->natom) j=model->natom;
       else j=vacatmperpass;
       for(i=0;i<j;i++) vacread[vacatmpass++]=i;
    } else {
       int tmpi=0;
       for(i=vaclastmol;i<model->nmol;i++) {
          tmpi += model->mol[i].natom;
          if(tmpi<=vacatmperpass) {
            for(j=0;j<model->mol[i].natom;j++) vacread[vacatmpass++]=model->mol[i].atm[j];
            vacmolpass++;
          } else {
            break; //break i-loop
          }
       }
       vactmolpassed+=vacmolpass;
       vaclastmol=i;
    }
    vactatmpassed +=vacatmpass;

    vacE=0;
    vacV=0;
    vacP=0;
    vacT=new double *[model->ngrp];  for(i=0;i<model->ngrp;i++) vacT[i]=new double [vactype];
    vacDF=new double *[model->ngrp]; for(i=0;i<model->ngrp;i++) vacDF[i]=new double [vactype];
    double trdf; //translation and rotational degrees of freedom removed
    if(model->ctl.in_trj_flag==usrtrj) model->periodic=1;
    if (!model->periodic){
          if(model->natom==1) trdf=3;
          else if(model->natom==2) trdf=5;
          else trdf=6;
          /* trdf=3; this setup is for lammps, which does not remove rotational degrees of freedom*/
    } else {
       trdf=3;
    }
    trjT=0;

    //determine the degree of freedom in each group
    for(i=0;i<model->ngrp;i++) {
        for(j=0;j<vactype;j++) vacT[i][j]=0;
        if(cmol==0) vacDF[i][0]=3.0*model->grp[i].natom;
        else {
          vacDF[i][vtrans]=0;
          vacDF[i][vrotat]=0;
          vacDF[i][vangul]=0;
          vacDF[i][vimvib]=0;
          vacDF[i][vtotal]=0;
          for(j=0;j<model->grp[i].nmol;j++) {
              int tmpi=model->mol[ model->grp[i].mol[j] ].natom; /*number of atoms in this molecule*/
              int tmpj=0; 
              if(model->grp[i].linear>0) model->mol[ model->grp[i].mol[j] ].linear=tmpj=1;
              switch(tmpi) {//check number of atoms in each molecule
                 case 1: //monoatomic
                      vacDF[i][vtrans]+=3.0;
                      vacDF[i][vrotat]+=0.0;
                      vacDF[i][vangul]+=0.0;
                      vacDF[i][vimvib]+=0.0;
                      vacDF[i][vtotal]+=3.0;
                      break;
                 case 2: //diatomic
                      vacDF[i][vtrans]+=3.0;
                      vacDF[i][vrotat]+=2.0;
                      vacDF[i][vangul]+=2.0;
                      vacDF[i][vimvib]+=1.0;
                      vacDF[i][vtotal]+=6.0;
                      break;
                 default: //polyatomic
                      vacDF[i][vtrans]+=3.0;
                      vacDF[i][vrotat]+=(3.0-tmpj); 
                      vacDF[i][vangul]+=(3.0-tmpj); 
                      vacDF[i][vimvib]+=(3.0*tmpi-6.0+tmpj);
                      vacDF[i][vtotal]+=(3.0*tmpi);
                      break;
              }
           }
        }
    }
    //now correct the degrees of freedom for stretching constraints
    if(model->ctl.ana_vac_const[0]=='G'||model->ctl.ana_vac_const[0]=='g') ;
    else model->grp[model->ngrp-1].constraint= atoi(model->ctl.ana_vac_const);
    for(i=0;i<model->ngrp;i++) {
       if(cmol==0) {
         vacDF[i][0]-= model->grp[i].constraint;
       } else {
         vacDF[i][vimvib]-= model->grp[i].constraint;
         vacDF[i][vtotal]-= model->grp[i].constraint;
       }
    }

    //now specify the rotational symmetry for each group
    if(model->ctl.ana_vac_rotsym[0]=='G'||model->ctl.ana_vac_rotsym[0]=='g') ;
    else {
      for(i=0;i<model->ngrp;i++) model->grp[i].rotsym = atoi(model->ctl.ana_vac_rotsym);
    }

    //trdf will be evenly distributed to each degrees of freedom
    for(i=0;i<model->ngrp;i++) {
        for(j=0;j<vactype;j++) { 
            vacDF[i][j]-= trdf*vacDF[i][j]/vacDF[model->ngrp-1][vactype-1];
            //cout<<"DF group "<<i<<" type "<<j<<" "<<vacDF[i][j]<<endl; 
        }
        model->grp[i].eng=0; //clear group energy
    }

    if(cmol) {
       pI=new double [model->nmol][3];
       for(i=0;i<model->nmol;i++) {
          for(k=0;k<3;k++) pI[i][k]+=model->mol[i].pI[k];
       }
    }

    if(model->ctl.ana_vac_2pt==5) {
      PvacInt=new STAT **[model->ngrp];
      vacTInt=new double *[model->ngrp];
      vacDFInt=new double *[model->ngrp];
      for(i=0;i<model->ngrp;i++) {
          int tmp=3*model->mol[(model->grp[i].mol[0])].natom-6;
          vacTInt[i]=new double [tmp];
          vacDFInt[i]=new double [tmp];
          Pvac[i]=new STAT *[tmp];
          for(j=0;j<tmp;j++) Pvac[i][j]=new STAT [tot_N];
      }
 
    }
    cout<<flush;
    return 0;
}

int ANALYSIS::ana_vac()
{
    int i,k,id;

    if(0) {
      for(i=0;i<model->natom;i++) {
          printf("step %d atom %d v %8.6f %8.6f %8.6f c %8.6f %8.6f %8.6f\n",nframe,id,model->atom[id].vel[0],model->atom[id].vel[1],model->atom[id].vel[2],model->atom[id].pv[0],model->atom[id].pv[1],model->atom[id].pv[2]);
      }
    }

    if(cmol==0) {
       for(i=0;i<vacatmpass;i++) {
          id=vacread[i];
          for(k=0;k<3;k++) vacvv[0][nframe][i][k]=model->atom[id].vel[k];
       }
    } else {
       model->cal_vcomp(); //perform velocity decomposition: vtrn, vrot, vvib
       for(i=0;i<vacatmpass;i++) {
          id=vacread[i];
          for(k=0;k<3;k++) {
              vacvv[vtrans][nframe][i][k]=model->atom[id].vt[k]; 
              vacvv[vrotat][nframe][i][k]=model->atom[id].vr[k]; 
              vacvv[vimvib][nframe][i][k]=model->atom[id].vv[k]; 
              vacvv[vtotal][nframe][i][k]=model->atom[id].vel[k]; 
          }
       } 
       id=vaclastmol-vacmolpass; //first molecular id
       for(i=id;i<vaclastmol;i++) { //i is the molecular id
          for(k=0;k<3;k++) 
              vacvv[vangul][nframe][i-id][k]=model->mol[i].anguv[k];
       }
       for(i=0;i<model->nmol;i++) {
          for(k=0;k<3;k++) pI[i][k]+=model->mol[i].pI[k];
       }
    }

    if(0) {
        for(i=0;i<vacatmpass;i++) {
           id=vacread[i];
           //printf("step %d atom %d vtrans %8.6f %8.6f %8.6f\n",nframe,id,vacvv[vtrans][nframe][i][0],vacvv[vtrans][nframe][i][1],vacvv[vtrans][nframe][i][2]);
           //printf("step %d atom %d vrotat %8.6f %8.6f %8.6f\n",nframe,id,vacvv[vrotat][nframe][i][0],vacvv[vrotat][nframe][i][1],vacvv[vrotat][nframe][i][2]);
           //printf("step %d atom %d vimvib %8.6f %8.6f %8.6f\n",nframe,id,vacvv[vimvib][nframe][i][0],vacvv[vimvib][nframe][i][1],vacvv[vimvib][nframe][i][2]);
           //printf("step %d atom %d vtotal %8.6f %8.6f %8.6f\n",nframe,id,vacvv[vtotal][nframe][i][0],vacvv[vtotal][nframe][i][1],vacvv[vtotal][nframe][i][2]);
           printf("step %d atom %d mass %f vtotal %8.6f %8.6f %8.6f\n",nframe,id,model->atom[id].mass,vacvv[vtotal][nframe][i][0],vacvv[vtotal][nframe][i][1],vacvv[vtotal][nframe][i][2]);
        }
    }

    int g1,g2;
    for(i=0;i<model->ngrp;i++) {
        for(k=0;k<model->grp[i].natom;k++) {
            id= model->grp[i].atm[k];
            model->grp[i].eng += model->atom[id].eng;
        }
    }

    trjT+=model->prp.T;
    vacE+=model->prp.Et;
    vacV+=model->prp.V;
    vacP+=model->prp.P;

    return 0;
}

int ANALYSIS::out_vac(ostream *outf) //mass weighted vac (vac_driver approach)
{ 
    int i,j,id,k,l,g1,g2,tp;
    char null[1024]; 
    double dst;
    STAT *atmvac;

    ofstream outvac("vac.vac",ios::out);
    ofstream outpwr("vac.pwr",ios::out);
    ofstream out3n("vac.3n",ios::out);
    ofstream outthermo("vac.thermo",ios::out);

    tatom=model->atom;
    model->atom=vacatom;
    if(cmol) {
       for(i=0;i<model->nmol;i++) {
           for(j=0;j<model->mol[i].natom;j++) model->mol[i].atom[j] = & vacatom[ model->mol[i].atm[j] ];
       }
    }
    for(i=0;i<model->trj.ntrjf;i++) {
        model->trj.strj[i].atom=model->atom;
        model->trj.strj[i].cell=&vaccell;
        model->trj.strj[i].prp =&vacprp;
    }

    trjT/=nframe; //average temperature
    vacE*=(caltoj/nframe); //average total energy
    vacV/=nframe; //average volume
    vacP/=nframe; //average pressure

    for(i=0;i<model->ngrp;i++) {
        if(model->grp[i].eng==0) model->grp[i].eng=vacE*vacDF[i][vactype-1]/vacDF[model->ngrp-1][vactype-1]; //group energy estimated from DF
        else model->grp[i].eng *= (caltoj/nframe); //average group energy in kJ/mol
    }

    //check MD energy consistency
    if(fabs((vacE-model->grp[ model->ngrp-1 ].eng)/vacE)>0.001) {
       cout<<" Warning: MD energy ("<<vacE<<" kJ/mol) differs from sum of atomic energy ("<<model->grp[ model->ngrp-1 ].eng<<" kJ/mol) by "<<vacE-model->grp[ model->ngrp-1 ].eng<<" kJ/mol"<<endl;
    }

    if(cmol) {
       for(i=0;i<model->nmol;i++) {
           for(k=0;k<3;k++) pI[i][k]/=nframe; //average principle moment of inertia
       }
       if(0) {
        for(i=0;i<model->nmol;i++) printf("%.2e %.2e %.2e\n",pI[i][0],pI[i][1],pI[i][2]);
       }
    }
    
    int tot_N=vacnsteps+vacmaxf; //tot_N=nsteps+cor_length=vacnsteps+fraction*vacnsteps=vacnsteps+ana_vac_corlen*vacnsteps=vacnsteps+vacmaxf
    int nused=tot_N/2;
    fftw_complex *in,*out;
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * tot_N);
    out= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * tot_N);
    double *power_spectrum;
    fftw_plan p1;
    fftw_plan p2;
    p1 = fftw_plan_dft_1d(tot_N,in,out,FFTW_BACKWARD,FFTW_ESTIMATE); //backward(+1)
    p2 = fftw_plan_dft_1d(tot_N,in,out,FFTW_FORWARD,FFTW_ESTIMATE);  //forward(-1)

    while(1) {
       vacipass++;
       printf(" Pass %d of %d: atoms %d/%d ",vacipass,vacnpass,vactatmpassed,model->natom);
       if(cmol) printf("mols %d/%d ",vactmolpassed,model->nmol);
       cout<<endl;
       //determine VAC
       for(tp=0;tp<vactype;tp++) {
          if(tp!=vangul) { //translation, rotation, vibration, total
            for(i=0;i<vacatmpass;i++) {
               atmvac=new STAT [tot_N];
               id=vacread[i];
               for(k=0;k<3;k++) {
                  for(j=0;j<vacnsteps;j++) {
                     in[j][0]= vacvv[tp][j][i][k];  //velocity in k direction of atom i in the j_th step
                     in[j][1]= 0;
                  }
                  for(;j<tot_N;j++) in[j][0]=in[j][1]=0;
                  fftw_execute(p1); // do backward(+1) fft
     
                  for(j=0;j<tot_N;j++) {
                      in[j][0]= (out[j][0]*out[j][0] + out[j][1]*out[j][1])/tot_N; //work(j)=work(j)*dconjg(work(j))/tot_N
                      in[j][1]= 0;
                  }
                  fftw_execute(p2); // do forward(-1) fft
     
                  for(j=0;j<tot_N;j++) {
                     atmvac[j].add_data(model->atom[id].mass*out[j][0]/vacnsteps);
                  }
               }
               
               for(g1=0;g1<model->ngrp;g1++) {
                  for(g2=0;g2<model->grp[g1].natom;g2++) {
                     if(id==model->grp[g1].atm[g2]) {
                        for(j=0;j<tot_N;j++)  Pvac[g1][tp][j].add_data(atmvac[j].sum);
                     }
                  }
               }
               delete [] atmvac;
            }
          } else { //angular velocity
            id=vaclastmol-vacmolpass; //first molecular id
            for(i=id;i<vaclastmol;i++) {  //i is the molecular id
               atmvac=new STAT [tot_N];
               
               for(k=0;k<3;k++) {
                  for(j=0;j<vacnsteps;j++) {
                     in[j][0]= vacvv[tp][j][i-id][k];  //velocity in k direction of molecule i in the j_th step
                     in[j][1]= 0;
                  }
                  for(;j<tot_N;j++) in[j][0]=in[j][1]=0;
                  fftw_execute(p1); // do backward(+1) fft
     
                  for(j=0;j<tot_N;j++) {
                      in[j][0]= (out[j][0]*out[j][0] + out[j][1]*out[j][1])/tot_N; //work(j)=work(j)*dconjg(work(j))/tot_N
                      in[j][1]= 0;
                  }
                  fftw_execute(p2); // do forward(-1) fft
     
                  for(j=0;j<tot_N;j++) {
                     atmvac[j].add_data(1.0*out[j][0]/vacnsteps); //use unity for mass, angular vel has been weighted by sqrt(I).
                  }
               }
               
               for(g1=0;g1<model->ngrp;g1++) {
                  for(g2=0;g2<model->grp[g1].nmol;g2++) {
                     if(i==model->grp[g1].mol[g2]) {
                        for(j=0;j<tot_N;j++)  Pvac[g1][tp][j].add_data(atmvac[j].sum);
                     }
                  }
               }
               delete [] atmvac;
            }
          }
       }

       if(vactatmpassed<model->natom) {  //read the trj file again if not all atom data were loaded
          vacatmpass =0;
          vacmolpass =0;
          //determine which atoms to read
          if(cmol==0) {
             j= model->natom - vactatmpassed;
             if(j>vacatmperpass) j=vacatmperpass;
             for(i=0;i<j;i++) vacread[vacatmpass++]=vactatmpassed+i;
          } else {
             int tmpi=0;
             for(i=vaclastmol;i<model->nmol;i++) {
                tmpi += model->mol[i].natom;
                if(tmpi<=vacatmperpass) {
                  for(j=0;j<model->mol[i].natom;j++) vacread[vacatmpass++]=model->mol[i].atm[j];
                  vacmolpass++;
                } else {
                  break; //break i-loop
                }
             }
             vactmolpassed+=vacmolpass;
             vaclastmol=i;
          }

          //read the trj file again
          for(i=0;i<vactype;i++) for(j=0;j<vacnsteps;j++) for(k=0;k<vacatmperpass;k++) for(l=0;l<3;l++) vacvv[i][j][k][l]=0;
          l=0;
          for(j=model->ctl.ana_iframe;j<=model->ctl.ana_fframe;j+=model->ctl.ana_sframe) {
              model->trj.rd_frame(j);
              if(cmol==0) {
                 for(i=0;i<vacatmpass;i++) {
                    id=vacread[i];
                    for(k=0;k<3;k++) vacvv[0][l][i][k]=model->atom[id].vel[k];
                 }
              } else {
                 model->cal_vcomp();
                 for(i=0;i<vacatmpass;i++) {
                    id=vacread[i];
                    for(k=0;k<3;k++) {
                        vacvv[vtrans][l][i][k]=model->atom[id].vt[k]; 
                        vacvv[vrotat][l][i][k]=model->atom[id].vr[k]; 
                        vacvv[vimvib][l][i][k]=model->atom[id].vv[k]; 
                        vacvv[vtotal][l][i][k]=model->atom[id].vel[k]; 
                    }
                 }
                 id=vaclastmol-vacmolpass; //first molecular id
                 for(i=id;i<vaclastmol;i++) {
                    for(k=0;k<3;k++) 
                        vacvv[vangul][l][i-id][k]=model->mol[i].anguv[k];
                 }
              }
              l++;
          }
          vactatmpassed +=vacatmpass;
       } else break;
    }
    
    sprintf(null," iframe %d, fframe %d, step %d, frames used %d, max correlation length %f ps, time interval %f ps",model->ctl.ana_iframe,model->ctl.ana_fframe,model->ctl.ana_sframe,nframe,vacmaxf*vacdtime,vacdtime);
    *outf<<null<<endl;
    sprintf(null," Velocity autocorrelation function (g/mol*A^2/ps^2) analysis");
    *outf<<null<<endl;
    outvac<<null<<endl;
    sprintf(null," %13s","time(ps)"); *outf<<null; outvac<<null;
    for(j=0;j<model->ngrp;j++) {  
        if(c2pt==0&&cmol==0) { sprintf(null,"  %6s[G%03d]","VAC",j+1); *outf<<null; outvac<<null;}
        else {
          if(cmol) sprintf(null,"  %6s[G%03d]  %6s[G%03d]  %6s[G%03d]  %6s[G%03d]","VACcmt",j+1,"VACcmo",j+1,"VACimv",j+1,"VACcmr",j+1); *outf<<null; outvac<<null;
          sprintf(null,"  %6s[G%03d]","VACtot",j+1); *outf<<null; outvac<<null;
        }
    }
    *outf<<endl; outvac<<endl;
    for(i=0;i<vacmaxf;i++) { 
       sprintf(null," %13.3f",vacdtime*i); *outf<<null; outvac<<null;
       for(j=0;j<model->ngrp;j++) {
          for(k=0;k<vactype;k++) {
             Pvac[j][k][i].cal_avg();
             sprintf(null," %13.3f",Pvac[j][k][i].sum); *outf<<null; outvac<<null;
             //Pvac[j][k][i].sum is the sum of vac of all atoms 
             //Pvac[j][k][i].npt is the number of atoms
          }
       }
       //sprintf(null," %d",Pvac[0][i].npt); *outf<<null;
       *outf<<endl; outvac<<endl;
    }

    //determine the real temperature of each group
    //real temp from df*kT/2=sum(mv^2/2), sum(mv^2)=Pvac[j][0]
    for(j=0;j<model->ngrp;j++) {
        for(i=0;i<vactype;i++) {
           if((int)(vacDF[j][i])==0)  vacT[j][i] = 0;
           else vacT[j][i] = Pvac[j][i][0].sum/(0.1*R*vacDF[j][i]);
        }
    }
    i=model->ngrp-1; j=vactype-1;
    if(fabs(vacT[i][j]-trjT)/trjT > 0.0001) {
      cout<<" Trajectory temperature "<<trjT<<" and velocity temperature "<<vacT[i][j]<<" differ by "<<(vacT[i][j]-trjT)/trjT*100.0<<"%"<<endl;
      cout<<" Check the degrees of freedom calculation in the vac code"<<endl;
    } else cout<<" Temperature consistent"<<endl;
  
    //determine the power spectrum
    double pwrfreq=1.0E10/(vacdtime*tot_N*vlight);  //in cm-1
    for(j=0;j<model->ngrp;j++) {
       for(k=0;k<vactype;k++) {
          if(vacT[j][k]<=0) for(i=0;i<tot_N;i++) in[i][0]=in[i][1]=0;
          else for(i=0;i<tot_N;i++)  {in[i][0]=Pvac[j][k][i].sum*20.0/(R*vacT[j][k]); in[i][1]=0; } // in[] is dimensionless 
          fftw_execute(p1);
          for(i=0;i<nused;i++) {
             Pvac[j][k][i].max=out[i][0]*vacdtime*vlight*1.0E-10; //pwr(i)=dble(work(i))*time_step*cspeed, out[] dimensionless, Pvac[j][k][i].max is in unit of cm
          }
          Pvac[j][k][0].min = Pvac[j][k][0].max*pwrfreq*0.5;  //use min to store PWR integration
          for(i=1;i<nused-1;i++) Pvac[j][k][i].min = Pvac[j][k][i-1].min + Pvac[j][k][i].max*pwrfreq;
          Pvac[j][k][i].min = Pvac[j][k][i-1].min + Pvac[j][k][i].max*pwrfreq*0.5;
       }
    }
    fftw_destroy_plan(p1);
    fftw_destroy_plan(p2);
    fftw_free(in); 
    fftw_free(out);

    k=model->ngrp;
    double f2pt[2][k],K2pt[2][k],fmf[2][k];
    double pwr[2][k][nused];
    double tmpg,tmps,tmp;
    double wep[k],wsp[k],wspd[k],wap[k],wcvp[k],wsehs[k],wsr[k],war[k],y[2][k],hsdf[2][k],ttdf[2][k];

    if(vacV<=0) {
       printf("\nstlin original volume %lf \n",vacV);
       vacV= 9.9668*model->natom;
       printf("\nstlin estimate volume to be %lf from water 9.9668 A3/atom\n",vacV);
    }

    //calculation of partial molar volume using volume ratio in the grp file
    double vacPMV[k],vscale=0;
    for(j=0;j<k-1;j++) vscale+= (model->grp[j].nmol*model->grp[j].pmvr);
    vscale=vacV/vscale; //scaling factor for partial molar volume

    for(j=0;j<k-1;j++) vacPMV[j]= model->grp[j].nmol*model->grp[j].pmvr*vscale;
    vacPMV[k-1]=vacV; //total volume of system;

    if(c2pt==1) {
       for(j=0;j<model->ngrp;j++) {
          int ngmol= model->grp[j].nmol; //number of molecules in the group of interest
          int ntmol= model->nmol; //number of molecules in the whole system
          double gmass= model->grp[j].mass; //total mass of atoms (molecules) in the group of interest
          double rT[3];  //rotation temperatrues
          double rs=model->grp[j].rotsym*1.0; //rotational symmetry number
          vacV=vacPMV[j]; //now vacV is the partial box volume for each group
          rT[0]=rT[1]=rT[2]=0.0;
          printf("\nGroup %d nmol %d mass %f",j+1,ngmol,gmass);

          if(cmol) {
            for(i=0;i<model->grp[j].nmol;i++) {
               id=model->grp[j].mol[i];
               tmp=(h*h/(8.0*PI*PI*pI[id][0]*kb)); if(tmp<0) tmp=0; rT[0]+=tmp;
               tmp=(h*h/(8.0*PI*PI*pI[id][1]*kb)); if(tmp<0) tmp=0; rT[1]+=tmp;
               tmp=(h*h/(8.0*PI*PI*pI[id][2]*kb)); if(tmp<0) tmp=0; rT[2]+=tmp;
               if(tmp<0) tmp=0;
               //printf("\nmol %d I %e %e %e Trot3 %e\n",id+1,pI[id][0],pI[id][1],pI[id][2],tmp); 
            }
            rT[0]/=model->grp[j].nmol;  
            rT[1]/=model->grp[j].nmol;  
            rT[2]/=model->grp[j].nmol;  
            if(abs(model->grp[j].linear)) rT[2]=-999;
            printf("\n rotational temperatures %e %e %e K\n",rT[0],rT[1],rT[2]);
          }

          for(i=0;i<nused;i++) pwr[0][j][i]=Pvac[j][vtrans][i].max; //translation, in cm
          if(cmol) for(i=0;i<nused;i++) pwr[1][j][i]=Pvac[j][vangul][i].max; //rotation, in cm
          K2pt[0][j]=(pwr[0][j][0]/vlight*1E-2)/ngmol*sqrt(PI*Na*kb*vacT[j][vtrans]/(gmass*1e-3/ngmol))*2.0/9.0*pow(ngmol/vacV,1.0/3.0)*1E10*pow(6.0/PI,2.0/3.0); //translation
          printf("\nstlin trans s0=%lf T=%lf mass=%lf V=%lf ngmol=%d ntmol=%d K=%lf",pwr[0][j][0],vacT[j][vtrans],gmass/ngmol,vacV,ngmol,ntmol,K2pt[0][j]);
          f2pt[0][j]=search2PT(K2pt[0][j]);
          fmf[0][j]=f2pt[0][j];
          if(c2ptmf) twoPTmf(pwr[0][j],pwr[0][j][0],nused,pwrfreq,ngmol,f2pt[0][j],fmf[0][j]); //determine the value of fmf

          double ry[2],hs_sigma[2];
          ry[0]=pow(f2pt[0][j]/K2pt[0][j],1.5);

          if (f2pt[0][j]==fmf[0][j]) ry[0]=pow(f2pt[0][j]/K2pt[0][j],1.5);  
          else ry[0]=pow(fmf[0][j]/K2pt[0][j],1.5);  //stlin check why change to fmf[0][j]?

          hs_sigma[0]=pow(ry[0]*6.0/PI/ngmol*vacV,1.0/3.0);
          hsdf[0][j]=HSDF(pwr[0][j],nused,pwrfreq,ngmol,f2pt[0][j],fmf[0][j]); /*determine hard sphere degrees of freedom hsdf*/
          ttdf[0][j]=TTDF(pwr[0][j],nused,pwrfreq); /*determine total degrees of freedom ttdf*/
          y[0][j]=ry[0]*(hsdf[0][j]/ttdf[0][j]);
          if(rT[0]>0) {
            K2pt[1][j]=(pwr[1][j][0]/vlight*1E-2)/ngmol*sqrt(PI*Na*kb*vacT[j][vangul]/(gmass*1e-3/ngmol))*2.0/9.0*pow(ngmol/vacV,1.0/3.0)*1E10*pow(6.0/PI,2.0/3.0); //rotation
            //K2pt[1][j]=K2pt[0][j]; //tie f_trans=f_rot
            f2pt[1][j]=search2PT(K2pt[1][j]);
            fmf[1][j]=f2pt[1][j];
            if(c2ptmf) twoPTmf(pwr[1][j],pwr[1][j][0],nused,pwrfreq,ngmol,f2pt[1][j],fmf[1][j]);
            if (f2pt[1][j]==fmf[1][j]) ry[1]=pow(f2pt[1][j]/K2pt[1][j],1.5);
            else ry[1]=pow(fmf[1][j]/K2pt[1][j],1.5);  //stlin check why change to fmf[0][j]?
            hs_sigma[1]=pow(ry[1]*6.0/PI/ngmol*vacV,1.0/3.0);
            hsdf[1][j]=HSDF(pwr[1][j],nused,pwrfreq,ngmol,f2pt[1][j],fmf[1][j]); /*determine hard sphere degrees of freedom hsdf*/
            ttdf[1][j]=TTDF(pwr[1][j],nused,pwrfreq); /*determine total degrees of freedom ttdf*/
            y[1][j]=ry[1]*(hsdf[1][j]/ttdf[1][j]);
          } else K2pt[1][j]=f2pt[1][j]=ry[1]=hs_sigma[1]=hsdf[1][j]=ttdf[1][j]=y[1][j]=0;

          printf("\nstlin rotat s0=%lf T=%lf mass=%lf V=%lf ngmol=%d ntmol=%d K=%lf",pwr[1][j][0],vacT[j][vangul],gmass/ngmol,vacV,ngmol,ntmol,K2pt[1][j]);
          printf("\nstlin trans fhs fhc packfrac hsdf tdf  %lf %lf %lf %lf %lf",f2pt[0][j],fmf[0][j],hsdf[0][j]/ttdf[0][j],hsdf[0][j],ttdf[0][j]);
          printf("\nstlin rotat frr fhc packfrac hsdf tdf  %lf %lf %lf %lf %lf",f2pt[1][j],fmf[1][j],hsdf[1][j]/ttdf[1][j],hsdf[1][j],ttdf[1][j]);
          if(y[0][j]>0.74) f2pt[0][j]=f2pt[1][j]=0; /*packing fraction too large, ignore partition*/
          HSweighting(&wep[j],&wap[j],&wcvp[j],&wsehs[j],&wsp[j],&wspd[j],y[0][j],gmass/ngmol,hsdf[0][j]/3.0,vacT[j][vtrans],vacT[j][vangul],vacV,&wsr[j],&war[j],rT,rs);
   
          printf("\nstlin constant K                     %lf %lf",K2pt[0][j],K2pt[1][j]);
          printf("\nstlin fludicity                      %lf %lf",f2pt[0][j],f2pt[1][j]);
          printf("\nstlin fludicity mf                   %lf %lf",fmf[0][j],fmf[1][j]);
          printf("\nstlin refernce Hard Sphere packing   %lf %lf",ry[0],ry[1]);
          printf("\nstlin Hard Sphere diameter           %lf %lf",hs_sigma[0],hs_sigma[1]);
          printf("\nstlin Hard Sphere packing fraction   %lf %lf",y[0][j],y[1][j]);
          printf("\nstlin Effective N for hard sphere    %lf %lf",ngmol*(hsdf[0][j]/(3.0*ngmol-3.0)),ngmol*(hsdf[1][j]/(3.0*ngmol-3.0)));
          printf("\nstlin Hard Sphere reduced density    %lf %lf",y[0][j]*6.0/PI,y[1][j]*6.0/PI);
          printf("\nstlin True Hard Sphere dof           %lf %lf",hsdf[0][j],hsdf[1][j]);
          printf("\nstlin temperature rotT3              %lf %lf %lf %lf",vacT[j][vangul],rT[0],rT[1],rT[2]);
          printf("\nstlin Weighting Stran Srot           %lf %lf\n",wsp[j],wsr[j]);

          //Perform 2PT decomposition

          for(i=0;i<nused;i++) {
              twoPT(&tmpg,&tmps,pwr[0][j][0],pwr[0][j][i],pwrfreq*i,ngmol,f2pt[0][j],fmf[0][j]);
              Pvac[j][vtrans][i].a=tmpg; //use a to store hard sphere gas contribution
              Pvac[j][vtrans][i].b=tmps; //use b to store solid contribution
              if(rT[0]>0) twoPT(&tmpg,&tmps,pwr[1][j][0],pwr[1][j][i],pwrfreq*i,ngmol,f2pt[1][j],fmf[1][j]); else tmpg=tmps=0;
              if(cmol) { 
                   Pvac[j][vangul][i].a=tmpg; //use a to store hard sphere gas contribution
                   Pvac[j][vangul][i].b=tmps; //use b to store solid contribution
              }
          }

          //Perform DoS integration
          Pvac[j][vtrans][0].SEa = Pvac[j][vtrans][0].a*pwrfreq*0.5;  //use SEa to store PWR integration of hs 
          Pvac[j][vtrans][0].SEb = Pvac[j][vtrans][0].b*pwrfreq*0.5;  //use SEb to store PWR integration of solid 
          if(cmol) {
          Pvac[j][vangul][0].SEa = Pvac[j][vangul][0].a*pwrfreq*0.5;  //use SEa to store PWR integration of hs 
          Pvac[j][vangul][0].SEb = Pvac[j][vangul][0].b*pwrfreq*0.5;  //use SEb to store PWR integration of solid 
          }

          for(i=1;i<nused-1;i++) { 
                Pvac[j][vtrans][i].SEa = Pvac[j][vtrans][i-1].SEa + Pvac[j][vtrans][i].a*pwrfreq;
                Pvac[j][vtrans][i].SEb = Pvac[j][vtrans][i-1].SEb + Pvac[j][vtrans][i].b*pwrfreq;
                if(cmol) {
                Pvac[j][vangul][i].SEa = Pvac[j][vangul][i-1].SEa + Pvac[j][vangul][i].a*pwrfreq;
                Pvac[j][vangul][i].SEb = Pvac[j][vangul][i-1].SEb + Pvac[j][vangul][i].b*pwrfreq;
                }
          }
          Pvac[j][vtrans][i].SEa = Pvac[j][vtrans][i-1].SEa + Pvac[j][vtrans][i].a*pwrfreq*0.5;
          Pvac[j][vtrans][i].SEb = Pvac[j][vtrans][i-1].SEb + Pvac[j][vtrans][i].b*pwrfreq*0.5;
          if(cmol) {
          Pvac[j][vangul][i].SEa = Pvac[j][vangul][i-1].SEa + Pvac[j][vangul][i].a*pwrfreq*0.5;
          Pvac[j][vangul][i].SEb = Pvac[j][vangul][i-1].SEb + Pvac[j][vangul][i].b*pwrfreq*0.5;
          }
       }
    }

    sprintf(null,"\n Power Spectrum (Density of State) (cm)");
    *outf<<null<<endl; outpwr<<null<<endl;
    sprintf(null," %13s","freq(cm-1)"); *outf<<null; outpwr<<null;
    for(j=0;j<model->ngrp;j++) {  
        if(c2pt==0&&cmol==0) { sprintf(null,"  %6s[G%03d]","PWR",j+1); *outf<<null; outpwr<<null; }
        else {
          if(c2pt) { 
             sprintf(null,"  %6s[G%03d]  %6s[G%03d]","PWR_hs",j+1,"PWR_s",j+1); *outf<<null; outpwr<<null; 
             sprintf(null,"  %6s[G%03d]  %6s[G%03d]","PWR_fr",j+1,"PWR_s",j+1); *outf<<null; outpwr<<null; }
          if(cmol) { sprintf(null,"  %6s[G%03d]  %6s[G%03d]  %6s[G%03d]  %6s[G%03d]","PWRcmt",j+1,"PWRcmo",j+1,"PWRimv",j+1,"PWRcmr",j+1); *outf<<null; outpwr<<null; }
          sprintf(null,"  %6s[G%03d]","PWRtot",j+1); *outf<<null; outpwr<<null;
        }
    }
    *outf<<endl; outpwr<<endl;
    for(i=0;i<nused;i++) { 
       sprintf(null," %13.4f",pwrfreq*i); *outf<<null; outpwr<<null;
       for(j=0;j<model->ngrp;j++) {
          if(c2pt) {  k=0;
            sprintf(null," %13.4f %13.4f",Pvac[j][k][i].a,Pvac[j][k][i].b); *outf<<null; outpwr<<null;
            if(cmol) {
            k=vangul;
            sprintf(null," %13.4f %13.4f",Pvac[j][k][i].a,Pvac[j][k][i].b); *outf<<null; outpwr<<null;
            }
          } 
          for(k=0;k<vactype;k++) {
              sprintf(null," %13.4f",Pvac[j][k][i].max); *outf<<null; outpwr<<null;
          }
       }
       *outf<<endl; outpwr<<endl;
    }

    sprintf(null,"\n Integration of Power Spectrum (cm)");
    *outf<<null<<endl; out3n<<null<<endl;

    sprintf(null," %13s","freq(cm-1)"); *outf<<null; out3n<<null;
    for(j=0;j<model->ngrp;j++) {  
        if(c2pt==0&&cmol==0) { sprintf(null,"  %6s[G%03d]","INT",j+1); *outf<<null; out3n<<null; }
        else {
          if(c2pt) { 
            sprintf(null,"  %6s[G%03d]  %6s[G%03d]","INT_hs",j+1,"INT_s",j+1); *outf<<null; out3n<<null;
            sprintf(null,"  %6s[G%03d]  %6s[G%03d]","INT_fr",j+1,"INT_s",j+1); *outf<<null; out3n<<null;}
          if(cmol) { sprintf(null,"  %6s[G%03d]  %6s[G%03d]  %6s[G%03d]  %6s[G%03d]","INTcmt",j+1,"INTcmo",j+1,"INTimv",j+1,"INTcmr",j+1); *outf<<null; out3n<<null; }
          sprintf(null,"  %6s[G%03d]","INTtot",j+1); *outf<<null; out3n<<null;
        }
    }
    *outf<<endl; out3n<<endl;
    for(i=0;i<nused;i++) { 
       sprintf(null," %13.4f",pwrfreq*i); *outf<<null; out3n<<null;
       for(j=0;j<model->ngrp;j++) {
          if(c2pt) {  k=0;
            sprintf(null," %13.4f %13.4f",Pvac[j][k][i].SEa,Pvac[j][k][i].SEb); *outf<<null; out3n<<null;
            if(cmol) {
            k=vangul;
            sprintf(null," %13.4f %13.4f",Pvac[j][k][i].SEa,Pvac[j][k][i].SEb); *outf<<null; out3n<<null;
            }
          } 
          for(k=0;k<vactype;k++) {
             sprintf(null," %13.4f",Pvac[j][k][i].min); *outf<<null; out3n<<null;
          }
       }
       *outf<<endl; out3n<<endl;
    }

    //vibrational analysis
    double thermo[14][model->ngrp][vactype];
    double scaled_temp; // kb*vacT/(h*vlight*100); //scaled_temp= k*T/h/cspeed, unit 1/cm 
    double u;                                //scaled frequency hcv/kT= hc pwrfreq*i /kT
    double weq,wec,wsq,wsc,waq,wac,wcvq,wcvc;
    double tmpTT;
    weq=wec=wsq=wsc=waq=wac=wcvq=wcvc=0;
    //property for zero frequency
    for(k=0;k<vactype;k++) { 
       i=0;
       for(j=0;j<model->ngrp;j++) { 
          tmpTT=vacT[j][k]; 
          //tmpTT=vacT[j][vtotal]; //now use 2pt temperature for property evaluation
          scaled_temp = kb*tmpTT/(h*vlight*100); 
          if(vacT[j][k]>0) u=pwrfreq*0.5/scaled_temp; 
          else u=0;
          if(c2pt&&k==vtrans) {
             tmpg=Pvac[j][k][i].a; //gas contribution
             tmps=Pvac[j][k][i].b; //solid contribution
          } else if(c2pt&&k==vangul) { 
             tmpg=Pvac[j][k][i].a; //gas contribution
             tmps=Pvac[j][k][i].b; //solid contribution
//adhoc      tmpg+=tmps; tmps=0;
          } else { 
             tmpg=0;                 //gas-like fraction, from 2pt
             tmps=Pvac[j][k][i].max; //solid-like fraction
             if(tmps!=tmps) tmps=Pvac[j][k][i].max=u=0;
          }
          thermo[0][j][k] =0.5*(tmps+tmpg);             //degrees of freedom
          thermo[1][j][k] =0.5*(tmps*0.0);              //ZPE
          thermo[2][j][k] =0.5*(tmps*1.0+tmpg*wep[j]);     //Eq
          thermo[3][j][k] =0.5*(tmps*1.0+tmpg*wep[j]);     //Ec
          thermo[4][j][k] =0.5*(tmps*sqweighting(u)+tmpg*wsp[j]);     //Sq
          thermo[5][j][k] =0.5*(tmps*scweighting(u)+tmpg*wsp[j]);     //Sc
          thermo[6][j][k] =0.5*(tmps*(1-sqweighting(u))+tmpg*wap[j]); //Aq
          thermo[7][j][k] =0.5*(tmps*(1-scweighting(u))+tmpg*wap[j]); //Ac
          thermo[8][j][k] =0.5*(tmps*1.0+tmpg*wcvp[j]);    //Cvq
          thermo[9][j][k] =0.5*(tmps*1.0+tmpg*wcvp[j]);    //Cvc
          thermo[13][j][k]=0.5*(tmps*sqweighting(u)+tmpg*wspd[j]);   //Sdq
          if(c2pt&&k==vangul) {
             thermo[4][j][k] += 0.5*tmpg*(wsr[j]-wsp[j]);    //Sq
             thermo[5][j][k] += 0.5*tmpg*(wsr[j]-wsp[j]);    //Sc
             thermo[6][j][k] += 0.5*tmpg*(war[j]-wap[j]);    //Aq
             thermo[7][j][k] += 0.5*tmpg*(war[j]-wap[j]);    //Ac
          }
       }
    
       for(i=1; i<nused; i++) {
          for(j=0;j<model->ngrp;j++) {
             tmpTT=vacT[j][k];
             //tmpTT=vacT[j][vtotal]; //now use 2pt temperature for property evaluation
             scaled_temp = kb*tmpTT/(h*vlight*100); 
             if(scaled_temp>0) {
               u=pwrfreq*i/scaled_temp; 
               weq=u/2.0+u/(exp(u)-1.0);
               wec=1.0;
               wsq=u/(exp(u)-1)-log(1-exp(-u));
               wsc=1-log(u);
               waq=log((1-exp(-u))/exp(-u/2));
               wac=log(u);
               wcvq=u*u*exp(u)/(1-exp(u))/(1-exp(u));
               wcvc=1.0;
             } else u=weq=wec=wsq=wsc=waq=wac=wcvq=wcvc=0;

             if(c2pt&&k==vtrans) {
                tmpg=Pvac[j][k][i].a; //gas contribution
                tmps=Pvac[j][k][i].b; //solid contribution
             } else if(c2pt&&k==vangul) { 
                tmpg=Pvac[j][k][i].a; //gas contribution
                tmps=Pvac[j][k][i].b; //solid contribution
//adhoc         tmpg+=tmps; tmps=0;
             } else { 
                tmpg=0;                 //gas-like fraction, from 2pt
                tmps=Pvac[j][k][i].max; //solid-like fraction
                if(tmps!=tmps) tmps=u=weq=wec=wsq=wsc=waq=wac=wcvq=wcvc=Pvac[j][k][i].max=0;
             }
             if(i==nused-1) { tmpg*=0.5; tmps*=0.5; }
             thermo[0][j][k] +=(tmps+tmpg);            //degrees of freedom
             thermo[1][j][k] +=(tmps*u*0.5);           //ZPE
             thermo[2][j][k] +=(tmps*weq+tmpg*wep[j]);    //Eq
             thermo[3][j][k] +=(tmps*wec+tmpg*wep[j]);    //Ec
             thermo[4][j][k] +=(tmps*wsq+tmpg*wsp[j]);    //Sq
             thermo[5][j][k] +=(tmps*wsc+tmpg*wsp[j]);    //Sc
             thermo[6][j][k] +=(tmps*waq+tmpg*wap[j]);    //Aq
             thermo[7][j][k] +=(tmps*wac+tmpg*wap[j]);    //Ac
             thermo[8][j][k] +=(tmps*wcvq+tmpg*wcvp[j]);  //Cvq
             thermo[9][j][k] +=(tmps*wcvc+tmpg*wcvp[j]);  //Cvc
             thermo[13][j][k]+=(tmps*wsq+tmpg*wspd[j]);   //Sdq
             if(c2pt&&k==vangul) {
                thermo[4][j][k] += tmpg*(wsr[j]-wsp[j]);    //Sq
                thermo[5][j][k] += tmpg*(wsr[j]-wsp[j]);    //Sc
                thermo[6][j][k] += tmpg*(war[j]-wap[j]);    //Aq
                thermo[7][j][k] += tmpg*(war[j]-wap[j]);    //Ac
             }
          }
       }
    
       for(j=0;j<model->ngrp;j++) {
          tmpTT=vacT[j][k];                        //added 2013.3.28 bug found by Kuan-Yu Yeh
          thermo[0][j][k]*= pwrfreq;               //total degrees of freedom
          thermo[1][j][k]*= pwrfreq*tmpTT*R*1.0e-3; //zero point energy kJ/mol/SimBox
          thermo[2][j][k]*= pwrfreq*tmpTT*R*1.0e-3; //quantum energy    kJ/mol/SimBox
          thermo[3][j][k]*= pwrfreq*tmpTT*R*1.0e-3; //classical energy  kJ/mol/SimBox
          thermo[4][j][k]*= pwrfreq*R;             //quantum entropy   J/(K mol)/SimBox
          thermo[5][j][k]*= pwrfreq*R;             //classical entropy J/(K mol)/SimBox
          thermo[6][j][k]*= pwrfreq*tmpTT*R*1.0e-3; //quantum Helmholtz free energy   kJ/mol/SimBox
          thermo[7][j][k]*= pwrfreq*tmpTT*R*1.0e-3; //classical Helmholtz free energy kJ/mol/SimBox
          thermo[8][j][k]*= pwrfreq*R; //quantum constant volume heat capacity J/(K mol)/SimBox
          thermo[9][j][k]*= pwrfreq*R; //classical constant volume heat capacity J/(K mol)/SimBox
          //thermo[10][j][k]= vacE*vacDF[j][k]/vacDF[model->ngrp-1][vactype-1]; //MD total (strain) energy  kJ/mol/SimBox
          thermo[10][j][k]= model->grp[j].eng*vacDF[j][k]/vacDF[j][vactype-1]; //MD total (strain) energy  kJ/mol/SimBox
          thermo[11][j][k]= thermo[10][j][k]-thermo[3][j][k];   //reference energy Eo kJ/mol/SimBox
          thermo[2][j][k]+= thermo[11][j][k]; //Eq= Eq'+Eo
          thermo[3][j][k]+= thermo[11][j][k]; //Ec= Ec'+Eo
          thermo[6][j][k]+= thermo[11][j][k]; //Aq= Aq'+Eo
          thermo[7][j][k]+= thermo[11][j][k]; //Ac= Ac'+Eo
          thermo[12][j][k]= Pvac[j][k][0].max;//DoS at zero frequency     cm
          thermo[13][j][k]*= pwrfreq*R; //quantum entropy distinguishable particles  J/(K mol)/SimBox
       }
    }
    
    if(1 && c2pt && cmol ) { //Apply 2PT to total system property
      for(j=0;j<model->ngrp;j++) {
         i=vrotat;
         for(l=0;l<14;l++) thermo[l][j][i]=0;
         //ZPE,Eq,Ec,Sq,Sc,Aq,Ac,Cvq,Cvc,Emd,Eo,So
         for(k=0;k<i;k++) {
            for(l=0;l<14;l++) thermo[l][j][i]+=thermo[l][j][k];
         }
         vacT[j][i]=vacT[j][vtotal];
      }
    }

    int show2pt=1;
    sprintf(null,"\n Calculation of Thermodynamic Properties");
    *outf<<null<<endl; outthermo<<null<<endl;
    sprintf(null," %20s","property"); *outf<<null; outthermo<<null;
    for(j=0;j<model->ngrp;j++) {  
        if(c2pt==0&&cmol==0) { sprintf(null," %4s[G%03d]","GRP",j+1); *outf<<null; outthermo<<null; }
        else {
          if(c2pt&&show2pt) { sprintf(null," %4s[G%03d] %4s[G%03d]","Mgas",j+1,"Msol",j+1); *outf<<null; outthermo<<null;}
          if(cmol) { 
             sprintf(null," %4s[G%03d]","Mcmt"); *outf<<null; outthermo<<null;
             if(c2pt&&show2pt) { sprintf(null," %4s[G%03d] %4s[G%03d]","Mgas",j+1,"Msol",j+1); *outf<<null; outthermo<<null;}
             sprintf(null," %4s[G%03d] %4s[G%03d] %4s[G%03d]","Mcmr",j+1,"Mimv",j+1,"M2pt",j+1); 
             *outf<<null; outthermo<<null;
          }
          sprintf(null," %4s[G%03d]","Mtot",j+1); *outf<<null; outthermo<<null;
        }
    }
    *outf<<endl; outthermo<<endl;
    sprintf(null," %20s","nmolecules"); *outf<<null; outthermo<<null;
    for(j=0;j<model->ngrp;j++) {
        k=0; i=model->grp[j].nmol;
        if(c2pt&&show2pt) { 
           double x=hsdf[0][j]/thermo[0][j][k];
           sprintf(null," %10.2f %10.2f",x*i,(1-x)*i); *outf<<null; outthermo<<null; 
        }
        sprintf(null," %10d",i); *outf<<null; outthermo<<null; k=1; 
        if(cmol&&c2pt&&show2pt) { 
           double x=hsdf[1][j]/thermo[0][j][k];
           sprintf(null," %10.2f %10.2f",x*i,(1-x)*i); *outf<<null; outthermo<<null; 
        }
        for(;k<vactype;k++) {
            sprintf(null," %10d",i); *outf<<null; outthermo<<null;
        }
    }
    *outf<<endl; outthermo<<endl;
    sprintf(null," %20s","natom"); *outf<<null; outthermo<<null;
    for(j=0;j<model->ngrp;j++) {
        k=0; i=model->grp[j].natom; 
        if(c2pt&&show2pt) { 
           double x=hsdf[0][j]/thermo[0][j][k];
           sprintf(null," %10.2f %10.2f",x*i,(1-x)*i); *outf<<null; outthermo<<null; 
        }
        sprintf(null," %10d",i); *outf<<null; outthermo<<null; k=1;
        if(cmol&&c2pt&&show2pt) { 
           double x=hsdf[1][j]/thermo[0][j][k];
           sprintf(null," %10.2f %10.2f",x*i,(1-x)*i); *outf<<null; outthermo<<null; 
        }
        for(;k<vactype;k++) {
            sprintf(null," %10d",i); *outf<<null; outthermo<<null;
        }
    }
    *outf<<endl; outthermo<<endl;
    sprintf(null," %20s","dof"); *outf<<null; outthermo<<null;
    for(j=0;j<model->ngrp;j++) {
        k=0;
        if(c2pt&&show2pt) { sprintf(null," %10.2f %10.2f",hsdf[0][j],thermo[0][j][k]-hsdf[0][j]); *outf<<null; outthermo<<null; }
        sprintf(null," %10.2f",thermo[0][j][k]); *outf<<null; outthermo<<null; k=1;
        if(cmol&&c2pt&&show2pt) { sprintf(null," %10.2f %10.2f",hsdf[1][j],thermo[0][j][k]-hsdf[1][j]); *outf<<null; outthermo<<null;} 
        for(;k<vactype;k++) {
            sprintf(null," %10.2f",thermo[0][j][k]); *outf<<null; outthermo<<null;
        }
    }
    *outf<<endl; outthermo<<endl;
    sprintf(null," %20s","temperature_____(K)"); *outf<<null; outthermo<<null;
    for(j=0;j<model->ngrp;j++) {
        k=0;
        if(c2pt&&show2pt) { sprintf(null," %10.2f %10.2f",vacT[j][k],vacT[j][k]); *outf<<null; outthermo<<null; }
        sprintf(null," %10.2f",vacT[j][k]); *outf<<null; outthermo<<null; k=1;
        if(cmol&&c2pt&&show2pt) { sprintf(null," %10.2f %10.2f",vacT[j][k],vacT[j][k]); *outf<<null; outthermo<<null; } 
        for(;k<vactype;k++) {
            sprintf(null," %10.2f",vacT[j][k]); *outf<<null; outthermo<<null;
        }
    }
    *outf<<endl; outthermo<<endl;
    sprintf(null," %20s","pressure______(GPa)"); *outf<<null; outthermo<<null;
    for(j=0;j<model->ngrp;j++) {
        k=0;
        if(c2pt&&show2pt) { sprintf(null," %10.2f %10.2f",vacP,vacP); *outf<<null; outthermo<<null; }
        sprintf(null," %10.2f",vacP); *outf<<null; outthermo<<null; k=1;
        if(cmol&&c2pt&&show2pt) { sprintf(null," %10.2f %10.2f",vacP,vacP); *outf<<null; outthermo<<null; }
        for(;k<vactype;k++) {
            sprintf(null," %10.2f",vacP); *outf<<null; outthermo<<null;
        }
    }
    *outf<<endl; outthermo<<endl;
    sprintf(null," %20s","volume_________(A3)"); *outf<<null; outthermo<<null;
    for(j=0;j<model->ngrp;j++) {
        k=0;
        if(c2pt&&show2pt) { sprintf(null," %10.2f %10.2f",vacPMV[j],vacPMV[j]); *outf<<null; outthermo<<null; }
        sprintf(null," %10.2f",vacPMV[j]); *outf<<null; outthermo<<null; k=1;
        if(cmol&&c2pt&&show2pt) { sprintf(null," %10.2f %10.2f",vacPMV[j],vacPMV[j]); *outf<<null; outthermo<<null; }
        for(;k<vactype;k++) {
            sprintf(null," %10.2f",vacPMV[j]); *outf<<null; outthermo<<null;
        }
    }
    *outf<<endl; outthermo<<endl;
    sprintf(null," %20s","ZPE_(kJ/mol/SimBox)"); *outf<<null; outthermo<<null;
    for(j=0;j<model->ngrp;j++) {
        k=0;
        if(c2pt&&show2pt) { sprintf(null," %10.2f %10.2f",0.0,thermo[1][j][k]); *outf<<null; outthermo<<null; } 
        sprintf(null," %10.2f",thermo[1][j][k]); *outf<<null; outthermo<<null; k=1;
        if(cmol&&c2pt&&show2pt) { sprintf(null," %10.2f %10.2f",0.0,thermo[1][j][k]); *outf<<null; outthermo<<null; }
        for(;k<vactype;k++) {
            sprintf(null," %10.2f",thermo[1][j][k]); *outf<<null; outthermo<<null;
        }
    }
    *outf<<endl; outthermo<<endl;
    sprintf(null," %20s","Emd_(kJ/mol/SimBox)"); *outf<<null; outthermo<<null;
    for(j=0;j<model->ngrp;j++) {
        k=0;
        if(c2pt&&show2pt) { 
           double x=hsdf[0][j]/thermo[0][j][k];
           sprintf(null," %10.2f %10.2f",x*thermo[10][j][k],(1-x)*thermo[10][j][k]); *outf<<null; outthermo<<null;
        }
        sprintf(null," %10.2f",thermo[10][j][k]); *outf<<null; outthermo<<null; k=1;
        if(cmol&&c2pt&&show2pt) { 
           double x=hsdf[1][j]/thermo[0][j][k];
           sprintf(null," %10.2f %10.2f",x*thermo[10][j][k],(1-x)*thermo[10][j][k]); *outf<<null; outthermo<<null;
        } 
        for(;k<vactype;k++) {
            sprintf(null," %10.2f",thermo[10][j][k]); *outf<<null; outthermo<<null;
        }
    }
    *outf<<endl; outthermo<<endl;
    sprintf(null," %20s","Eo__(kJ/mol/SimBox)"); *outf<<null; outthermo<<null;
    for(j=0;j<model->ngrp;j++) {
        k=0;
        if(c2pt&&show2pt) { 
           double x=(hsdf[0][j]/thermo[0][j][k])*thermo[10][j][k]-hsdf[0][j]*wep[j]; 
           sprintf(null," %10.2f %10.2f",x,thermo[11][j][k]-x); *outf<<null; outthermo<<null;
        }
        sprintf(null," %10.2f",thermo[11][j][k]); *outf<<null; outthermo<<null; k=1;
        if(cmol&&c2pt&&show2pt) {
           double x=(hsdf[1][j]/thermo[0][j][k])*thermo[10][j][k]-hsdf[1][j]*wep[j];
           sprintf(null," %10.2f %10.2f",x,thermo[11][j][k]-x); *outf<<null; outthermo<<null;
        }
        for(;k<vactype;k++) {
            sprintf(null," %10.2f",thermo[11][j][k]); *outf<<null; outthermo<<null;
        }
    }
    *outf<<endl; outthermo<<endl;
    sprintf(null," %20s","Eq__(kJ/mol/SimBox)"); *outf<<null; outthermo<<null;
    for(j=0;j<model->ngrp;j++) {
        k=0;
        if(c2pt&&show2pt) { 
           double x=hsdf[0][j]/thermo[0][j][k]*thermo[10][j][k];  //gas_fraction*Eo
           sprintf(null," %10.2f %10.2f",x,thermo[2][j][k]-x); *outf<<null; outthermo<<null;
        } 
        sprintf(null," %10.2f",thermo[2][j][k]); *outf<<null; outthermo<<null; k=1;
        if(cmol&&c2pt&&show2pt) {
           double x=hsdf[1][j]/thermo[0][j][k]*thermo[10][j][k];  //gas_fraction*Eo
           sprintf(null," %10.2f %10.2f",x,thermo[2][j][k]-x); *outf<<null; outthermo<<null;
        }
        for(;k<vactype;k++) {
            sprintf(null," %10.2f",thermo[2][j][k]); *outf<<null; outthermo<<null;
        }
    }
    *outf<<endl; outthermo<<endl;
    sprintf(null," %20s","Ec__(kJ/mol/SimBox)"); *outf<<null; outthermo<<null;
    for(j=0;j<model->ngrp;j++) {
        k=0;
        if(c2pt&&show2pt) { 
           double x=hsdf[0][j]/thermo[0][j][k]*thermo[10][j][k];
           sprintf(null," %10.2f %10.2f",x,thermo[3][j][k]-x); *outf<<null; outthermo<<null;
        } 
        sprintf(null," %10.2f",thermo[3][j][k]); *outf<<null; outthermo<<null; k=1;
        if(cmol&&c2pt&&show2pt) {
           double x=hsdf[1][j]/thermo[0][j][k]*thermo[10][j][k];
           sprintf(null," %10.2f %10.2f",x,thermo[3][j][k]-x); *outf<<null; outthermo<<null;
        }
        for(;k<vactype;k++) {
            sprintf(null," %10.2f",thermo[3][j][k]); *outf<<null; outthermo<<null;
        }
    }
    *outf<<endl; outthermo<<endl;
    sprintf(null," %20s","Sq_(J/mol_K/SimBox)"); *outf<<null; outthermo<<null;
    for(j=0;j<model->ngrp;j++) {
        k=0;
        if(c2pt&&show2pt) { 
           double x=hsdf[0][j]*wsp[j];
           sprintf(null," %10.2f %10.2f",x,thermo[4][j][k]-x); *outf<<null; outthermo<<null;
        }
        sprintf(null," %10.2f",thermo[4][j][k]); *outf<<null; outthermo<<null; k=1;
        if(cmol&&c2pt&&show2pt) {
           double x=hsdf[1][j]*wsr[j];
           sprintf(null," %10.2f %10.2f",x,thermo[4][j][k]-x); *outf<<null; outthermo<<null;
        }
        for(;k<vactype;k++) {
            sprintf(null," %10.2f",thermo[4][j][k]); *outf<<null; outthermo<<null;
        }
    }
    *outf<<endl; outthermo<<endl;
    sprintf(null," %20s","Sc_(J/mol_K/SimBox)"); *outf<<null; outthermo<<null;
    for(j=0;j<model->ngrp;j++) {
        k=0;
        if(c2pt&&show2pt) { 
           double x=hsdf[0][j]*wsp[j];
           sprintf(null," %10.2f %10.2f",x,thermo[5][j][k]-x); *outf<<null; outthermo<<null;
        }
        sprintf(null," %10.2f",thermo[5][j][k]); *outf<<null; outthermo<<null; k=1;
        if(cmol&&c2pt&&show2pt) {
           double x=hsdf[1][j]*wsr[j];
           sprintf(null," %10.2f %10.2f",x,thermo[5][j][k]-x); *outf<<null; outthermo<<null;
        }
        for(;k<vactype;k++) {
            sprintf(null," %10.2f",thermo[5][j][k]); *outf<<null; outthermo<<null;
        }
    }
    *outf<<endl; outthermo<<endl;
    sprintf(null," %20s","Aq__(kJ/mol/SimBox)"); *outf<<null; outthermo<<null;
    for(j=0;j<model->ngrp;j++) {
        k=0;
        if(c2pt&&show2pt) { 
           double x=hsdf[0][j]/thermo[0][j][k]*thermo[10][j][k]-vacT[j][k]*hsdf[0][j]*wsp[j]*1E-3;
           sprintf(null," %10.2f %10.2f",x,thermo[6][j][k]-x); *outf<<null; outthermo<<null;
        } 
        sprintf(null," %10.2f",thermo[6][j][k]); *outf<<null; outthermo<<null; k=1;
        if(cmol&&c2pt&&show2pt) {
           double x=hsdf[1][j]/thermo[0][j][k]*thermo[10][j][k]-vacT[j][k]*hsdf[1][j]*wsr[j]*1E-3;
           sprintf(null," %10.2f %10.2f",x,thermo[6][j][k]-x); *outf<<null; outthermo<<null;
        }
        for(;k<vactype;k++) {
            sprintf(null," %10.2f",thermo[6][j][k]); *outf<<null; outthermo<<null;
        }
    }
    *outf<<endl; outthermo<<endl;
    sprintf(null," %20s","Ac__(kJ/mol/SimBox)"); *outf<<null; outthermo<<null;
    for(j=0;j<model->ngrp;j++) {
        k=0;
        if(c2pt&&show2pt) { 
           double x=hsdf[0][j]/thermo[0][j][k]*thermo[10][j][k]-vacT[j][k]*hsdf[0][j]*wsp[j]*1E-3;
           sprintf(null," %10.2f %10.2f",x,thermo[7][j][k]-x); *outf<<null; outthermo<<null;
        }
        sprintf(null," %10.2f",thermo[7][j][k]); *outf<<null; outthermo<<null; k=1;
        if(cmol&&c2pt&&show2pt) {
           double x=hsdf[1][j]/thermo[0][j][k]*thermo[10][j][k]-vacT[j][k]*hsdf[1][j]*wsr[j]*1E-3;
           sprintf(null," %10.2f %10.2f",x,thermo[7][j][k]-x); *outf<<null; outthermo<<null;
        }
        for(;k<vactype;k++) {
            sprintf(null," %10.2f",thermo[7][j][k]); *outf<<null; outthermo<<null;
        }
    }
    *outf<<endl; outthermo<<endl;
    sprintf(null," %20s","Cvq(J/mol_K/SimBox)"); *outf<<null; outthermo<<null;
    for(j=0;j<model->ngrp;j++) {
        k=0;
        if(c2pt&&show2pt) { 
           double x=hsdf[0][j]*wcvp[j];
           sprintf(null," %10.2f %10.2f",x,thermo[8][j][k]-x); *outf<<null; outthermo<<null;
        } 
        sprintf(null," %10.2f",thermo[8][j][k]); *outf<<null; outthermo<<null; k=1;
        if(cmol&&c2pt&&show2pt) {
           double x=hsdf[1][j]*wcvp[j];
           sprintf(null," %10.2f %10.2f",x,thermo[8][j][k]-x); *outf<<null; outthermo<<null;
        }
        for(;k<vactype;k++) {
            sprintf(null," %10.2f",thermo[8][j][k]); *outf<<null; outthermo<<null;
        }
    }
    *outf<<endl; outthermo<<endl;
    sprintf(null," %20s","Cvc(J/mol_K/SimBox)"); *outf<<null; outthermo<<null;
    for(j=0;j<model->ngrp;j++) {
        k=0;
        if(c2pt&&show2pt) { 
           double x=hsdf[0][j]*wcvp[j];
           sprintf(null," %10.2f %10.2f",x,thermo[9][j][k]-x); *outf<<null; outthermo<<null;
        } 
        sprintf(null," %10.2f",thermo[9][j][k]); *outf<<null; outthermo<<null; k=1;
        if(cmol&&c2pt&&show2pt) {
           double x=hsdf[1][j]*wcvp[j];
           sprintf(null," %10.2f %10.2f",x,thermo[9][j][k]-x); *outf<<null; outthermo<<null;
        }
        for(;k<vactype;k++) {
            sprintf(null," %10.2f",thermo[9][j][k]); *outf<<null; outthermo<<null;
        }
    }
    *outf<<endl; outthermo<<endl;
    sprintf(null," %20s","S(0)(cm/mol/SimBox)"); *outf<<null; outthermo<<null;
    for(j=0;j<model->ngrp;j++) {
        k=0;
        if(c2pt&&show2pt) { 
           double x=Pvac[j][k][0].a;
           sprintf(null," %10.2f %10.2f",x,thermo[12][j][k]-x); *outf<<null; outthermo<<null;
        } 
        sprintf(null," %10.2f",thermo[12][j][k]); *outf<<null; outthermo<<null; k=1;
        if(cmol&&c2pt&&show2pt) {
           double x=Pvac[j][k][0].a;
           sprintf(null," %10.2f %10.2f",x,thermo[12][j][k]-x); *outf<<null; outthermo<<null;
        }
        for(;k<vactype;k++) {
            sprintf(null," %10.2f",thermo[12][j][k]); *outf<<null; outthermo<<null;
        }
    }
    *outf<<endl; outthermo<<endl;
    sprintf(null," %20s","Diffus(cm2/s_,_1/s)"); *outf<<null; outthermo<<null;
    for(j=0;j<model->ngrp;j++) {
        k=0;
        if(c2pt&&show2pt) { 
           double x=(Pvac[j][k][0].a*R*vacT[j][k])/(12.0*vlight*model->grp[j].mass)*1.0E5;
           sprintf(null," %10.4e %10.4e",x,0.0); *outf<<null; outthermo<<null;
        }
        double x=(thermo[12][j][k]*R*vacT[j][k])/(12.0*vlight*model->grp[j].mass)*1.0E5;
        sprintf(null," %10.4e",x); *outf<<null; outthermo<<null; k=1;
        if(cmol&&c2pt&&show2pt) {
           double x=(Pvac[j][k][0].a*R*vacT[j][k])/(12.0*vlight*model->grp[j].mass)*1.0E5;
           sprintf(null," %10.4e %10.4e",x,0.0); *outf<<null; outthermo<<null;
        }
        for(;k<vactype;k++) {
            double x=(thermo[12][j][k]*R*vacT[j][k])/(12.0*vlight*model->grp[j].mass)*1.0E5;
            sprintf(null," %10.4e",x); *outf<<null; outthermo<<null;
        }
    }
    *outf<<endl; outthermo<<endl;

    if(c2pt) {
    sprintf(null," %20s","fluidicity_________"); *outf<<null; outthermo<<null;
    for(j=0;j<model->ngrp;j++) {
        k=0;
        if(c2pt&&show2pt) { 
           sprintf(null," %10.8f %10.8f",f2pt[0][j],0.0); *outf<<null; outthermo<<null;
        } 
        sprintf(null," %10.8f",f2pt[0][j]); *outf<<null; outthermo<<null; k=1;
        if(cmol&&c2pt&&show2pt) {
           sprintf(null," %10.8f %10.8f",f2pt[1][j],0.0); *outf<<null; outthermo<<null;
        }
        for(;k<vactype;k++) {
            double x=0.0;
            if(k==vtrans) x=f2pt[0][j];
            else if(k==vangul) x=f2pt[1][j];
            sprintf(null," %10.8f",x); *outf<<null; outthermo<<null;
        }
    }
    *outf<<endl; outthermo<<endl;
    }

    if(c2ptmf) {
    sprintf(null," %20s","fluidicitymf_______"); *outf<<null; outthermo<<null;
    for(j=0;j<model->ngrp;j++) {
        k=0;
        if(c2pt&&show2pt) {
           sprintf(null," %10.8f %10.8f",fmf[0][j],0.0); *outf<<null; outthermo<<null;
        }
        sprintf(null," %10.8f",fmf[0][j]); *outf<<null; outthermo<<null; k=1;
        if(cmol&&c2pt&&show2pt) {
           sprintf(null," %10.8f %10.8f",fmf[1][j],0.0); *outf<<null; outthermo<<null;
        }
        for(;k<vactype;k++) {
            double x=0.0;
            if(k==vtrans) x=fmf[0][j];
            else if(k==vangul) x=fmf[1][j];
            sprintf(null," %10.8f",x); *outf<<null; outthermo<<null;
        }
    }
    *outf<<endl; outthermo<<endl;
    }

    sprintf(null," %20s","RotSymNum__________"); *outf<<null; outthermo<<null;
    for(j=0;j<model->ngrp;j++) {
        k=0;
        if(c2pt&&show2pt) {
           sprintf(null," %10d %10d",model->grp[j].rotsym,model->grp[j].rotsym); *outf<<null; outthermo<<null;
        }
        sprintf(null," %10d",model->grp[j].rotsym); *outf<<null; outthermo<<null; k=1;
        if(cmol&&c2pt&&show2pt) {
           sprintf(null," %10d %10d",model->grp[j].rotsym,model->grp[j].rotsym); *outf<<null; outthermo<<null;
        }
        for(;k<vactype;k++) {
            sprintf(null," %10d",model->grp[j].rotsym); *outf<<null; outthermo<<null;
        }
    }
    *outf<<endl; outthermo<<endl;
    sprintf(null," %20s","LinearMoleculeFlag_"); *outf<<null; outthermo<<null;
    for(j=0;j<model->ngrp;j++) {
        k=0;
        if(c2pt&&show2pt) {
           sprintf(null," %10d %10d",model->grp[j].linear,model->grp[j].linear); *outf<<null; outthermo<<null;
        }
        sprintf(null," %10d",model->grp[j].linear); *outf<<null; outthermo<<null; k=1;
        if(cmol&&c2pt&&show2pt) {
           sprintf(null," %10d %10d",model->grp[j].linear,model->grp[j].linear); *outf<<null; outthermo<<null;
        }
        for(;k<vactype;k++) {
            sprintf(null," %10d",model->grp[j].linear); *outf<<null; outthermo<<null;
        }
    }
    *outf<<endl; outthermo<<endl;
    sprintf(null," %20s","PMVRatio___________"); *outf<<null; outthermo<<null;
    for(j=0;j<model->ngrp;j++) {
        k=0;
        if(c2pt&&show2pt) {
           sprintf(null," %10.4e %10.4e",model->grp[j].pmvr,model->grp[j].pmvr); *outf<<null; outthermo<<null;
        }
        sprintf(null," %10.8f",model->grp[j].pmvr); *outf<<null; outthermo<<null; k=1;
        if(cmol&&c2pt&&show2pt) {
           sprintf(null," %10.4e %10.4e",model->grp[j].pmvr,model->grp[j].pmvr); *outf<<null; outthermo<<null;
        }
        for(;k<vactype;k++) {
            sprintf(null," %10.8f",model->grp[j].pmvr); *outf<<null; outthermo<<null;
        }
    }
    *outf<<endl; outthermo<<endl;


    model->atom=tatom;
    if(cmol) {
       for(i=0;i<model->nmol;i++) {
           for(j=0;j<model->mol[i].natom;j++) model->mol[i].atom[j] = & tatom[ model->mol[i].atm[j] ];
       }
    }
    for(i=0;i<model->trj.ntrjf;i++) {
        model->trj.strj[i].atom=model->atom;
        model->trj.strj[i].cell=&(model->cell);
        model->trj.strj[i].prp =&(model->prp);
    }

    outvac.close();
    outpwr.close();
    out3n.close();
    outthermo.close();

    return 0;
}

