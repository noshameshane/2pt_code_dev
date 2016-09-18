#include <iostream>
#include <fstream>
#include <cmath>
#include "structure.h"
#include "utility.h"

int AGRID::set_grid(double insgrid)
{
    if(sgrid==insgrid) return 0; //grids have already been set
    sgrid=insgrid;
    int i,j,k,cc,dd,li;
    free_grid();
    ngrid=(int)(radius/sgrid)+1;  
    vgrid=(char ***) calloc(ngrid,sizeof(char **));
    for(i=0;i<ngrid;i++) {
       vgrid[i]=(char **) calloc(ngrid,sizeof(char *));
       for(j=0;j<ngrid;j++) {
          vgrid[i][j]=(char *) calloc(ngrid,sizeof(char));
          for(k=0;k<ngrid;k++) vgrid[i][j][k]=0;
       }
    }
    double dist2,r2,s2;
    r2=radius*radius;
    s2=sgrid*sgrid;
    dd=cc=0;
    switch(4) {
      case 1: //brute force, finding occupied grids
      for(i=0;i<ngrid;i++) {
         for(j=0;j<ngrid;j++) {
            for(k=0;k<ngrid;k++) {
               dist2=(i*i+j*j+k*k)*s2; dd++;
               if(dist2<r2) { vgrid[i][j][k]=1; cc++; }
            }
         }
      }
      break; 
      case 2: // apply symmetry
      for(i=0;i<ngrid;i++) {
         for(j=0;j<ngrid;j++) {
            for(k=0;k<=j;k++) {
               dist2=(i*i+j*j+k*k)*s2; dd++;
               if(dist2<r2) { vgrid[i][j][k]=vgrid[i][k][j]=1; (j==k?cc++:cc+=2); }
            }
         }
      }
      break;
      case 3: // apply symmetry + certain voids
      for(i=0;i<ngrid;i++) {
         for(j=0;j<ngrid;j++) {
            for(k=0;k<=j;k++) {
               dist2=(i*i+j*j+k*k)*s2; dd++;
               if(dist2<r2) { vgrid[i][j][k]=vgrid[i][k][j]=1; (j==k?cc++:cc+=2); }
               else k=j; //these will certainly be empty
            }
         }
      }
      break;
      case 4: // apply symmetry + certain voids + certain occupied
      li=ngrid-1;
      for(i=ngrid-1;i>=0;i--) {
         for(j=0;j<ngrid;j++) {
            for(k=0;k<=j;k++) {
               if(vgrid[li][j][k]==1) { vgrid[i][j][k]=vgrid[i][k][j]=1; (j==k?cc++:cc+=2); }
               else { 
                 dist2=(i*i+j*j+k*k)*s2; dd++;
                 if(dist2<r2) { vgrid[i][j][k]=vgrid[i][k][j]=1; (j==k?cc++:cc+=2); }
                 else k=j; //these will certainly be empty
               }
            }
         }
         li=i;
      }
      break;
    }
    //for(i=0;i<ngrid;i++) for(j=0;j<ngrid;j++) for(k=0;k<ngrid;k++) printf("%d %d %d %d \n",i,j,k,(int)vgrid[i][j][k]);
    //printf("radius %f gs %f ngrid %d tot %d occupied %d fraction %f dist calc %d\n",radius,sgrid,ngrid,ngrid*ngrid*ngrid,cc,cc*1.0/(ngrid*ngrid*ngrid),dd);

    return 0;
}

int AGRID::free_grid()
{
    if(vgrid!=NULL) {
       int i,j;
       for(i=0;i<ngrid;i++) {
          for(j=0;j<ngrid;j++) free(vgrid[i][j]);
          free(vgrid[i]);
       }
       free(vgrid); vgrid=NULL;
    }
    return 0;
}

int MOLECULE::init_atom()
{
    if(natom==0) return 1;
    if(atm!=NULL) delete [] atm;
    atm=new int [natom];
    atom=new ATOM *[natom];
    return 0;
}

int MOLECULE::cal_mass()
{
    int i;
    mass=0;
    for(i=0;i<natom;i++) mass += atom[i]->mass;
    return 0;
}

int MOLECULE::cal_chgmu()
{
    int i,k;
    chg=mu[0]=mu[1]=mu[2]=0;
    for(i=0;i<natom;i++) {
        chg += atom[i]->chg;
        for(k=0;k<3;k++) mu[k] += (atom[i]->chg)*(atom[i]->pv[k])*eA2debye;
    }
    return 0;
}

int MOLECULE::find_cm()
{
    int i,k;
    mass=pv[0]=pv[1]=pv[2]=0;
    for(i=0;i<natom;i++) {
        mass += atom[i]->mass;
        for(k=0;k<3;k++) pv[k] += (atom[i]->mass)*(atom[i]->pv[k]);
    }
    for(k=0;k<3;k++) pv[k]/=mass;
    return 0;
}

int MOLECULE::cal_vc()
{
// translation = center of mass motion
// M cmvel_k = sum_i sum_k ( mi vi_k )
//
// angular velocity
// omega = sum_i ( ri x vi )
// vrot = omega x ri
//
// vibration
// vib = vi - cmvel - vrot

    int i,j,k;
    double m;

    if(natom==1) {
      i=0;
      for(k=0;k<3;k++) {
         atom[i]->vt[k]=atom[i]->vel[k];
         atom[i]->vr[k]=atom[i]->vv[k]=anguv[k]=0;
         pI[k]=-1;
      }
      return 0;
    }

    for(i=0;i<natom;i++) for(k=0;k<3;k++)  atom[i]->vt[k]=atom[i]->vr[k]=atom[i]->vv[k]=0; 
     
    for(i=0;i<3;i++) {
       omega[i]=angmom[i]=0;
       for(k=0;k<3;k++) inertia[i][k]=0;
       pv[i]=vel[i]=0;
    }

    if(0) {
      printf("atomic coordinates (A)\n");
      for(i=0;i<natom;i++) printf("%d %8.4f %8.4f %8.4f\n",i+1,atom[i]->pv[0],atom[i]->pv[1],atom[i]->pv[2]);
      printf("atomic velocities (A/ps)\n");
      for(i=0;i<natom;i++) printf("%d %8.4f %8.4f %8.4f\n",i+1,atom[i]->vel[0],atom[i]->vel[1],atom[i]->vel[2]);
    }
    
    //determine cm pv and vel
    for(i=0;i<natom;i++) {
        for(k=0;k<3;k++) {
            m=atom[i]->mass;
            pv[k]  += m*(atom[i]->pv[k]);
            vel[k] += m*(atom[i]->vel[k]);
        }
    }
    for(k=0;k<3;k++) { pv[k]/=mass; vel[k]/=mass;}

    double relpv[natom][3],relvel[natom][3];
    //set the center of mass to the origin
    for(i=0;i<natom;i++) {
        for(k=0;k<3;k++){
            relpv[i][k] = atom[i]->pv[k]-pv[k];
            relvel[i][k]= atom[i]->vel[k]-vel[k];
            atom[i]->vt[k]=vel[k];  //translational velocity
        }
    }

    if(0) {
      printf("translation velocities\n");
      printf("cm %8.4f %8.4f %8.4f\n",vel[0],vel[1],vel[2]);
      for(i=0;i<natom;i++) printf("%d %8.4f %8.4f %8.4f\n",i+1,atom[i]->vt[0],atom[i]->vt[1],atom[i]->vt[2]);
    }

    //calculate angular momentum 
    for(i=0;i<natom;i++) { // L = sum_i mi (ri x vi) 
        m=atom[i]->mass;
        angmom[0] += m*(relpv[i][1]*relvel[i][2]-relpv[i][2]*relvel[i][1]);
        angmom[1] += m*(relpv[i][2]*relvel[i][0]-relpv[i][0]*relvel[i][2]);
        angmom[2] += m*(relpv[i][0]*relvel[i][1]-relpv[i][1]*relvel[i][0]);
    }

    //calculate inertia tensor of molecule
    for(i=0;i<natom;i++) {
        m=atom[i]->mass;
        inertia[0][0] += m*(pow(relpv[i][1],2.0) + pow(relpv[i][2],2.0));
        inertia[1][1] += m*(pow(relpv[i][0],2.0) + pow(relpv[i][2],2.0));
        inertia[2][2] += m*(pow(relpv[i][0],2.0) + pow(relpv[i][1],2.0));
        inertia[0][1] -= m* relpv[i][0] * relpv[i][1];
        inertia[0][2] -= m* relpv[i][0] * relpv[i][2];
        inertia[1][2] -= m* relpv[i][1] * relpv[i][2];
    }
    inertia[1][0] = inertia[0][1];
    inertia[2][0] = inertia[0][2];
    inertia[2][1] = inertia[1][2];

    if(0) {
      printf("\nInertia Tensor\n");
      for(i=0;i<3;i++) {
          for(k=0;k<3;k++) printf("%16.8lf",inertia[i][k]);
              printf("\n");
      }
      for(i=0;i<natom;i++) {
         printf("atom %d mass %f pv %f %f %f\n",i+1,atom[i]->mass,atom[i]->pv[0],atom[i]->pv[1],atom[i]->pv[2]);
      }
    }

    //calculate principle moments of inertia
    double evl[3],**evt,pomega[3];
    evt=new double *[3]; for(i=0;i<3;i++) evt[i]=new double [3];
    for(i=0;i<3;i++) evl[i]=evt[i][0]=evt[i][1]=evt[i][2]=pomega[i]=0;
    jacobi(inertia,3,evl,evt); //inertia tensor,dimension,eigenvalues,eigenvectors
    eigsrt(evl,evt,3); //sort eigenvalues
    if( linear ) evl[2]=0; //stlin 2010/08/04, added to correct for flexible CO2
    for(i=0;i<3;i++) pI[i]=evl[i]/(Na*1e23); //principle moment of inertia
    //for(i=0;i<3;i++) rotT[i]=h*h*Na*1e23/(8.0*PI*PI*evl[i]*kb); rotational temperatures

    if(0) {
      printf("\nPrinciple moments\n");
      for(i=0;i<3;i++) printf("%16.8lf vector %16.8lf %16.8lf %16.8lf\n",evl[i],evt[0][i],evt[1][i],evt[2][i]);
      printf("\n");
      printf("\nRotational temperatures\n");
      for(i=0;i<3;i++) printf("I %16.8lf (g/mol*A2) theta %16.8lf K\n",evl[i],h*h/(8.0*PI*PI*pI[i]*kb));
      printf("\n");
    }

    //calculate angular velocities along principle axises
    for(i=0;i<3;i++) pomega[i]=0;
    for(i=0;i<3;i++) for(k=0;k<3;k++) { if(evl[i]>0) pomega[i]+=angmom[k]*evt[k][i]/evl[i]; }

    //calculate inertia weighted angular velocities
    for(i=0;i<3;i++) {
        anguv[i]=omega[i]=0;
        for(k=0;k<3;k++) {
           if(evl[k]>0) {
              //angular velocity
              omega[i] += pomega[k]*evt[i][k];
              //angular vel weighted by principle moments of inertia 
              anguv[i] += pomega[k]*evt[i][k]*sqrt(evl[k]);
           }

        }
    }

    if(0) {
      printf("angular velocities along principle axises\n");
      printf("w %8.4f %8.4f %8.4f\n",pomega[0],pomega[1],pomega[2]);
      for(i=0;i<3;i++) printf("%8.4f v %8.4f %8.4f %8.4f\n",pomega[i],pomega[i]*evt[0][i],pomega[i]*evt[1][i],pomega[i]*evt[2][i]);
    }

    //calculate velocity due to rotation
    for(i=0;i<natom;i++) { //vr = w x r
        atom[i]->vr[0]=(omega[1]*relpv[i][2]-omega[2]*relpv[i][1]);
        atom[i]->vr[1]=(omega[2]*relpv[i][0]-omega[0]*relpv[i][2]);
        atom[i]->vr[2]=(omega[0]*relpv[i][1]-omega[1]*relpv[i][0]);
    }

    if(0) {
      printf("angular velocities\n");
      printf("w %8.4f %8.4f %8.4f\n",omega[0],omega[1],omega[2]);
      printf("weighted angular velocities\n");
      printf("w %8.4f %8.4f %8.4f\n",anguv[0],anguv[1],anguv[2]);
      printf("rotational velocities\n");
      for(i=0;i<natom;i++) printf("%d %8.4f %8.4f %8.4f\n",i+1,atom[i]->vr[0],atom[i]->vr[1],atom[i]->vr[2]);
    }

    //calculate vibrational velocities
    for(i=0;i<natom;i++) {
       for(k=0; k<3; k++) atom[i]->vv[k]= relvel[i][k]-(atom[i]->vr[k]);
    }

    //calculate mass weighted velocities
    if(iccnt!=NULL) {
       double **Zhp,**Zhm,**Zdev; //the Zmatrix [natom][3]
       double **B,**G,**Gi,**Mi;       //[3natom][3natom]
       double vvib[3*natom],qvib[3*natom];  //vibrational velocity
       double det;

       new_2d_matrix(Zhp,natom,3);
       new_2d_matrix(Zhm,natom,3);
       new_2d_matrix(Zdev,natom,3);
       new_2d_matrix(B,3*natom,3*natom);  //
       new_2d_matrix(G,3*natom,3*natom);  //Wilson G matrix
       new_2d_matrix(Gi,3*natom,3*natom); //inverse of G matrix
       new_2d_matrix(Mi,3*natom,3*natom); //inverse of mass tensor

       for(i=0;i<3*natom;i++) {
          k=i/3;
          B[i][i]=G[i][i]=0;
          Mi[i][i]=1.0/atom[k]->mass;
          vvib[i]=atom[k]->vv[i%3];
          for(j=i+1;j<3*natom;j++) B[i][j]=B[j][i]=G[i][j]=G[j][i]=Mi[i][j]=Mi[j][i]=0;
       }
       
       if(0) {
          cal_Zmat(relpv,Zhp);
          cout<<"Zmatrix"<<endl;
          for(i=0;i<natom;i++) printf("%4d %10.4f %10.4f %10.4f\n",i+1,Zhp[i][0],Zhp[i][1],Zhp[i][2]);
       }

       //Mass matrix

       double dh=1e-6;
       j=0;
       for(i=0;i<natom;i++) {
          for(k=0; k<3; k++) {
              relpv[i][k]+=dh;
              cal_Zmat(relpv,Zhp);
              relpv[i][k]-=2*dh;
              cal_Zmat(relpv,Zhm);
              relpv[i][k]+=dh;
              matrix_sub(Zhp,Zhm,natom,3,Zdev);
              matrix_scale(Zdev,natom,3,1.0/(2*dh));
              B[0][j]=Zdev[1][0];
              B[1][j]=Zdev[2][0];
              B[2][j]=Zdev[2][1];
              int ii,jj,kk; jj=3;
              for(ii=3;ii<natom;ii++) for(kk=0;kk<3;kk++) B[jj++][j]=Zdev[ii][kk]; 
              j++;
          }
       }

       for(i=0;i<natom;i++) {
           double mm=atom[i]->mass; 
           B[3*natom-6][3*i  ]=mm/sqrt(mass);
           B[3*natom-5][3*i+1]=mm/sqrt(mass);
           B[3*natom-4][3*i+2]=mm/sqrt(mass);
           B[3*natom-3][3*i  ]=0;
           B[3*natom-3][3*i+1]=-mm*relpv[i][2]/sqrt(evl[2]);
           B[3*natom-3][3*i+2]= mm*relpv[i][1]/sqrt(evl[2]);
           B[3*natom-2][3*i  ]= mm*relpv[i][2]/sqrt(evl[1]);
           B[3*natom-2][3*i+1]=0;
           B[3*natom-2][3*i+2]=-mm*relpv[i][0]/sqrt(evl[1]);
           B[3*natom-1][3*i  ]=-mm*relpv[i][1]/sqrt(evl[0]);
           B[3*natom-1][3*i+1]= mm*relpv[i][0]/sqrt(evl[0]);
           B[3*natom-1][3*i+2]=0;
       }

       if(0) {
           cout<<"MMMMMMMMMMMMM"<<endl;
           for(i=0;i<3*natom;i++) {
               for(j=0;j<3*natom;j++) printf("%10.4f ",B[i][j]);
               //for(j=0;j<3*natom;j++) printf("%10.4f ",M[i][j]);
               cout<<endl;
           } 
       }

       matrix_product(B,3*natom,3*natom,vvib,3*natom,qvib); //qvib=B*vvib
       matrix_BMBt(B,3*natom,Mi,G); //G=B*Mi*Bt
       matrix_inverse(G,3*natom,Gi,det);
       matrix_sqrt(Gi,3*natom,G);
       matrix_product(G,3*natom,3*natom,qvib,3*natom,mqvib); //mqvib=G*qvib

       if(linear) k=5; else k=6;
       for(i=0;i<3*natom-k;i++) { 
          moi[i]=mqvib[i]*mqvib[i]/(qvib[i]*qvib[i]);  //moment of inertia
       }
       for(;i<3*natom;i++) moi[i]=0;

       if(0) {
         printf("%5s %10s %10s\n","ID","qvib","moi");
         for(i=0;i<3*natom;i++)
            printf("%5d %10.4f %10.4f\n",i+1,mqvib[i],moi[i]);
       }
       delete_2d_matrix(Zhp,natom,3);
       delete_2d_matrix(Zhm,natom,3);
       delete_2d_matrix(Zdev,natom,3);
       delete_2d_matrix(B,3*natom,3*natom);
       delete_2d_matrix(G,3*natom,3*natom);
       delete_2d_matrix(Gi,3*natom,3*natom);
       delete_2d_matrix(Mi,3*natom,3*natom);
       
    }

    if(0) {
      printf("vibration velocities\n");
      for(i=0;i<natom;i++) printf("%d %8.4f %8.4f %8.4f\n",i+1,atom[i]->vv[0],atom[i]->vv[1],atom[i]->vv[2]);
    }

    if(0) {
      printf("translational velocities\n");
      for(i=0;i<natom;i++) printf("%d %8.4f %8.4f %8.4f\n",i+1,atom[i]->vt[0],atom[i]->vt[1],atom[i]->vt[2]);
      printf("rotational velocities\n");
      for(i=0;i<natom;i++) printf("%d %8.4f %8.4f %8.4f\n",i+1,atom[i]->vr[0],atom[i]->vr[1],atom[i]->vr[2]);
      printf("vibrational velocities\n");
      for(i=0;i<natom;i++) printf("%d %8.4f %8.4f %8.4f\n",i+1,atom[i]->vv[0],atom[i]->vv[1],atom[i]->vv[2]);
    }
    for(i=0;i<3;i++) delete [] evt[i];
    delete [] evt;

    return 0;
}

int MOLECULE::cal_vc1()
{
// translation = center of mass motion
// M cmvel_k = sum_i sum_k ( mi vi_k )
//
// rotation
// L = sum_i ( mi ri x vi ) = I omega
// vrot = omega x ri
//
// vibration
// vib = vi - cmvel - vrot

    int i,j,k;
    double m;

    if(natom==1) {
      i=0;
      for(k=0;k<3;k++) {
         atom[i]->vt[k]=atom[i]->vel[k];
         atom[i]->vr[k]=atom[i]->vv[k]=anguv[k]=0;
         pI[k]=-1;
      }
      return 0;
    }

    for(i=0;i<natom;i++) for(k=0;k<3;k++)  atom[i]->vt[k]=atom[i]->vr[k]=atom[i]->vv[k]=0; 
     
    for(i=0;i<3;i++) {
       omega[i]=angmom[i]=0;
       for(k=0;k<3;k++) inertia[i][k]=0;
       pv[i]=vel[i]=0;
    }

    if(0) {
      printf("atomic coordinates (A)\n");
      for(i=0;i<natom;i++) printf("%d %8.4f %8.4f %8.4f\n",i+1,atom[i]->pv[0],atom[i]->pv[1],atom[i]->pv[2]);
      printf("atomic velocities (A/ps)\n");
      for(i=0;i<natom;i++) printf("%d %8.4f %8.4f %8.4f\n",i+1,atom[i]->vel[0],atom[i]->vel[1],atom[i]->vel[2]);
    }
    
    //determine cm pv and vel
    for(i=0;i<natom;i++) {
        for(k=0;k<3;k++) {
            m=atom[i]->mass;
            pv[k]  += m*(atom[i]->pv[k]);
            vel[k] += m*(atom[i]->vel[k]);
        }
    }
    for(k=0;k<3;k++) { pv[k]/=mass; vel[k]/=mass;}

    double relpv[natom][3],relvel[natom][3];
    //set the center of mass to the origin
    for(i=0;i<natom;i++) {
        for(k=0;k<3;k++){
            relpv[i][k] = atom[i]->pv[k]-pv[k];
            relvel[i][k]= atom[i]->vel[k]-vel[k];
            atom[i]->vt[k]=vel[k];  //translational velocity
        }
    }

    if(0) {
      printf("translation velocities\n");
      printf("cm %8.4f %8.4f %8.4f\n",vel[0],vel[1],vel[2]);
      for(i=0;i<natom;i++) printf("%d %8.4f %8.4f %8.4f\n",i+1,atom[i]->vt[0],atom[i]->vt[1],atom[i]->vt[2]);
    }

    //calculate inertia tensor of molecule
    for(i=0;i<natom;i++) {
        m=atom[i]->mass;
        inertia[0][0] += m*(pow(relpv[i][1],2.0) + pow(relpv[i][2],2.0));
        inertia[1][1] += m*(pow(relpv[i][0],2.0) + pow(relpv[i][2],2.0));
        inertia[2][2] += m*(pow(relpv[i][0],2.0) + pow(relpv[i][1],2.0));
        inertia[0][1] -= m* relpv[i][0] * relpv[i][1];
        inertia[0][2] -= m* relpv[i][0] * relpv[i][2];
        inertia[1][2] -= m* relpv[i][1] * relpv[i][2];
    }
    inertia[1][0] = inertia[0][1];
    inertia[2][0] = inertia[0][2];
    inertia[2][1] = inertia[1][2];

    if(0) {
      printf("\nInertia Tensor\n");
      for(i=0;i<3;i++) {
          for(k=0;k<3;k++) printf("%16.8lf",inertia[i][k]);
              printf("\n");
      }
      for(i=0;i<natom;i++) {
         printf("atom %d mass %f pv %f %f %f\n",i+1,atom[i]->mass,atom[i]->pv[0],atom[i]->pv[1],atom[i]->pv[2]);
      }
    }
    //Calculate inverse of inertia tensor
    double inv_inr[3][3];
    double det;
    det= inertia[0][0]*(inertia[2][2]*inertia[1][1]-inertia[2][1]*inertia[1][2])
        -inertia[1][0]*(inertia[2][2]*inertia[0][1]-inertia[2][1]*inertia[0][2])
        +inertia[2][0]*(inertia[1][2]*inertia[0][1]-inertia[1][1]*inertia[0][2]);
    det=1.0/det;
    inv_inr[0][0]= det*(inertia[2][2]*inertia[1][1]-inertia[2][1]*inertia[1][2]);
    inv_inr[1][1]= det*(inertia[2][2]*inertia[0][0]-inertia[2][0]*inertia[0][2]);
    inv_inr[2][2]= det*(inertia[1][1]*inertia[0][0]-inertia[1][0]*inertia[0][1]);
    inv_inr[0][1]=-det*(inertia[2][2]*inertia[1][0]-inertia[2][0]*inertia[1][2]);
    inv_inr[0][2]= det*(inertia[2][1]*inertia[1][0]-inertia[2][0]*inertia[1][1]);
    inv_inr[1][2]=-det*(inertia[2][1]*inertia[0][0]-inertia[2][0]*inertia[0][1]);
    inv_inr[1][0]= inv_inr[0][1];
    inv_inr[2][0]= inv_inr[0][2];
    inv_inr[2][1]= inv_inr[1][2];

    if(0) {
      printf("\nInverse of Inertia Tensor\n");
      for(i=0;i<3;i++) {
          for(k=0;k<3;k++) printf("%16.8lf",inv_inr[i][k]);
              printf("\n");
      }
    }
    //calculate angular momentum
    for(i=0;i<natom;i++) { // L = sum_i mi (ri x vi) 
        m=atom[i]->mass;
        angmom[0] += m*(relpv[i][1]*relvel[i][2]-relpv[i][2]*relvel[i][1]);
        angmom[1] += m*(relpv[i][2]*relvel[i][0]-relpv[i][0]*relvel[i][2]);
        angmom[2] += m*(relpv[i][0]*relvel[i][1]-relpv[i][1]*relvel[i][0]);
    }

    //calculate angular velocity
    for(i=0;i<3;i++) { // w = I^-1 L
        for(k=0;k<3;k++) omega[i]+=inv_inr[i][k]*angmom[k];
    }

    //calculate velocity due to rotation
    for(i=0;i<natom;i++) { //vr = w x r
        atom[i]->vr[0]=(omega[1]*relpv[i][2]-omega[2]*relpv[i][1]);
        atom[i]->vr[1]=(omega[2]*relpv[i][0]-omega[0]*relpv[i][2]);
        atom[i]->vr[2]=(omega[0]*relpv[i][1]-omega[1]*relpv[i][0]);
    }

    if(0) {
      printf("angular momentum\n");
      printf("w %8.4f %8.4f %8.4f\n",angmom[0],angmom[1],angmom[2]);
      printf("angular velocities\n");
      printf("w %8.4f %8.4f %8.4f\n",omega[0],omega[1],omega[2]);
      printf("rotational velocities\n");
      for(i=0;i<natom;i++) printf("%d %8.4f %8.4f %8.4f\n",i+1,atom[i]->vr[0],atom[i]->vr[1],atom[i]->vr[2]);
    }

    //calculate vibrational velocities
    for(i=0;i<natom;i++) {
       for(k=0; k<3; k++) atom[i]->vv[k]= relvel[i][k]-(atom[i]->vr[k]);
    }

    if(0) {
      printf("vibration velocities\n");
      for(i=0;i<natom;i++) printf("%d %8.4f %8.4f %8.4f\n",i+1,atom[i]->vv[0],atom[i]->vv[1],atom[i]->vv[2]);
    }

    //calculate principle moments of inertia
    double evl[3],**evt,pomega[3];
    evt=new double *[3]; for(i=0;i<3;i++) evt[i]=new double [3];
    for(i=0;i<3;i++) evl[i]=evt[i][0]=evt[i][1]=evt[i][2]=pomega[i]=0;
    jacobi(inertia,3,evl,evt); //inertia tensor,dimension,eigenvalues,eigenvectors
    eigsrt(evl,evt,3); //sort eigenvalues
    for(i=0;i<3;i++) pI[i]=evl[i]/(Na*1e23); //principle moment of inertia
    //for(i=0;i<3;i++) rotT[i]=h*h*Na*1e23/(8.0*PI*PI*evl[i]*kb); rotational temperatures

    if(0) {
      printf("\nPrinciple moments\n");
      for(i=0;i<3;i++) printf("%16.8lf vector %16.8lf %16.8lf %16.8lf\n",evl[i],evt[0][i],evt[1][i],evt[2][i]);
      printf("\n");
      printf("\nRotational temperatures\n");
      for(i=0;i<3;i++) printf("I %16.8lf (g/mol*A2) theta %16.8lf K\n",evl[i],h*h/(8.0*PI*PI*pI[i]*kb));
      printf("\n");
    }

    //calculate angular velocities along principle axises
    for(i=0;i<3;i++) pomega[i]=0;
    for(i=0;i<3;i++) for(k=0;k<3;k++) pomega[i]+=angmom[k]*evt[k][i]/evl[i];

    //calculate inertia weighted angular velocities
    for(i=0;i<3;i++) {
        if(1) // weighted by principle moments of inertia 
           anguv[i]= pomega[0]*evt[i][0]*sqrt(evl[0])
                   + pomega[1]*evt[i][1]*sqrt(evl[1])
                   + pomega[2]*evt[i][2]*sqrt(evl[2]);
        else
           anguv[i]= pomega[0]*evt[i][0]
                   + pomega[1]*evt[i][1]
                   + pomega[2]*evt[i][2];
    }

    if(0) {
      printf("angular velocities along principle axises\n");
      printf("w %8.4f %8.4f %8.4f\n",pomega[0],pomega[1],pomega[2]);
      for(i=0;i<3;i++) printf("%8.4f v %8.4f %8.4f %8.4f\n",pomega[i],pomega[i]*evt[0][i],pomega[i]*evt[1][i],pomega[i]*evt[2][i]);
    }

    if(0) {
      printf("translational velocities\n");
      for(i=0;i<natom;i++) printf("%d %8.4f %8.4f %8.4f\n",i+1,atom[i]->vt[0],atom[i]->vt[1],atom[i]->vt[2]);
      printf("rotational velocities\n");
      for(i=0;i<natom;i++) printf("%d %8.4f %8.4f %8.4f\n",i+1,atom[i]->vr[0],atom[i]->vr[1],atom[i]->vr[2]);
      printf("vibrational velocities\n");
      for(i=0;i<natom;i++) printf("%d %8.4f %8.4f %8.4f\n",i+1,atom[i]->vv[0],atom[i]->vv[1],atom[i]->vv[2]);
    }
    for(i=0;i<3;i++) delete [] evt[i];
    delete [] evt;

    return 0;
}

int MOLECULE::cal_Zmat(double (*relpv)[3],double **Zmat)
{
    int i,a2,a3,a4;
    Zmat[0][0]=Zmat[0][1]=Zmat[0][2]=0;
    i=1; a2=iccnt[i][0];
    Zmat[i][0]=dist(relpv[i],relpv[a2]);  Zmat[i][1]=Zmat[i][2]=0;
    i=2; a2=iccnt[i][0]; a3=iccnt[i][1]; 
    Zmat[i][0]=dist(relpv[i],relpv[a2]);  Zmat[i][1]=angl(relpv[i],relpv[a2],relpv[a3]); Zmat[i][2]=0;
    for(i=3;i<natom;i++) {
        a2=iccnt[i][0]; a3=iccnt[i][1];  a4=iccnt[i][2];
        Zmat[i][0]=dist(relpv[i],relpv[a2]);
        Zmat[i][1]=angl(relpv[i],relpv[a2],relpv[a3]);
        Zmat[i][2]=dihe(relpv[i],relpv[a2],relpv[a3],relpv[a4]);
    }
    return 0;
}


int MOLECULE::ini_iccnt()
{
    if(natom==0) return 1;
    if(iccnt!=NULL) delete [] iccnt;
    iccnt=new int [natom][5];

    if(mqvib=NULL) delete [] mqvib;
    mqvib=new double [3*natom];

    if(moi=NULL) delete [] moi;
    moi=new double [3*natom];

    return 0;
}


int GROUP::init_atom()
{
   if(natom==0) return 1;
   if(atm!=NULL) delete [] atm;
   atm=new int [natom];
   return 0;
}

int CELL::s2r2H(double *s2r)
{
  H[0][0]=s2r[0];
  H[1][1]=s2r[1];
  H[2][2]=s2r[2];
  H[1][0]=s2r[5];
  H[2][0]=s2r[4];
  H[2][1]=s2r[3];
  H[0][1]=0;
  H[0][2]=0;
  H[1][2]=0;

  //printf("H    %8.2f %8.2f %8.2f\n",H[0][0],H[0][1],H[0][2]);
  //printf("H    %8.2f %8.2f %8.2f\n",H[1][0],H[1][1],H[1][2]);
  //printf("H    %8.2f %8.2f %8.2f\n",H[2][0],H[2][1],H[2][2]);
}

int CELL::cal_Hinv()
{
  //inverse of H matrix
  Hinv[0][0]=1/H[0][0];
  Hinv[1][1]=1/H[1][1];
  Hinv[2][2]=1/H[2][2];
  Hinv[1][0]=-H[1][0]/H[0][0]/H[1][1];
  Hinv[2][0]=(H[1][0]*H[2][1]-H[1][1]*H[2][0])/H[0][0]/H[1][1]/H[2][2];
  Hinv[2][1]=-H[2][1]/H[1][1]/H[2][2];
  Hinv[0][1]=0;
  Hinv[1][2]=0;
  Hinv[0][2]=0;

  //printf("Hinv %8.5f %8.5f %8.5f\n",Hinv[0][0],Hinv[0][1],Hinv[0][2]);
  //printf("Hinv %8.5f %8.5f %8.5f\n",Hinv[1][0],Hinv[1][1],Hinv[1][2]);
  //printf("Hinv %8.5f %8.5f %8.5f\n",Hinv[2][0],Hinv[2][1],Hinv[2][2]);
}

int CELL::H2others()
{
  //  the H matrix
  //  h00  h01 h02     a1 b1 c1   h00   0   0
  //  h10  h11 h12     a2 b2 c2   h10 h11   0
  //  h20  h21 h22     a3 b3 c3   h20 h21 h22

   double fac=180.0/acos(-1.0);
// a=(H[0][0],H[1][0],H[2][0])=(s2r[0],s2r[5],s2r[4])
// b=(0      ,H[1][1],H[2][1])=(0     ,s2r[1],s2r[3])
// c=(0      ,0      ,H[2][2])=(0     ,0     ,s2r[2])
   a[0]=H[0][0]; a[1]=H[1][0]; a[2]=H[2][0];
   b[0]=0;       b[1]=H[1][1]; b[2]=H[2][1];
   c[0]=0;       c[1]=0;       c[2]=H[2][2];
   la = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
   lb = sqrt(b[1]*b[1] + b[2]*b[2]);
   lc = c[2];
   alpha = fac * acos( b[2] / lb );
   beta =  fac * acos( a[2] / la );
   gamma = fac * acos( (a[1]*b[1] + a[2]*b[2])/la/lb );

   //center of box
   cb[0]=0.5*(a[0]+b[0]+c[0]);
   cb[1]=0.5*(a[1]+b[1]+c[1]);
   cb[2]=0.5*(a[2]+b[2]+c[2]);

   //inverse of H matrix
   cal_Hinv();

   cal_volume();
   return 0;
}

int CELL::labc2H()
{
   //given la,lb,lc,alpha,beta,gamma, get H2
   double aa,bb,cc;
   double fac=acos(-1.0)/180.0;
   aa=alpha*fac; bb=beta*fac; cc=gamma*fac;
   //place c along Z
   c[0]=0; c[1]=0; c[2]=lc;
   //place b in the yz plane
   b[0]=0; b[1]=lb*sin(aa); b[2]=lb*cos(aa);
   //find vector for a
   a[2]=la*cos(bb);
   a[1]=(la*lb*cos(cc)-a[2]*b[2])/b[1];  //la*(cos(cc)-cos(aa)*cos(bb))/sin(aa);
   a[0]=sqrt(la*la-a[1]*a[1]-a[2]*a[2]);

   H[0][0]=a[0]; H[1][0]=a[1]; H[2][0]=a[2];
   H[0][1]=b[0]; H[1][1]=b[1]; H[2][1]=b[2];
   H[0][2]=c[0]; H[1][2]=c[1]; H[2][2]=c[2];

   return 0;
}

int CELL::cal_volume()
{
   volume=0;
   volume=fabs(H[0][0]*H[1][1]*H[2][2]);
   return 0;
}

int CELL::Map2UnitCell(double *ori,double *shift)
{ //maps ori+shift to within the unit cell

  // cm[3x1] = H[3x3] n[3x1]
  //cal_Hinv(); //inverse of H matrix

  double n[3];
  n[0]=Hinv[0][0]*ori[0]+Hinv[0][1]*ori[1]+Hinv[0][2]*ori[2];
  n[1]=Hinv[1][0]*ori[0]+Hinv[1][1]*ori[1]+Hinv[1][2]*ori[2];
  n[2]=Hinv[2][0]*ori[0]+Hinv[2][1]*ori[1]+Hinv[2][2]*ori[2];

  //cout<<"init  n "<<n[0]<<" "<<n[1]<<" "<<n[2]<<endl;
  while(n[0]>1) n[0]-=1;
  while(n[0]<0) n[0]+=1;
  while(n[1]>1) n[1]-=1;
  while(n[1]<0) n[1]+=1;
  while(n[2]>1) n[2]-=1;
  while(n[2]<0) n[2]+=1;

  //cout<<"final n "<<n[0]<<" "<<n[1]<<" "<<n[2]<<endl;
  shift[0]=H[0][0]*n[0]+H[0][1]*n[1]+H[0][2]*n[2]-ori[0];
  shift[1]=H[1][0]*n[0]+H[1][1]*n[1]+H[1][2]*n[2]-ori[1];
  shift[2]=H[2][0]*n[0]+H[2][1]*n[1]+H[2][2]*n[2]-ori[2];

  return 0;
}

int CELL::mindst(double *in,double *pv)
{ //find the shortest vector of in

  // cm[3x1] = H[3x3] n[3x1]
  //cal_Hinv(); //inverse of H matrix

  double n[3],ori[3];
  ori[0]=in[0]+cb[0]; //move to the center of box
  ori[1]=in[1]+cb[1]; //move to the center of box
  ori[2]=in[2]+cb[2]; //move to the center of box
  n[0]=Hinv[0][0]*ori[0]+Hinv[0][1]*ori[1]+Hinv[0][2]*ori[2];
  n[1]=Hinv[1][0]*ori[0]+Hinv[1][1]*ori[1]+Hinv[1][2]*ori[2];
  n[2]=Hinv[2][0]*ori[0]+Hinv[2][1]*ori[1]+Hinv[2][2]*ori[2];

  //cout<<"init  n "<<n[0]<<" "<<n[1]<<" "<<n[2]<<endl;
  while(n[0]>1) n[0]-=1;
  while(n[0]<0) n[0]+=1;
  while(n[1]>1) n[1]-=1;
  while(n[1]<0) n[1]+=1;
  while(n[2]>1) n[2]-=1;
  while(n[2]<0) n[2]+=1;

  //cout<<"final n "<<n[0]<<" "<<n[1]<<" "<<n[2]<<endl;
  pv[0]=H[0][0]*n[0]+H[0][1]*n[1]+H[0][2]*n[2]-cb[0];
  pv[1]=H[1][0]*n[0]+H[1][1]*n[1]+H[1][2]*n[2]-cb[1];
  pv[2]=H[2][0]*n[0]+H[2][1]*n[1]+H[2][2]*n[2]-cb[2];

  return 0;
}

int CELL::out_cell(ostream *outf)
{
   char null[1024];
   *outf<<"The H matrix"<<endl;
   sprintf(null,"%10.4f %10.4f %10.4f",H[0][0],H[0][1],H[0][2]);
   *outf<<null<<endl;
   sprintf(null,"%10.4f %10.4f %10.4f",H[1][0],H[1][1],H[1][2]);
   *outf<<null<<endl;
   sprintf(null,"%10.4f %10.4f %10.4f",H[2][0],H[2][1],H[2][2]);
   *outf<<null<<endl;
   *outf<<"origin"<<endl;
   sprintf(null,"o: (%10.4f,%10.4f,%10.4f) ",o[0],o[1],o[2]);
   *outf<<null<<endl;
   *outf<<"The a,b,c vectors"<<endl;
   sprintf(null,"a: (%10.4f,%10.4f,%10.4f) |a|=%10.4f",a[0],a[1],a[2],la);
   *outf<<null<<endl;
   sprintf(null,"b: (%10.4f,%10.4f,%10.4f) |b|=%10.4f",b[0],b[1],b[2],lb);
   *outf<<null<<endl;
   sprintf(null,"c: (%10.4f,%10.4f,%10.4f) |c|=%10.4f",c[0],c[1],c[2],lc);
   *outf<<null<<endl;
   sprintf(null,"alpha %10.4f beta %10.4f gamma %10.4f",alpha,beta,gamma);
   *outf<<null<<endl;
   cal_volume();
   sprintf(null,"volume %10.4f A3",volume);
   *outf<<null<<endl;

   return 0;
}

int CELL::set_grid(double insgrid)
{

    int tg[3];
    tg[0]=(int)(la/insgrid)+1;
    tg[1]=(int)(lb/insgrid)+1;
    tg[2]=(int)(lc/insgrid)+1;
    if(insgrid==sgrid && tg[0]==ngrid[0] &&
        tg[1]==ngrid[1] && tg[2]==ngrid[2] )  return 0; //grids have already been set
    sgrid=insgrid;
    int i,j,k;
    free_grid();

    for(i=0;i<3;i++) ngrid[i]=tg[i];
    vgrid=new char **[ ngrid[0] ];
    for(i=0;i<ngrid[0];i++) {
       vgrid[i]=new char *[ ngrid[1] ];
       for(j=0;j<ngrid[1];j++) {
          vgrid[i][j]=new char [ ngrid[2] ];
          for(k=0;k<ngrid[2];k++) vgrid[i][j][k]=0;
       }
    }
    //for(i=0;i<3;i++) printf("%d ngrid %d\n",i+1,ngrid[i]);
    return 0;
}

int CELL::free_grid()
{
    if(vgrid!=NULL) {
       int i,j;
       for(i=0;i<ngrid[0];i++) {
          for(j=0;j<ngrid[1];j++) delete [] vgrid[i][j];
          delete [] vgrid[i];
       }
       delete [] vgrid; vgrid=NULL;
    }
    return 0;
}

int BOND::cal_len()
{
    len=0;
    double v[3]; //v = R(C1-C2)
    int k;
    double ip=0;
    for(k=0;k<3;k++) {
       v[k]=(atom[0]->pv[k])-(atom[1]->pv[k]);
       ip+=(v[k]*v[k]);
    }
    len=sqrt(ip);
    return 0;
}

int ANGLE::cal_theta()
{
    theta=0;

    // C1-C2-C3
    // v[0] = R(C1-C2)  
    // v[1] = R(C3-C2) 
    double v[2][3];  //vectors of the three bonds
    int k;
    double ip00,ip11,ip01;

    ip00=ip11=ip01=0.0;

    for(k=0;k<3;k++) {
      v[0][k]=(atom[0]->pv[k])-(atom[1]->pv[k]);
      v[1][k]=(atom[2]->pv[k])-(atom[1]->pv[k]);
      ip00+=(v[0][k]*v[0][k]);
      ip11+=(v[1][k]*v[1][k]);
      ip01+=(v[0][k]*v[1][k]);
    }

    costh=ip01/sqrt(ip00*ip11);
    theta=acos(costh);  // 0<theta<pi
    return 0;
}

int TORSION::cal_theta()
{
    theta=0;

    // C1-C2-C3-C4
    // v[0] = R(C2-C1) F
    // v[1] = R(C3-C2) G
    // v[2] = R(C4-C3) H
    // n    = v[1]xv[2] = GxH
    double v[3][3];  //vectors of the three bonds
    double n[3];
    int i,j,k;
    double ip00,ip11,ip22,ip01,ip02,ip12,ip01x2;
    double x,y;  //x: component of v[0] on projected (perpendicular to v[1]) v[2]
                 //y: component of v[0] on the normal vector of C2-C3-C4 plane

    ip00=ip11=ip22=ip01=ip02=ip12=ip01x2=0.0;

    for(k=0;k<3;k++) {
      v[0][k]=(atom[1]->pv[k])-(atom[0]->pv[k]);
      v[1][k]=(atom[2]->pv[k])-(atom[1]->pv[k]);
      v[2][k]=(atom[3]->pv[k])-(atom[2]->pv[k]);
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

    //x=-v[0]*v[2]+(v[0]*v[1])(v[1]*v[2])= v[0]*(-v[2]+(v[1]*v[2])v[1])
    x=-ip02/sqrt(ip00*ip22)+(ip01*ip12)/(ip11*sqrt(ip00*ip22));
    //y=v[0]*(v[1]xv[2])
    y=-ip01x2/sqrt(ip00*ip11*ip22);
    //printf("%f %f %f %f %f %f %f \n",ip00,ip11,ip22,ip01,ip02,ip12,ip01x2);
    //printf("%e %e \n",x,y);
    
    costh=x/sqrt(x*x+y*y);
    theta=acos(costh);
    if( y<0 ) theta*=-1;

    return 0;
}

int TORSION::cal_theta1()
{
    theta=0;

    // C1-C2-C3-C4
    // v[0] = R(C2-C1) F
    // v[1] = R(C3-C2) G
    // v[2] = R(C4-C3) H
    // n[0] = v[0]xv[1] = FxG
    // n[1] = v[1]xv[2] = GxH
    double v[3][3];  //vectors of the three bonds
    double n[2][3];
    int i,j,k;
    double n0,n1,ip;

    n0=n1=ip=0;

    for(k=0;k<3;k++) {
      v[0][k]=(atom[1]->pv[k])-(atom[0]->pv[k]);
      v[1][k]=(atom[2]->pv[k])-(atom[1]->pv[k]);
      v[2][k]=(atom[3]->pv[k])-(atom[2]->pv[k]);
    }

    i=0; j=1; k=0;
    n[k][0]=v[i][1]*v[j][2]-v[i][2]*v[j][1];
    n[k][1]=v[i][2]*v[j][0]-v[i][0]*v[j][2];
    n[k][2]=v[i][0]*v[j][1]-v[i][1]*v[j][0];

    i=1; j=2; k=1;
    n[k][0]=v[i][1]*v[j][2]-v[i][2]*v[j][1];
    n[k][1]=v[i][2]*v[j][0]-v[i][0]*v[j][2];
    n[k][2]=v[i][0]*v[j][1]-v[i][1]*v[j][0];

    for(k=0;k<3;k++) { n0+=n[0][k]*n[0][k]; n1+=n[1][k]*n[1][k]; ip+=n[0][k]*n[1][k]; }
    printf("%d %d %d %d\n",atm[0],atm[1],atm[2],atm[3]);
    for(k=0;k<3;k++) printf("v%d %f %f %f\n",k,v[k][0],v[k][1],v[k][2]);
    printf("n0 %f n1 %f ip %f\n",n0,n1,ip);

    costh=ip/sqrt(n0*n1);
    theta=acos(costh);

    return 0;
}

int INVERSION::cal_theta()
{
    theta=0; // angle between C4 and the plan of C2-C1-C3
             // in Cerius2 inverstion selection C1 C2 C3 C4 or C1 C3 C2 C4
    // C2-C1-C3
    //     |
    //    C4
    // v[0] = R(C4-C1)    r[0]=length of v[0]
    // v[1] = R(C3-C1)    r[1]=length of v[1]
    // v[2] = R(C2-C1)    r[2]=length of v[2]
    // n    = v[1]xv[2]
    double v[3][3];  //vectors of the three bonds
    double n[3];
    int i,j,k;
    double ip00,ip01x2,ipnn;
    double sinth; //sin(theta) 

    ip00=ip01x2=ipnn=0.0;

    for(k=0;k<3;k++) {
      v[0][k]=(atom[3]->pv[k])-(atom[0]->pv[k]);
      v[1][k]=(atom[2]->pv[k])-(atom[0]->pv[k]);
      v[2][k]=(atom[1]->pv[k])-(atom[0]->pv[k]);
      ip00+=(v[0][k]*v[0][k]);
    }

    i=1; j=2;
    n[0]=v[i][1]*v[j][2]-v[i][2]*v[j][1];
    n[1]=v[i][2]*v[j][0]-v[i][0]*v[j][2];
    n[2]=v[i][0]*v[j][1]-v[i][1]*v[j][0];

    for(k=0;k<3;k++) { ip01x2+=(v[0][k]*n[k]); ipnn+=(n[k]*n[k]); }

    //sinth=v[0]*(v[1]xv[2])
    sinth=ip01x2/sqrt(ip00*ipnn);
    //printf("n  %f %f %f\n",n[0]/sqrt(ipnn),n[1]/sqrt(ipnn),n[2]/sqrt(ipnn));
    //printf("v0 %f %f %f\n",v[0][0]/sqrt(ip00),v[0][1]/sqrt(ip00),v[0][2]/sqrt(ip00));
    //printf("%f %f %f %f %f %f %f \n",ip00,ip11,ip22,ip01,ip02,ip12,ip01x2);
    //printf("%f \n",sinth);
    
    costh=sqrt(1-sinth*sinth);
    theta=acos(costh);
    if( sinth<0 ) theta*=-1;

    return 0;
}

int RING::init_atm()
{
    if(size==0) return 1;
    if(atm!=NULL) delete [] atm;
    atm=new int [size];
    //if(sbond!=NULL) delete [] sbond;
    //sbond=new SBOND [size];

    return 0;
}


