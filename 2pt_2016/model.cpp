#include "model.h"

int MODEL::rd_strt()
{
   strcpy(name,ctl.in_strt);  //set the name of the model to the input filename
   switch(ctl.in_strt_flag) {
      case none:
        break;
      case strtgro:
        rd_gro(name);
        element2prp(); //assign atom mass,color,radius,crad based on element type
        dst2bond();   //from atom separation distance to bond connectivity
        bond2mol();   //find molecules
        bond2cnt();   //find connectivity
        cnt2bod();    //find bond order
        cnt2valence();//find bond, angle, torsion, inversion from connectivity
        farneigh();   //find the farthermost connected atom
        headtail();   //find heand and tail atoms
        break;
      case strtg96:
        rd_g96(name);
        timing.gettimetotal(&cout);
        element2prp(); //assign atom mass,color,radius,crad based on element type
        timing.gettimetotal(&cout);
        cout<<"element2prp done"<<endl;
        //dst2bond();   //from atom separation distance to bond connectivity
        mol2bond();   //from atom separation distance within a molecule to bond connectivity
        timing.gettimetotal(&cout);
        cout<<"mol2bond done"<<endl;
        bond2mol();   //find molecules
        timing.gettimetotal(&cout);
        cout<<"bond2mol done"<<endl;
        bond2cnt();   //find connectivity
        timing.gettimetotal(&cout);
        cnt2bod();    //find bond order
        timing.gettimetotal(&cout);
        cnt2valence();//find bond, angle, torsion, inversion from connectivity
        timing.gettimetotal(&cout);
        farneigh();   //find the farthermost connected atom
        headtail();   //find heand and tail atoms
        timing.gettimetotal(&cout);
        break;
      case strtlmp:
        rd_lmpdata(name);
        bond2mol();   //find molecules
        bond2cnt();   //find connectivity
        cnt2bod();    //find bond order
        farneigh();   //find the farthermost connected atom
        headtail();   //find heand and tail atoms
        mass2types(); //find element names and fftypes
        element2prp(); //assign atom mass,color,radius,crad based on element type
        break;
      case strtcpmd:
        rd_cpmdin(name);
        element2prp(); //assign atom mass,color,radius,crad based on element type
        dst2bond();
        bond2mol();   //find molecules
        bond2cnt();   //find connectivity
        cnt2bod();    //find bond order
        cnt2valence();//find bond, angle, torsion, inversion from connectivity
        farneigh();   //find the farthermost connected atom
        headtail();   //find heand and tail atoms
        break;
   }
   //rd_grp("none");
   cal_mass(); //calculate total mass of the model
   //cal_multipole();  //calculate multipole moments of the model
   ck_range(); //check atom range for nonperiodic systems
   cell.out_cell(&cout);
   return 0;
}

int MODEL::rd_lmpdata(char *indata)
{
    cout<<"Reading lammps data file "<<indata<<endl;
    ifstream inf(indata,fstream::in);

    if(!inf.is_open()) {
      cout<<" Error: Lammps data file "<<indata<<" cannot be opened."<<endl;
      return filenotopen;
    }
    
    char null[1024],*charptr;
    int i,j;  
    char delims[] = {" ;,\t"};

    inf.getline(null,1024); inf>>ws>>ws;
    inf.getline(null,1024); 
    while(null[0]>0) {
        charptr = strtok(null, delims);
        j=atoi(charptr);
        charptr = strtok(NULL, delims);
        if(strcmp(charptr,"atoms")==0) { natom=j; init_atom(); }
        else if(strcmp(charptr,"bonds")==0) { nbond=j; init_bond(); }
        else if(strcmp(charptr,"angles")==0) { nangle=j; init_angle(); }
        else if(strcmp(charptr,"dihedrals")==0) { ntor=j; init_tor(); }
        else if(strcmp(charptr,"impropers")==0) { nimp=j; init_imp(); }
        else if(strcmp(charptr,"inversions")==0) { ninv=j; init_inv(); }
        inf.getline(null,1024); //inf>>ws; stlin 2011.09.05
    }

    printf("stlin natom %d nbond %d nangle %d ntorsion %d nimproper %d ninversion %d\n",natom,nbond,nangle,ntor,nimp,ninv);

    inf.getline(null,1024);
    while(null[0]>0) {
        charptr = strtok(null, delims);
        j=atoi(charptr);
        charptr = strtok(NULL, delims);
        if(strcmp(charptr,"atom")==0) ff.nvdw=j;
        else if(strcmp(charptr,"off")==0) ff.noffvdw=j;
        else if(strcmp(charptr,"bond")==0) ff.nbond=j;
        else if(strcmp(charptr,"angle")==0) ff.nangle=j;
        else if(strcmp(charptr,"dihedral")==0) ff.ntor=j;
        else if(strcmp(charptr,"inversion")==0) ff.ninv=j;
        else if(strcmp(charptr,"improper")==0) ff.nimp=j;

        inf.getline(null,1024); 
    }
    
    ff.natmt=ff.nvdw;                 
    printf("stlin nvdw %d noff %d nbond %d nangle %d ntor %d ninv %d nimp %d\n",ff.nvdw,ff.noffvdw,ff.nbond,ff.nangle,ff.ntor,ff.ninv,ff.nimp);

    ff.init_atmt();
    ff.init_ff(&ff.vdw,ff.nvdw,1); 
    ff.init_ff(&ff.offvdw,ff.noffvdw,2);
    ff.init_ff(&ff.bond,ff.nbond,2);
    ff.init_ff(&ff.angle,ff.nangle,3);
    ff.init_ff(&ff.tor,ff.ntor,4);
    ff.init_ff(&ff.inv,ff.ninv,4);
    ff.init_ff(&ff.imp,ff.nimp,4);

    while(null[0]>0) {inf.getline(null,1024);}; 

    //box information
    periodic=3;
    inf>>cell.o[0]>>cell.a[0];  cell.a[1]=cell.a[2]=0.0; cell.a[0]-=cell.o[0]; cell.la=cell.a[0];
    inf.getline(null,1024);
    inf>>cell.o[1]>>cell.b[1];  cell.b[0]=cell.b[2]=0.0; cell.b[1]-=cell.o[1]; cell.lb=cell.b[1];
    inf.getline(null,1024);
    inf>>cell.o[2]>>cell.c[2];  cell.c[0]=cell.c[1]=0.0; cell.c[2]-=cell.o[2]; cell.lc=cell.c[2];
    cell.alpha=cell.beta=cell.gamma=90.0;
    cell.labc2H();
    cell.H2others();
    //cell.out_cell(&cout);

    do{inf>>null;}while(strcmp(null,"Masses")!=0); inf>>null;
    for(i=0;i<ff.natmt;i++) {
        inf>>ff.atmt[i].mass>>null;
    }

    //while(strcmp(null,"Nonbond")!=0 ) {inf>>null;}; inf>>null>>null;
    while(strcmp(null,"Nonbond")!=0 && strcmp(null,"Pair")!=0) {inf>>null;}; inf>>null>>null;
    for(i=0;i<ff.nvdw;i++) {
        inf.getline(null,1024);
        charptr = strtok(null, delims);
        j=0;
        while(charptr!=NULL) {
           ff.vdw[i].param[j++]=atof(charptr);
           //cout<<j<<" "<<charptr<<" "<<ff.vdw[i].param[j-1]<<endl;
           charptr = strtok(NULL, delims);
        }
        SWAP(&ff.vdw[i].param[0],&ff.vdw[i].param[1]); //Ro:param[0], Do:param[1]
        ff.vdw[i].nparam=ff.vdw[i].nmax=j;
        ff.vdw[i].fatmtid[0]=i;
        inf>>null;
    }

    if(ff.noffvdw) { while(strcmp(null,"Off")!=0) {inf>>null;}; inf>>null>>null>>null; }
    for(i=0;i<ff.noffvdw;i++) {
        inf>>null;
        ff.offvdw[i].fatmtid[0]=atoi(null)-1;
        inf>>null;
        ff.offvdw[i].fatmtid[1]=atoi(null)-1;
        inf.getline(null,1024);
        charptr = strtok(null, delims);
        j=0;
        while(charptr!=NULL) {
           ff.offvdw[i].param[j++]=atof(charptr);
           //cout<<j<<" "<<charptr<<" "<<ff.offvdw[i].param[j-1]<<endl;
           charptr = strtok(NULL, delims);
        }
        SWAP(&ff.offvdw[i].param[0],&ff.offvdw[i].param[1]); 
        ff.offvdw[i].nparam=ff.offvdw[i].nmax=j;
        inf>>null;
    }
    ff.init_vdwij();

    if(ff.nbond) { while(strcmp(null,"Bond")!=0) {inf>>null;}; inf>>null>>null; }
    for(i=0;i<ff.nbond;i++) {
        inf.getline(null,1024);
        charptr = strtok(null, delims);
        j=0;
        while(charptr!=NULL) {
           ff.bond[i].param[j++]=atof(charptr);
           //cout<<i<<" "<<j<<" "<<charptr<<" "<<ff.bond[i].param[j-1]<<endl;
           charptr = strtok(NULL, delims);
        }
        ff.bond[i].nparam=j;
        inf>>null;
    }

    if(ff.nangle) { while(strcmp(null,"Angle")!=0) {inf>>null;}; inf>>null>>null; }
    for(i=0;i<ff.nangle;i++) {
        inf.getline(null,1024);
        charptr = strtok(null, delims);
        j=0;
        while(charptr!=NULL) {
           ff.angle[i].param[j++]=atof(charptr);
           //cout<<i<<" "<<j<<" "<<charptr<<" "<<ff.angle[i].param[j-1]<<endl;
           charptr = strtok(NULL, delims);
        }
        ff.angle[i].nparam=j;
        inf>>null;
    }

    if(ff.ntor) { while(strcmp(null,"Dihedral")!=0) {inf>>null;}; inf>>null>>null; }
    for(i=0;i<ff.ntor;i++) {
        inf.getline(null,1024);
        charptr = strtok(null, delims);
        j=0;
        while(charptr!=NULL) {
           ff.tor[i].param[j++]=atof(charptr);
           //cout<<i<<" "<<j<<" "<<charptr<<" "<<ff.tor[i].param[j-1]<<endl;
           charptr = strtok(NULL, delims);
        }
        ff.tor[i].nparam=j;
        inf>>null;
    }

    if(ff.ninv) { while(strcmp(null,"Inversion")!=0) {inf>>null;}; inf>>null>>null; }
    for(i=0;i<ff.ninv;i++) {
        inf.getline(null,1024);
        charptr = strtok(null, delims);
        j=0;
        while(charptr!=NULL) {
           ff.inv[i].param[j++]=atof(charptr);
           //cout<<i<<" "<<j<<" "<<charptr<<" "<<ff.inv[i].param[j-1]<<endl;
           charptr = strtok(NULL, delims);
        }
        ff.inv[i].nparam=j;
        inf>>null;
    }

    while(strcmp(null,"Atoms")!=0) {inf>>null;};
    for(i=0;i<natom;i++) {
      inf>>atom[i].id>>atom[i].mol>>atom[i].ffid>>atom[i].chg>>atom[i].pv[0]>>atom[i].pv[1]>>atom[i].pv[2];
      atom[i].id--;
      atom[i].ffid--;
      atom[i].mol=-1;
      atom[i].mass=ff.atmt[atom[i].ffid].mass;
      //printf("stlin id %d chg %f mass %f pv %f %f %f\n",atom[i].id,atom[i].chg,atom[i].mass,atom[i].pv[0],atom[i].pv[1],atom[i].pv[2]);
    }

    inf>>null;
    if(strcmp(null,"Velocities")==0) {
      for(i=0;i<natom;i++) {
         inf>>null>>atom[i].vel[0]>>atom[i].vel[1]>>atom[i].vel[2];
         atom[i].vel[0] *= 1000.0;
         atom[i].vel[1] *= 1000.0;
         atom[i].vel[2] *= 1000.0;
      }
    }

    if(nbond) {
       while(strcmp(null,"Bonds")!=0) inf>>null;
       for(i=0;i<nbond;i++) {
         inf>>bond[i].id>>bond[i].ffid>>bond[i].atm[0]>>bond[i].atm[1];
         bond[i].id--;
         bond[i].ffid--;
         bond[i].atm[0]--; bond[i].atm[1]--;
         bond[i].atom[0]=&atom[bond[i].atm[0]];
         bond[i].atom[1]=&atom[bond[i].atm[1]];
         //printf("lll %d \n",bond[i].atom[0]->id,bond[i].atom[1]->id);
         //printf("stlin id %d atoms %d %d %d %d\n",bond[i].id,bond[i].atm[0],bond[i].atom[0]->id,bond[i].atm[1],bond[i].atom[1]->id);
//here
       }
    }

    if(nangle) {
       do{inf>>null;}while(strcmp(null,"Angles")!=0);
       for(i=0;i<nangle;i++) {
         inf>>angle[i].id>>angle[i].ffid>>angle[i].atm[0]>>angle[i].atm[1]>>angle[i].atm[2];
         angle[i].id--;
         angle[i].ffid--;
         angle[i].atm[0]--; angle[i].atm[1]--; angle[i].atm[2]--;
         angle[i].atom[0]=&atom[angle[i].atm[0]];
         angle[i].atom[1]=&atom[angle[i].atm[1]];
         angle[i].atom[2]=&atom[angle[i].atm[2]];
         //printf("stlin id %d atoms %d %d %d \n",angle[i].id,angle[i].atm[0],angle[i].atm[1],angle[i].atm[2]);
       }
    }

    if(ntor) {
       do{inf>>null;}while(strcmp(null,"Dihedrals")!=0);
       for(i=0;i<ntor;i++) {
         inf>>tor[i].id>>tor[i].ffid>>tor[i].atm[0]>>tor[i].atm[1]>>tor[i].atm[2]>>tor[i].atm[3];
         inf.getline(null,1024); charptr = strtok(null, delims); if(charptr!=NULL) tor[i].deg=atoi(charptr);
         tor[i].id--;
         tor[i].ffid--;
         tor[i].atm[0]--; tor[i].atm[1]--; tor[i].atm[2]--; tor[i].atm[3]--;
         tor[i].atom[0]=&atom[tor[i].atm[0]];
         tor[i].atom[1]=&atom[tor[i].atm[1]];
         tor[i].atom[2]=&atom[tor[i].atm[2]];
         tor[i].atom[3]=&atom[tor[i].atm[3]];
         //printf("stlin id %d atoms %d %d %d %d\n",tor[i].id,tor[i].atm[0],tor[i].atm[1],tor[i].atm[2],tor[i].atm[3]);
       }
    }

    if(ninv) {
       do{inf>>null;}while(strcmp(null,"Inversions")!=0);
       for(i=0;i<ninv;i++) {
         inf>>inv[i].id>>inv[i].ffid>>inv[i].atm[0]>>inv[i].atm[1]>>inv[i].atm[2]>>inv[i].atm[3];
         inv[i].id--;
         inv[i].ffid--;
         inv[i].atm[0]--; inv[i].atm[1]--; inv[i].atm[2]--; inv[i].atm[3]--;
         inv[i].atom[0]=&atom[inv[i].atm[0]];
         inv[i].atom[1]=&atom[inv[i].atm[1]];
         inv[i].atom[2]=&atom[inv[i].atm[2]];
         inv[i].atom[3]=&atom[inv[i].atm[3]];
         //printf("stlin id %d atoms %d %d %d %d\n",inv[i].id,inv[i].atm[0],inv[i].atm[1],inv[i].atm[2],inv[i].atm[3]);
         //printf("stlin id %d atoms %d %d %d %d\n",inv[i].id,inv[i].atom[0]->id,inv[i].atom[1]->id,inv[i].atom[2]->id,inv[i].atom[3]->id);
       }
    }
    
    return 0;
}

int MODEL::out_lmpdata(ostream *outf)
{
    struct tm *newtime;
    time_t ltime;
    
    char null[1024];
    int i,j;  

    time(NULL);
    newtime = gmtime( &ltime );
    sprintf(null,"LAMMPS %4d-%2d-%2d data file for %s",newtime->tm_year+1900,
                  newtime->tm_mon+1,newtime->tm_mday,name);
    if (null[12]==32) null[12]=48; if (null[15]==32) null[15]=48;
    *outf<<null<<endl<<endl;

    sprintf(null,"%7d atoms",natom);
    *outf<<null<<endl;
    sprintf(null,"%7d bonds",nbond);
    *outf<<null<<endl;
    sprintf(null,"%7d angles",nangle);
    *outf<<null<<endl;
    sprintf(null,"%7d dihedrals",ntor);
    *outf<<null<<endl;
    sprintf(null,"%7d impropers",nimp);
    *outf<<null<<endl;
    if(ninv) {
    sprintf(null,"%7d inversions",ninv);
    *outf<<null<<endl; }
    *outf<<endl;
    sprintf(null,"%4d atom types",ff.nvdw);
    *outf<<null<<endl;
    if(ff.noffvdw) {
    sprintf(null,"%4d off diagonal types",ff.noffvdw);
    *outf<<null<<endl;
    }
    if(ff.nbond) {
        sprintf(null,"%4d bond types",ff.nbond);
        *outf<<null<<endl;
    }
    if(ff.nangle) {
        sprintf(null,"%4d angle types",ff.nangle);
        *outf<<null<<endl;
    }
    if(ff.ntor) {
        sprintf(null,"%4d dihedral types",ff.ntor);
        *outf<<null<<endl;
    }
    if(ff.ninv) {
        sprintf(null,"%4d inversion types",ff.ninv);
        *outf<<null<<endl;
    }
    if(ff.nimp) {
        sprintf(null,"%4d improper types",ff.nimp);
        *outf<<null<<endl;
    }
    *outf<<endl;

    //stlin 2012.2.18 modified for lammps triclinic cells
    double lx,ly,lz,xy,yz,xz,fac;
    fac = acos(-1.0)/180.0;
    lx = cell.la;
    xy = cell.lb*cos(cell.gamma*fac);
    xz = cell.lc*cos(cell.beta*fac);
    ly = sqrt( cell.lb*cell.lb - xy*xy);
    yz = (cell.lb*cell.lc*cos(cell.alpha*fac)-xy*xz)/ly;
    lz = sqrt( cell.lc*cell.lc - xz*xz - yz*yz);
    sprintf(null," %15.9f %15.9f xlo xhi",cell.o[0],lx);
    *outf<<null<<endl;
    sprintf(null," %15.9f %15.9f ylo yhi",cell.o[1],ly);
    *outf<<null<<endl;
    sprintf(null," %15.9f %15.9f zlo zhi",cell.o[2],lz);
    *outf<<null<<endl;
    sprintf(null," %15.9f %15.9f %15.9f xy xz yz",xy,xz,yz);
    *outf<<null<<endl<<endl;

    sprintf(null,"Masses");
    *outf<<null<<endl<<endl;
    for(i=0;i<ff.natmt;i++) {
      sprintf(null,"%4d %10f",i+1,ff.atmt[i].mass);
      *outf<<null<<endl;
    }
    *outf<<endl;

    sprintf(null,"Pair Coeffs");
    *outf<<null<<endl<<endl;
    for(i=0;i<ff.nvdw;i++) {
      sprintf(null,"%4d ",i+1);
      for(j=ff.vdw[i].nparam-1;j>-1;j--) {
         sprintf(null,"%s %13.10f ",null,ff.vdw[i].param[j]);
      }
      *outf<<null<<endl;
    }
    *outf<<endl;

    if(ff.noffvdw) {
       sprintf(null,"Off Diagonal Nonbond Coeffs"); *outf<<null<<endl<<endl;
       for(i=0;i<ff.noffvdw;i++) {
         sprintf(null,"%d %d ",ff.offvdw[i].fatmtid[0]+1,ff.offvdw[i].fatmtid[1]+1);
         for(j=0;j<ff.offvdw[i].nparam;j++) {
            sprintf(null,"%s%20.10f ",null,ff.offvdw[i].param[j]);
         }
         *outf<<null<<endl;
       }
       *outf<<endl;
    }

    if(ff.nbond) {
       sprintf(null,"Bond Coeffs"); *outf<<null<<endl<<endl;
       for(i=0;i<ff.nbond;i++) {
         sprintf(null,"%4d ",i+1);
         for(j=0;j<ff.bond[i].nparam;j++) {
            sprintf(null,"%s %13.10f ",null,ff.bond[i].param[j]);
         }
         *outf<<null<<endl;
       }
       *outf<<endl;
    }

    if(ff.nangle) {
       sprintf(null,"Angle Coeffs"); *outf<<null<<endl<<endl;
       for(i=0;i<ff.nangle;i++) {
         sprintf(null,"%4d ",i+1);
         for(j=0;j<ff.angle[i].nparam;j++) {
            sprintf(null,"%s %13.10f ",null,ff.angle[i].param[j]);
         }
         *outf<<null<<endl;
       }
       *outf<<endl;
    }

    if(ff.ntor) {
       sprintf(null,"Dihedral Coeffs"); *outf<<null<<endl<<endl;
       for(i=0;i<ff.ntor;i++) {
         sprintf(null,"%4d ",i+1);
         for(j=0;j<ff.tor[i].nparam;j++) {
            sprintf(null,"%s %13.10f ",null,ff.tor[i].param[j]);
         }
         *outf<<null<<endl;
       }
       *outf<<endl;
    }

    if(ff.ninv) {
       sprintf(null,"Inversion Coeffs"); *outf<<null<<endl<<endl;
       for(i=0;i<ff.ninv;i++) {
         sprintf(null,"%4d ",i+1);
         for(j=0;j<ff.inv[i].nparam;j++) {
            sprintf(null,"%s %13.10f ",null,ff.inv[i].param[j]);
         }
         *outf<<null<<endl;
       }
       *outf<<endl;
    }

    if(ff.nimp) {
       sprintf(null,"Improper Coeffs"); *outf<<null<<endl<<endl;
       for(i=0;i<ff.nimp;i++) {
         sprintf(null,"%4d ",i+1);
         for(j=0;j<ff.imp[i].nparam;j++) {
            sprintf(null,"%s %13.10f ",null,ff.imp[i].param[j]);
         }
         *outf<<null<<endl;
       }
       *outf<<endl;
    }

    sprintf(null,"Atoms");
    *outf<<null<<endl<<endl;
    for(i=0;i<natom;i++) {
       //sprintf(null,"%d %d %d %11f %15.9f %15.9f %15.9f",atom[i].id+1,0,atom[i].ffid+1,atom[i].chg,atom[i].pv[0]+atom[i].shift[0],atom[i].pv[1]+atom[i].shift[1],atom[i].pv[2]+atom[i].shift[2]);
       sprintf(null,"%7d %6d %3d %11f %15.9f %15.9f %15.9f",atom[i].id+1,atom[i].mol+1,atom[i].ffid+1,atom[i].chg,atom[i].pv[0],atom[i].pv[1],atom[i].pv[2]);
       *outf<<null<<endl;
   }
   *outf<<endl;

    sprintf(null,"Velocities");
    *outf<<null<<endl<<endl;
    for(i=0;i<natom;i++) {
       /* the unit of velocity in trajectory is A/ps; need to convert to A/fs */
       sprintf(null,"%d %15.9f %15.9f %15.9f",atom[i].id+1,atom[i].vel[0]*1E-3,atom[i].vel[1]*1E-3,atom[i].vel[2]*1E-3);
       *outf<<null<<endl;
   }
   *outf<<endl;

    if(nbond) {
       sprintf(null,"Bonds");
       *outf<<null<<endl<<endl;
       for(i=0;i<nbond;i++) {
          sprintf(null,"%6d %3d %6d %6d",bond[i].id+1,bond[i].ffid+1,bond[i].atm[0]+1,bond[i].atm[1]+1);
          *outf<<null<<endl;
      }
      *outf<<endl;
    }

    if(nangle) {
       sprintf(null,"Angles");
       *outf<<null<<endl<<endl;
       for(i=0;i<nangle;i++) {
          sprintf(null,"%6d %3d %6d %6d %6d",angle[i].id+1,angle[i].ffid+1,angle[i].atm[0]+1,angle[i].atm[1]+1,angle[i].atm[2]+1);
          *outf<<null<<endl;
      }
      *outf<<endl;
    }

    if(ntor) {
       sprintf(null,"Dihedrals");
       *outf<<null<<endl<<endl;
       for(i=0;i<ntor;i++) {
          sprintf(null,"%6d %3d %6d %6d %6d %6d %6d",tor[i].id+1,tor[i].ffid+1,tor[i].atm[0]+1,tor[i].atm[1]+1,tor[i].atm[2]+1,tor[i].atm[3]+1,tor[i].deg);
          *outf<<null<<endl;
      }
      *outf<<endl;
    }

    if(ninv) {
       sprintf(null,"Inversions");
       *outf<<null<<endl<<endl;
       for(i=0;i<ninv;i++) {
          sprintf(null,"%6d %3d %6d %6d %6d %6d",inv[i].id+1,inv[i].ffid+1,inv[i].atm[0]+1,inv[i].atm[1]+1,inv[i].atm[2]+1,inv[i].atm[3]+1);
          *outf<<null<<endl;
      }
      *outf<<endl;
    }

    if(nimp) {
       sprintf(null,"Impropers");
       *outf<<null<<endl<<endl;
       for(i=0;i<nimp;i++) {
          sprintf(null,"%6d %3d %6d %6d %6d %6d",imp[i].id+1,imp[i].ffid+1,imp[i].atm[0]+1,imp[i].atm[1]+1,imp[i].atm[2]+1,imp[i].atm[3]+1);
          *outf<<null<<endl;
      }
      *outf<<endl;
    }

    return 0;
}

int MODEL::rd_bgf(char *indata)
{
    cout<<"Reading BGF file "<<indata<<endl;
    ifstream inf(indata,fstream::in);
    if(!inf.is_open()) {
      cout<<" Error: BGF file "<<indata<<" cannot be opened."<<endl;
      return filenotopen;
    }
    
    char null[1024];
    int i,j;  
    char delims[] = {" ;,\t"};
    char *charptr;

    while(!inf.eof()) {
       inf.getline(null,1024);
       if (strncmp(null,"CRYSTX",6) == 0) break;
       else if ((strncmp(null,"HETATM",6) == 0) || (strncmp(null,"ATOM",4) == 0) ) break;
    }

    if(strncmp(null,"CRYSTX",6) == 0) {
       periodic = 3;
       charptr = strtok(null,delims);
       charptr = strtok(NULL, delims);
       cell.la=atof(charptr);  
       charptr = strtok(NULL, delims);
       cell.lb=atof(charptr);  
       charptr = strtok(NULL, delims);
       cell.lc=atof(charptr);  
       charptr = strtok(NULL, delims);
       cell.alpha=atof(charptr);  
       charptr = strtok(NULL, delims);
       cell.beta=atof(charptr);  
       charptr = strtok(NULL, delims);
       cell.gamma=atof(charptr);  
       while(!inf.eof()) {
          inf.getline(null,1024);
          if ((strncmp(null,"HETATM",6) == 0) || (strncmp(null,"ATOM",4) == 0) ) break;
       }
    } else {
       periodic = 0;
       cell.la=cell.lb=cell.lc=cell.alpha=cell.beta=cell.gamma=0;
    }
    //printf(" periodic %d a %f b %f c %f alpha %f beta %f gamma %f\n",periodic,cell.la,cell.lb,cell.lc,cell.alpha,cell.beta,cell.gamma);
    cell.labc2H();
    cell.H2others();

    if((strncmp(null,"HETATM",6) == 0) || (strncmp(null,"ATOM",4) == 0)) {
      natom=1;
      while(!inf.eof()) {
         inf.getline(null,1024);
         if ((strncmp(null,"HETATM",6) == 0) || (strncmp(null,"ATOM",4) == 0) ) natom++;
         else if(strncmp(null,"CONECT",6) == 0) break;
      }
    } else {
      natom=0;
      return 1;
    }
    init_atom();
    inf.clear();
    inf.seekg(0,ios::beg);
    inf>>null;
    while((strncmp(null,"HETATM",6) != 0) && (strncmp(null,"ATOM",4) != 0)) inf>>null;
    if(strncmp(null,"ATOM",4)==0) { inf.getline(null,1024); inf>>null;}
    int Symbsize=2;
    //get atomic info
    for(i=0;i<natom;i++) {
       atom[i].id=i;
       inf>>null>>atom[i].name>>null>>null>>null
          >>atom[i].pv[0]>>atom[i].pv[1]>>atom[i].pv[2]
          >>atom[i].fftype>>null>>null>>atom[i].chg>>null;
       //remove number in the atom name
        for(j=0;j<Symbsize;j++)  
          if(atom[i].name[j]>=48 && atom[i].name[j]<=57)  {
             for(;j<Symbsize;j++) atom[i].name[j]='\0';
          }
    }
    //get connectivity
    while(strncmp(null,"CONECT",6) != 0 && !inf.eof()) inf.getline(null,1024);
    i=0;
    while(!inf.eof()&&strncmp(null,"END",3)!=0) {
       if(strncmp(null,"CONECT",6)==0) {
          charptr = strtok(null,delims);
          charptr = strtok(NULL, delims);
          i=atoi(charptr)-1;
          atom[i].ncnt=0;
          charptr = strtok(NULL, delims);
          while(charptr!=NULL) {
              atom[i].bod[ atom[i].ncnt ] =1; //default bond order is 1
              atom[i].cnt[ atom[i].ncnt++ ]=atoi(charptr)-1;
              charptr = strtok(NULL, delims);
          }
       } else if(strncmp(null,"ORDER",5)==0) {
          charptr = strtok(null,delims);
          charptr = strtok(NULL, delims);
          if(i!=atoi(charptr)-1) {
             cout<<"Error: mismatch in rd_bgf cnect"<<endl;
             return 1;
          }
          for(j=0;j<atom[i].ncnt;j++) {
             charptr = strtok(NULL, delims);
             atom[i].bod[j]=atoi(charptr);
          }
       }
       inf.getline(null,1024);
    }
    if(0) {
      for(i=0;i<natom;i++) {
        printf("%6s %5d %5s RES A  444 %10.5f%10.5f%10.5f %-5s  1 0  %7.5f\n",
        "HETATM",i+1,atom[i].name,atom[i].pv[0],atom[i].pv[1],atom[i].pv[2],atom[i].fftype,atom[i].chg);
      }
      for(i=0;i<natom;i++) {
        printf("%6s%6d","CONECT",i+1);
        for(j=0;j<atom[i].ncnt;j++) printf("%6d",atom[i].cnt[j]+1);
        printf("\n");
        printf("%6s%6d","ORDER",i+1);
        for(j=0;j<atom[i].ncnt;j++) printf("%6d",atom[i].bod[j]);
        printf("\n");
      }
        
    }
    return 0;
}

int MODEL::out_bgf(ostream *outf)
{
    
    char null[1024],tmp[10];
    int i,j;  

    if(periodic) {
      sprintf(null,"XTLGRF 200");
    } else {
      sprintf(null,"BIOGRF 200");
    }
    *outf<<null<<endl;
    sprintf(null,"DESCRP %s","NULL");
    *outf<<null<<endl;
    sprintf(null,"REMARK %s","BGF file created by MDANALYSIS");
    *outf<<null<<endl;
    sprintf(null,"FORCEFIELD %s",name);
    *outf<<null<<endl;
    if(periodic) {
      sprintf(null,"PERIOD 111");
      *outf<<null<<endl;
      sprintf(null,"AXES   ZYX");
      *outf<<null<<endl;
      sprintf(null,"SGNAME P 1");
      *outf<<null<<endl;
      cell.H2others();
      sprintf(null,"CRYSTX %11.5f%11.5f%11.5f%11.5f%11.5f%11.5f",cell.la,cell.lb,cell.lc,cell.alpha,cell.beta,cell.gamma);
      *outf<<null<<endl;
      sprintf(null,"CELLS %5d%5d%5d%5d%5d%5d",-1,1,-1,1,-1,1);
      *outf<<null<<endl;
    } 
   
    sprintf(null,"FORMAT ATOM   (%s)","a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5,1x,a5, i3,i2,1x,f8.5");
    *outf<<null<<endl;
    for(i=0;i<natom;i++) {
      sprintf(tmp,"%s%d",atom[i].name,i+1); tmp[5]='\0';
      sprintf(null,"HETATM %5d %-5s RES A  444 %10.5f%10.5f%10.5f %-5s  1 0 %8.5f",i+1,tmp,atom[i].pv[0],atom[i].pv[1],atom[i].pv[2],atom[i].fftype,atom[i].chg);
      *outf<<null<<endl;
    }
   
    sprintf(null,"FORMAT CONECT (%s)","a6,12i6");
    *outf<<null<<endl;
    for(i=0;i<natom;i++) {
      sprintf(null,"CONECT%6d",i+1);
      for(j=0;j<atom[i].ncnt;j++) sprintf(null,"%s%6d",null,atom[i].cnt[j]+1);
      *outf<<null<<endl;
      sprintf(null,"ORDER %6d",i+1);
      for(j=0;j<atom[i].ncnt;j++) sprintf(null,"%s%6d",null,atom[i].bod[j]);
      *outf<<null<<endl;

    }
    sprintf(null,"END");
    *outf<<null<<endl;

    return 0;
}

int MODEL::rd_gro(char *indata)
{
    cout<<"Reading GROMACS gro file "<<indata<<endl;
    ifstream inf(indata,fstream::in);
    if(!inf.is_open()) {
      cout<<" Error: GROMACS gro file "<<indata<<" cannot be opened."<<endl;
      return filenotopen;
    }
    
    char null[1024],nxll[1024];
    int i,j;  
    char delims[] = {" ;,\t"};
    char *charptr;

    inf.getline(null,1024);
    inf>>natom;
    init_atom();
    inf.getline(null,1024);
    //get atomic info
    for(i=0;i<natom;i++) {
       atom[i].id=i;
       inf.getline(null,1024);
       strncpy(nxll,&null[10],5); nxll[5]='\0'; 
       for(j=0;j<5;j++) if(nxll[j]!=' ') break;
       strncpy(atom[i].name,&nxll[j],1);
       strncpy(atom[i].fftype,&nxll[j],5-j);
       strncpy(nxll,&null[20],8); nxll[8]='\0'; atom[i].pv[0]=atof(nxll); atom[i].pv[0]*=10;
       strncpy(nxll,&null[28],8); nxll[8]='\0'; atom[i].pv[1]=atof(nxll); atom[i].pv[1]*=10;
       strncpy(nxll,&null[36],8); nxll[8]='\0'; atom[i].pv[2]=atof(nxll); atom[i].pv[2]*=10;
       strncpy(nxll,&null[44],8); nxll[8]='\0'; atom[i].vel[0]=atof(nxll); atom[i].vel[0]*=10;
       strncpy(nxll,&null[52],8); nxll[8]='\0'; atom[i].vel[1]=atof(nxll); atom[i].vel[1]*=10;
       strncpy(nxll,&null[60],8); nxll[8]='\0'; atom[i].vel[2]=atof(nxll); atom[i].vel[2]*=10;
       atom[i].chg=0; //not available
    }
    inf>>cell.la>>cell.lb>>cell.lc;
    periodic = 3;
    cell.la*=10;
    cell.lb*=10;
    cell.lc*=10;
    cell.alpha=cell.beta=cell.gamma=90;
    
    //printf(" periodic %d a %f b %f c %f alpha %f beta %f gamma %f\n",periodic,cell.la,cell.lb,cell.lc,cell.alpha,cell.beta,cell.gamma);
    cell.labc2H();
    cell.H2others();

    return 0;
}

int MODEL::rd_g96(char *indata)
{
    cout<<"Reading GROMACS g96 file "<<indata<<endl;
    ifstream inf(indata,fstream::in);
    if(!inf.is_open()) {
      cout<<" Error: GROMACS g96 file "<<indata<<" cannot be opened."<<endl;
      return filenotopen;
    }

    char null[1024],nxll[1024];
    int i,j,k,l;
    char delims[] = {" ;,\t"};
    char *charptr;
    //get natom    
    natom=-1;
    do{
    inf.getline(null,1024);
    } while(strncmp(null,"POSITION",8));
    do{
    inf.getline(null,1024);
    natom++;
    } while(strncmp(null,"END",3));
    cout<<"natom "<<natom<<endl;
    init_atom();
    inf.seekg(0);

    //get atomic info
    do{
    inf.getline(null,1024);
    } while(strncmp(null,"POSITION",8));
    nmol=0;
    for(i=0;i<natom;i++) {
       atom[i].id=i;
       inf.getline(null,1024);
       strncpy(nxll,&null[0],5); nxll[5]='\0'; //molecule name in nxll
       atom[i].mol=atoi(nxll)-1;
       if(atom[i].mol>nmol) nmol=atom[i].mol+1;
       strncpy(nxll,&null[6],5); nxll[5]='\0'; //molecule name in nxll
       bool chk=0;
       for(j=0;j<5;j++)  if(nxll[j]!=' ') break;
       for(k=j;k<5;k++)  if(nxll[k]==' ') break;
       strncpy(atom[i].molname,&nxll[j],k-j);
       strncpy(nxll,&null[12],5); nxll[5]='\0'; //ff type
       for(j=0;j<5;j++) if(nxll[j]!=' ') break;
       for(k=0;k<5;k++)  {
           if(int(nxll[k])<=57 and int(nxll[k])>=48 and chk==0 ) {l=k;chk=1;}
           else if(int(nxll[k])<=122 and int(nxll[k])>=65 ) {chk=0;}
           else if (nxll[k]==' ') break;
       }
       if(nxll[j+1]>96&&nxll[j+1]<123) strncpy(atom[i].name,&nxll[j],2);
       else strncpy(atom[i].name,&nxll[j],1);
       strncpy(atom[i].fftype,&nxll[j],l-j);
       atom[i].fftype[l-j]='\0';
       /*if (ctl.addffid_flag==1) { if (strncmp(atom[i].fftype,ctl.ana_water_refatmtpname[0],2)==0) atom[i].ffid=ctl.ana_water_refatmtp[0]-1;
                                  else if (strncmp(atom[i].fftype,ctl.ana_water_refatmtpname[1],2)==0) atom[i].ffid=ctl.ana_water_refatmtp[1]-1;
                                  else atom[i].ffid=2;}*/
       strncpy(nxll,&null[24],15); nxll[15]='\0'; atom[i].pv[0]=atof(nxll); atom[i].pv[0]*=10;
       strncpy(nxll,&null[39],15); nxll[15]='\0'; atom[i].pv[1]=atof(nxll); atom[i].pv[1]*=10;
       strncpy(nxll,&null[54],15); nxll[15]='\0'; atom[i].pv[2]=atof(nxll); atom[i].pv[2]*=10;
       atom[i].chg=0; //not available
    }
    inf.getline(null,1024);
    inf.getline(null,1024);
    if(strcmp(null,"VELOCITY")==0){
       for(i=0;i<natom;i++) {
           inf.getline(null,1024);
           strncpy(nxll,&null[10],5); nxll[5]='\0';
           for(j=0;j<5;j++) if(nxll[j]!=' ') break;
           strncpy(nxll,&null[24],15); nxll[15]='\0'; atom[i].vel[0]=atof(nxll); atom[i].vel[0]*=10;
           strncpy(nxll,&null[39],15); nxll[15]='\0'; atom[i].vel[1]=atof(nxll); atom[i].vel[1]*=10;
           strncpy(nxll,&null[54],15); nxll[15]='\0'; atom[i].vel[2]=atof(nxll); atom[i].vel[2]*=10;
       }
       inf.getline(null,1024);
       inf.getline(null,1024);
    }
    inf>>cell.la>>cell.lb>>cell.lc;
    periodic = 3;
    cell.la*=10;
    cell.lb*=10;
    cell.lc*=10;
    cell.alpha=cell.beta=cell.gamma=90;
    printf(" periodic %d a %f b %f c %f alpha %f beta %f gamma %f\n",periodic,cell.la,cell.lb,cell.lc,cell.alpha,cell.beta,cell.gamma);
    cell.labc2H();
    cell.H2others();

    return 0;
}

int MODEL::out_g96(ostream *outf)
{
    int i;
    char null[1024];

    *outf<<"TITLE"<<endl;
    *outf<<"G96_structure"<<endl;
    *outf<<"END"<<endl;

    *outf<<"POSITION"<<endl;
    for(i=0;i<natom;i++) {
        sprintf(null,"%5d %-5s %-5s %6d%15.9f%15.9f%15.9f",atom[i].mol+1,atom[i].molname,atom[i].fftype,i+1,atom[i].pv[0]*0.1,atom[i].pv[1]*0.1,atom[i].pv[2]*0.1);
        *outf<<null<<endl;
    }
    *outf<<"END"<<endl;
    //cout<<"stlin c="<<c<<endl;

    *outf<<"VELOCITY"<<endl;
    for(i=0;i<natom;i++) {
        sprintf(null,"%5d %-5s %-5s %6d%15.9f%15.9f%15.9f",atom[i].mol+1,atom[i].molname,atom[i].fftype,i+1,atom[i].vel[0]*0.1,atom[i].vel[1]*0.1,atom[i].vel[2]*0.1);
        *outf<<null<<endl;
    }
    *outf<<"END"<<endl;

    *outf<<"BOX"<<endl;
    sprintf(null,"%15.9f%15.9f%15.9f",cell.la*0.1,cell.lb*0.1,cell.lc*0.1);
    *outf<<null<<endl;
    *outf<<"END"<<endl;

    return 0;
}

int MODEL::rd_cpmdin(char *indata)
{
    cout<<"Reading CPMD in file "<<indata<<endl;
    ifstream inf(indata,fstream::in);
    if(!inf.is_open()) {
      cout<<" Error: CPMD in file "<<indata<<" cannot be opened."<<endl;
      return filenotopen;
    }

    char null[1024],tmp[3];
    int i,j,k;
    char delims[] = {" ;,\t"};
    char *charptr;

    //find atoms
    natom=0;
    do{ inf>>null; } while (strcmp(null,"&ATOMS")!=0); 
    while(!inf.eof()) {
      if (strncmp(null,"*",1) == 0) {
         inf>>null>>null;
         natom+=atoi(null);
      } 
      inf.getline(null,1024);
    }
    //cout<<natom<<endl;
    init_atom();
    inf.clear();
    inf.seekg(0,ios::beg);
    do{ inf>>null; } while (strcmp(null,"&ATOMS")!=0);
    k=0;
    while(!inf.eof()) {
      if (strncmp(null,"*",1) == 0) {
         strncpy(tmp,&null[1],2); if(tmp[1]='_') tmp[1]='\0';
         inf>>null>>j;
         for(i=0;i<j;i++) {
            inf>>atom[k].pv[0]>>atom[k].pv[1]>>atom[k].pv[2];
            //atom[k].pv[0]*=Bohr2A;
            //atom[k].pv[1]*=Bohr2A;
            //atom[k].pv[2]*=Bohr2A;
            atom[k].id=k;
            strcpy(atom[k].name,tmp);
            k++;
         }
      }
      inf.getline(null,1024);
    }
    if(k!=natom) cout<<"Error in reading atom"<<endl;
    //for(i=0;i<natom;i++) {
    //   printf("%d %s %f %f %f\n",i+1,atom[i].name,atom[i].pv[0],atom[i].pv[1],atom[i].pv[2]);
    //}

    //read cell info
    inf.clear();
    inf.seekg(0,ios::beg);
    do{ inf>>null; } while (strcmp(null,"&SYSTEM")!=0);
    do{ inf>>null; } while (strcmp(null,"CELL")!=0);
    periodic=3; cell.o[0]=cell.o[1]=cell.o[2]=0;
    inf>>cell.la>>cell.lb>>cell.lc>>cell.alpha>>cell.beta>>cell.gamma;
    //cell.la*=Bohr2A; 
    cell.lb*=cell.la; cell.lc*=cell.la;
    double fac=180.0/acos(-1.0);
    cell.alpha=fac*acos(cell.alpha); cell.beta=fac*acos(cell.beta); cell.gamma=fac*acos(cell.gamma);
    cell.labc2H();
    cell.H2others();

    //read trj info
    inf.clear();
    inf.seekg(0,ios::beg);
    do{ inf>>null; } while (strcmp(null,"&CPMD")!=0);
    do{ inf>>null; } while (strcmp(null,"XYZ")!=0);
    inf>>xyzfreq; velfreq=xyzfreq;
    do{ inf>>null; } while (strcmp(null,"TIMESTEP")!=0);
    inf>>timestep; timestep*=0.0000241888428; //in ps
    do{ inf>>null; } while (strcmp(null,"STRESS")!=0);
    inf>>null>>strfreq; 
    engfreq=1;

    inf.close();
    return 0;
}

int MODEL::out_cpmdin(ostream *outf)
{
    return 0;
}

int MODEL::bond2mol()
{
    int i,j,k,ck;
    int tota,fstatm;
    int *tmatom=new int [natom];

    nmol=0;
    fstatm=0;
    for(i=0;i<natom;i++) atom[i].mol=-1;
    atom[fstatm].mol=nmol;
    tmatom[nmol]=1;
    tota=0;
    while(tmatom[nmol]) {
       MOLECULE tmol; tmol.natom=natom; tmol.init_atom();
       tmol.atm[0]=fstatm;
       for(i=0;i<tmatom[nmol];i++) {
          for(j=0;j<nbond;j++) {
             ck=-1;
             if(bond[j].atm[0]==tmol.atm[i]) ck=bond[j].atm[1];
             else if(bond[j].atm[1]==tmol.atm[i]) ck=bond[j].atm[0];
             if(ck>0) {
               //check if bonded atom in the list
               for(k=0;k<tmatom[nmol];k++) if(tmol.atm[k]==ck) break;
               if(k==tmatom[nmol]) { //add new atom to the list
                  atom[ck].mol=nmol;
                  tmol.atm[tmatom[nmol]++]=ck;
               }
             } 
          }
       }
       //printf("stlin mol %d natom %d first atom %d\n",nmol+1,tmatom[nmol],fstatm+1);
       tota+=tmatom[nmol];
       nmol++;
       for(i=0;i<natom;i++) if(atom[i].mol<0) break;
       if(i<natom) {
         tmatom[nmol]=1;
         fstatm=i;
         atom[fstatm].mol=nmol;
       } else tmatom[nmol]=0;
    }
    if(tota!=natom) {
      cout<<" Error: bond2mol failed natom="<<natom<<" tota="<<tota<<endl;
      delete [] tmatom;
      return 1;
    }
    init_mol();
    for(i=0;i<nmol;i++) {
       mol[i].natom=tmatom[i];
       mol[i].init_atom();
       mol[i].natom=0;
       mol[i].mass=0;
    }
    for(i=0;i<natom;i++) {
       j=atom[i].mol;
       mol[j].atm[ mol[j].natom ]=i;
       mol[j].atom[ mol[j].natom++ ]= & atom[i];
       mol[j].mass+=atom[i].mass;
    }  

    if(0) {
      for(i=0;i<nmol;i++) {
        for(j=0;j<mol[i].natom;j++) cout<<" "<<mol[i].atm[j];
        cout<<endl;
      }
    }
    delete [] tmatom;
    return 0;
}

int MODEL::bond2cnt()
{
    //clear connectivity
    int i,j,a1,a2;
    for(i=0;i<natom;i++) {
       atom[i].ncnt=0;
       for(j=0;j<maxcnt;j++) { atom[i].cnt[j]=-1; atom[i].bod[j]=1; }
    }
    for(i=0;i<nbond;i++) {
       a1=bond[i].atm[0];
       a2=bond[i].atm[1];
       atom[a1].cnt[ atom[a1].ncnt++ ]=a2;
       atom[a2].cnt[ atom[a2].ncnt++ ]=a1;
    }
    if(0) {
      for(i=0;i<natom;i++) {
        printf("%6s %5d %5s RES A  444 %10.5f%10.5f%10.5f %-5s  1 0  %7.5f\n",
        "HETATM",i+1,atom[i].name,atom[i].pv[0],atom[i].pv[1],atom[i].pv[2],atom[i].fftype,atom[i].chg);
      }
      for(i=0;i<natom;i++) {
        printf("%6s%6d","CONECT",i+1);
        for(j=0;j<atom[i].ncnt;j++) printf("%6d",atom[i].cnt[j]+1);
        printf("\n");
        printf("%6s%6d","ORDER",i+1);
        for(j=0;j<atom[i].ncnt;j++) printf("%6d",atom[i].bod[j]);
        printf("\n");
      }
    }
    return 0;
}

int MODEL::cnt2valence()
{   //bgf2lmp
    //MakeLists()  Building lists

    int i,j,p,q,k,id1,id2,id3;
    //find bond
    nbond=0;
    for(i=0;i<natom;i++) {
       for(j=0;j<atom[i].ncnt;j++) {
           if(atom[i].id<atom[atom[i].cnt[j]].id) nbond++; 
           //this avoids counting the same bond twice
       }
    }
    init_bond();
    k=0;
    for(i=0;i<natom;i++) {
       for(j=0;j<atom[i].ncnt;j++) {
           id1=atom[i].cnt[j];
           if(atom[i].id<atom[id1].id) {
              bond[k].id=k;
              bond[k].atm[0]=atom[i].id;
              bond[k].atm[1]=atom[id1].id;
              bond[k].atom[0]=&atom[i];
              bond[k].atom[1]=&atom[id1];
              k++;
           }
       }
    }

    //find angle
    nangle=0;
    for(i=0;i<natom;i++) {
       if(atom[i].ncnt>1) {
          nangle += (atom[i].ncnt*(atom[i].ncnt-1))/2;
       }
    }
    init_angle();
    k=0;
    for(i=0;i<natom;i++) { // id1 i id2
       if(atom[i].ncnt>1) {  //central atom
          for(p=0;p<atom[i].ncnt;p++) {
             id1=atom[i].cnt[p];
             for(q=p+1;q<atom[i].ncnt;q++) {
                id2=atom[i].cnt[q];
                if(atom[id1].id>atom[id2].id) { id2=id1; id1=atom[i].cnt[q]; }
                angle[k].id=k;
                angle[k].atm[0]=atom[id1].id;
                angle[k].atm[1]=atom[i].id;  
                angle[k].atm[2]=atom[id2].id;
                angle[k].atom[0]=&atom[id1];
                angle[k].atom[1]=&atom[i];
                angle[k].atom[2]=&atom[id2];
                k++;
             }
          }
       }
    }
   
    //find torsion
    ntor=0;
    for(i=0;i<natom;i++) { //p(id2) i j(id1) q(id3)
       if(atom[i].ncnt>1) {
          for(j=0;j<atom[i].ncnt;j++) {
             id1=atom[i].cnt[j];
             if((atom[id1].id>atom[i].id) && (atom[id1].ncnt>1)) {
                for(p=0;p<atom[i].ncnt;p++) {
                   id2=atom[i].cnt[p];
                   if( id2!=id1 ) {
                      for(q=0;q<atom[id1].ncnt;q++) {
                         id3=atom[id1].cnt[q];
                         if(id3 != i && id3 != id2) ntor++;
                      }
                   }
                }
             }
          }
       }
    }
    init_tor();
    k=0;
    for(i=0;i<natom;i++) { //p(id2) i j(id1) q(id3)
       if(atom[i].ncnt>1) {
          for(j=0;j<atom[i].ncnt;j++) {
             id1=atom[i].cnt[j];
             if(atom[id1].id>atom[i].id && atom[id1].ncnt>1) {
                for(p=0;p<atom[i].ncnt;p++) {
                   id2=atom[i].cnt[p];
                   if( id2!=id1 ) {
                      for(q=0;q<atom[id1].ncnt;q++) {
                         id3=atom[id1].cnt[q];
                         if(id3 != i && id3 != id2) {
                            tor[k].id=k;
                            tor[k].atm[0]=atom[id2].id;
                            tor[k].atm[1]=atom[i].id;
                            tor[k].atm[2]=atom[id1].id;
                            tor[k].atm[3]=atom[id3].id;
                            tor[k].atom[0]=&atom[id2];
                            tor[k].atom[1]=&atom[i];
                            tor[k].atom[2]=&atom[id1];
                            tor[k].atom[3]=&atom[id3];
                            tor[k].deg=(atom[i].ncnt-1)*(atom[id1].ncnt-1);
                            // The degeneracy of the dihedral is calculated as 
                            // the number of connections of atom i minus one, times
                            // the number of connections of atom id1 minus one.  
                            k++;
                         }
                      }
                   }
                }
             }
          }
       }
    }

    //find inversion
    ninv=0;
    for(i=0;i<natom;i++) {  // central atom
       if(atom[i].ncnt==3) ninv+=3;
    }
    init_inv();
    k=0;
    for(i=0;i<natom;i++) {  // central atom   i id1 id2 id3
       if(atom[i].ncnt==3) {
          id1=atom[i].cnt[0];
          id2=atom[i].cnt[1];
          id3=atom[i].cnt[2];
   
          inv[k].id=k;
          inv[k].atm[0]=atom[i].id;
          inv[k].atm[1]=atom[id1].id;
          inv[k].atm[2]=atom[id2].id;
          inv[k].atm[3]=atom[id3].id;
          inv[k].atom[0]=&atom[i];
          inv[k].atom[1]=&atom[id1];
          inv[k].atom[2]=&atom[id2];
          inv[k].atom[3]=&atom[id3];
          k++;

          inv[k].id=k;
          inv[k].atm[0]=atom[i].id;
          inv[k].atm[1]=atom[id2].id;
          inv[k].atm[2]=atom[id3].id;
          inv[k].atm[3]=atom[id1].id;
          inv[k].atom[0]=&atom[i];
          inv[k].atom[1]=&atom[id2];
          inv[k].atom[2]=&atom[id3];
          inv[k].atom[3]=&atom[id1];
          k++;

          inv[k].id=k;
          inv[k].atm[0]=atom[i].id;
          inv[k].atm[1]=atom[id3].id;
          inv[k].atm[2]=atom[id1].id;
          inv[k].atm[3]=atom[id2].id;
          inv[k].atom[0]=&atom[i];
          inv[k].atom[1]=&atom[id3];
          inv[k].atom[2]=&atom[id1];
          inv[k].atom[3]=&atom[id2];
          k++;
       }
    }
    //printf("nbond %d nangle %d ntor %d ninv %d\n",nbond,nangle,ntor,ninv);
    return 0;
}

int MODEL::cnt2bod()
{
    int i,j,k,id;
    int *ckf = new int [natom];
    
    //check BOD for H
    for(i=0;i<natom;i++) {
      if(strcmp(atom[i].name,"H")==0) {
        if(atom[i].ncnt>1) printf("Warning CNT (%d) of atom %d element H may be wrong\n",atom[i].ncnt,i+1);
        else {
          atom[i].bod[0]=1;
          id=atom[i].cnt[0];
          for(j=0;j<atom[id].ncnt;j++)
            if(atom[id].cnt[j]==i) {
               atom[id].bod[j]=1;
               break;
            }
        }
        ckf[i]=1;
      } else ckf[i]=0;
    }

    //check BOD for O
    for(i=0;i<natom;i++) {
      if(ckf[i]) continue;
      if(strcmp(atom[i].name,"O")==0) {
        if(atom[i].ncnt==1) {
          atom[i].bod[0]=2;
          id=atom[i].cnt[0];
          for(j=0;j<atom[id].ncnt;j++) 
            if(atom[id].cnt[j]==i) { 
               atom[id].bod[j]=2;
               break;
            }
        } else if(atom[i].ncnt==2) {
          for(k=0;k<atom[i].ncnt;k++) {
             atom[i].bod[k]=1;
             id=atom[i].cnt[k];
             for(j=0;j<atom[id].ncnt;j++)
               if(atom[id].cnt[j]==i) {
                  atom[id].bod[j]=1;
                  break;
               }
          }

        } else if(atom[i].ncnt==3) {
          for(k=0;k<atom[i].ncnt;k++) {
             atom[i].bod[k]=1;
             id=atom[i].cnt[k];
             for(j=0;j<atom[id].ncnt;j++)
               if(atom[id].cnt[j]==i) {
                  atom[id].bod[j]=1;
                  break;
               }
          }
        } else printf("Warning CNT (%d) of atom %d element O may be wrong\n",atom[i].ncnt,i+1);
        ckf[i]=1;
      } 
    }

    //check BOD for N
    for(i=0;i<natom;i++) {
      if(ckf[i]) continue;
      if(strcmp(atom[i].name,"N")==0) {
        if(atom[i].ncnt==1) {
          atom[i].bod[0]=3;
          id=atom[i].cnt[0];
          for(j=0;j<atom[id].ncnt;j++)
            if(atom[id].cnt[j]==i) {
               atom[id].bod[j]=3;
               break;
            }
        } else if(atom[i].ncnt==2) {
          if( (atom[i].bod[0]+atom[i].bod[1]) !=3) {
              if(ckf[atom[i].cnt[0]]&&ckf[atom[i].cnt[1]]) {
                printf("Warning failed to determine BOD of atom %d element N\n",i+1);
              } else if (ckf[atom[i].cnt[0]]) {
                int bo=3-atom[i].bod[0];
                atom[i].bod[1]=bo;
                id=atom[i].cnt[1];
                for(j=0;j<atom[id].ncnt;j++)
                  if(atom[id].cnt[j]==i) {
                     atom[id].bod[j]=bo;
                     break;
                  }
              } else if (ckf[atom[i].cnt[1]]) {
                int bo=3-atom[i].bod[1];
                atom[i].bod[0]=bo;
                id=atom[i].cnt[0];
                for(j=0;j<atom[id].ncnt;j++)
                  if(atom[id].cnt[j]==i) {
                     atom[id].bod[j]=bo;
                     break;
                  }
              } else {
                int a1,a2;
                if(strcmp(atom[ atom[i].cnt[0] ].name,"N")==0) { a1=0;a2=1; }
                else { a1=1;a2=0; } //pick any non N atom to be double bond
                atom[i].bod[a1]=1;
                atom[i].bod[a2]=2;
                id=atom[i].cnt[a1];
                for(j=0;j<atom[id].ncnt;j++)
                  if(atom[id].cnt[j]==i) {
                     atom[id].bod[j]=1;
                     break;
                  }
                id=atom[i].cnt[a2];
                for(j=0;j<atom[id].ncnt;j++)
                  if(atom[id].cnt[j]==i) {
                     atom[id].bod[j]=2;
                     break;
                  }
              }
          }
        } else if(atom[i].ncnt==3) {
          for(k=0;k<atom[i].ncnt;k++) {
             atom[i].bod[k]=1;
             id=atom[i].cnt[k];
             for(j=0;j<atom[id].ncnt;j++)
               if(atom[id].cnt[j]==i) {
                  atom[id].bod[j]=1;
                  break;
               }
          }
        } else printf("Warning CNT (%d) of atom %d element N may be wrong\n",atom[i].ncnt,i+1);
        ckf[i]=1;
      } 
    }

    //check BOD for C
    for(i=0;i<natom;i++) {
      if(ckf[i]) continue;
      if(strcmp(atom[i].name,"C")==0) {
        if(atom[i].ncnt==2) {
          //22,31
          if( (atom[i].bod[0]+atom[i].bod[1]) !=4) {
              if(ckf[atom[i].cnt[0]]&&ckf[atom[i].cnt[1]]) {
                printf("Warning failed to determine BOD of atom %d element C 2\n",i+1);
              } else if (ckf[atom[i].cnt[0]]) {
                int bo=4-atom[i].bod[0];
                atom[i].bod[1]=bo;
                id=atom[i].cnt[1];
                for(j=0;j<atom[id].ncnt;j++)
                  if(atom[id].cnt[j]==i) {
                     atom[id].bod[j]=bo;
                     break;
                  }
              } else if (ckf[atom[i].cnt[1]]) {
                int bo=4-atom[i].bod[1];
                atom[i].bod[0]=bo;
                id=atom[i].cnt[0];
                for(j=0;j<atom[id].ncnt;j++)
                  if(atom[id].cnt[j]==i) {
                     atom[id].bod[j]=bo;
                     break;
                  }
              } else {
                int a1,a2;
                if(strcmp(atom[ atom[i].cnt[0] ].name,"C")==0) { a1=0;a2=1; }
                else if(strcmp(atom[ atom[i].cnt[1] ].name,"C")==0) { a1=1;a2=0; } 
                else { a1=0;a2=1; }
                atom[i].bod[a1]=5-atom[atom[i].cnt[a1]].ncnt;
                atom[i].bod[a2]=4-atom[i].bod[a1];
                id=atom[i].cnt[a1];
                for(j=0;j<atom[id].ncnt;j++)
                  if(atom[id].cnt[j]==i) {
                     atom[id].bod[j]=atom[i].bod[a1];
                     break;
                  }
                id=atom[i].cnt[a2];
                for(j=0;j<atom[id].ncnt;j++)
                  if(atom[id].cnt[j]==i) {
                     atom[id].bod[j]=atom[i].bod[a2];
                     break;
                  }
              }
          }
        } else if(atom[i].ncnt==3) {
          //112
          if( (atom[i].bod[0]+atom[i].bod[1]+atom[i].bod[2]) !=4) {
              if(ckf[atom[i].cnt[0]]&&ckf[atom[i].cnt[1]]&&ckf[atom[i].cnt[2]]) {
                printf("Warning failed to determine BOD of atom %d element C 3\n",i+1);
              } else if (ckf[atom[i].cnt[0]]&&ckf[atom[i].cnt[1]]) {
                int bo=4-atom[i].bod[0]-atom[i].bod[1];
                atom[i].bod[2]=bo;
                id=atom[i].cnt[2];
                for(j=0;j<atom[id].ncnt;j++)
                  if(atom[id].cnt[j]==i) {
                     atom[id].bod[j]=bo;
                     break;
                  }
              } else if (ckf[atom[i].cnt[0]]&&ckf[atom[i].cnt[2]]) {
                int bo=4-atom[i].bod[0]-atom[i].bod[2];
                atom[i].bod[1]=bo;
                id=atom[i].cnt[1];
                for(j=0;j<atom[id].ncnt;j++)
                  if(atom[id].cnt[j]==i) {
                     atom[id].bod[j]=bo;
                     break;
                  }
              } else if (ckf[atom[i].cnt[1]]&&ckf[atom[i].cnt[2]]) {
                int bo=4-atom[i].bod[1]-atom[i].bod[2];
                atom[i].bod[0]=bo;
                id=atom[i].cnt[0];
                for(j=0;j<atom[id].ncnt;j++)
                  if(atom[id].cnt[j]==i) {
                     atom[id].bod[j]=bo;
                     break;
                  }
              } else if (ckf[atom[i].cnt[0]]) {
                int a1=0,a2=1,a3=2;
                int b1,b2;
                if(atom[i].bod[a1]==2) { atom[i].bod[a2]=atom[i].bod[a3]=1; 
                   id=atom[i].cnt[a2];
                   for(j=0;j<atom[id].ncnt;j++)
                     if(atom[id].cnt[j]==i) {
                        atom[id].bod[j]=1;
                        break;
                     }
                   id=atom[i].cnt[a3];
                   for(j=0;j<atom[id].ncnt;j++)
                     if(atom[id].cnt[j]==i) {
                        atom[id].bod[j]=1;
                        break;
                     }
                } else {
                   if(strcmp(atom[ atom[i].cnt[a2] ].name,"C")==0) { b1=a2;b2=a3; }
                   else if(strcmp(atom[ atom[i].cnt[a3] ].name,"C")==0) { b1=a3;b2=a2; }
                   else { b1=a2;b2=a3; }
                   atom[i].bod[b1]=5-atom[atom[i].cnt[b1]].ncnt; if(atom[i].bod[b1]>2) atom[i].bod[b1]=2;
                   atom[i].bod[b2]=3-atom[i].bod[b1];
                   id=atom[i].cnt[b1];
                   for(j=0;j<atom[id].ncnt;j++)
                     if(atom[id].cnt[j]==i) {
                        atom[id].bod[j]=atom[i].bod[b1];
                        break;
                     }
                   id=atom[i].cnt[b2];
                   for(j=0;j<atom[id].ncnt;j++)
                     if(atom[id].cnt[j]==i) {
                        atom[id].bod[j]=atom[i].bod[b2];
                        break;
                     }
                }
              } else if (ckf[atom[i].cnt[1]]) {
                int a1=1,a2=0,a3=2;
                int b1,b2;
                if(atom[i].bod[a1]==2) { atom[i].bod[a2]=atom[i].bod[a3]=1;
                   id=atom[i].cnt[a2];
                   for(j=0;j<atom[id].ncnt;j++)
                     if(atom[id].cnt[j]==i) {
                        atom[id].bod[j]=1;
                        break;
                     }
                   id=atom[i].cnt[a3];
                   for(j=0;j<atom[id].ncnt;j++)
                     if(atom[id].cnt[j]==i) {
                        atom[id].bod[j]=1;
                        break;
                     }
                } else {
                   if(strcmp(atom[ atom[i].cnt[a2] ].name,"C")==0) { b1=a2;b2=a3; }
                   else if(strcmp(atom[ atom[i].cnt[a3] ].name,"C")==0) { b1=a3;b2=a2; }
                   else { b1=a2;b2=a3; }
                   atom[i].bod[b1]=5-atom[atom[i].cnt[b1]].ncnt; if(atom[i].bod[b1]>2) atom[i].bod[b1]=2;
                   atom[i].bod[b2]=3-atom[i].bod[b1];
                   id=atom[i].cnt[b1];
                   for(j=0;j<atom[id].ncnt;j++)
                     if(atom[id].cnt[j]==i) {
                        atom[id].bod[j]=atom[i].bod[b1];
                        break;
                     }
                   id=atom[i].cnt[b2];
                   for(j=0;j<atom[id].ncnt;j++)
                     if(atom[id].cnt[j]==i) {
                        atom[id].bod[j]=atom[i].bod[b2];
                        break;
                     }
                }
              } else if (ckf[atom[i].cnt[2]]) {
                int a1=2,a2=0,a3=1;
                int b1,b2;
                if(atom[i].bod[a1]==2) { atom[i].bod[a2]=atom[i].bod[a3]=1;
                   id=atom[i].cnt[a2];
                   for(j=0;j<atom[id].ncnt;j++)
                     if(atom[id].cnt[j]==i) {
                        atom[id].bod[j]=1;
                        break;
                     }
                   id=atom[i].cnt[a3];
                   for(j=0;j<atom[id].ncnt;j++)
                     if(atom[id].cnt[j]==i) {
                        atom[id].bod[j]=1;
                        break;
                     }
                } else {
                   if(strcmp(atom[ atom[i].cnt[a2] ].name,"C")==0) { b1=a2;b2=a3; }
                   else if(strcmp(atom[ atom[i].cnt[a3] ].name,"C")==0) { b1=a3;b2=a2; }
                   else { b1=a2;b2=a3; }
                   atom[i].bod[b1]=5-atom[atom[i].cnt[b1]].ncnt; if(atom[i].bod[b1]>2) atom[i].bod[b1]=2;
                   atom[i].bod[b2]=3-atom[i].bod[b1];
                   id=atom[i].cnt[b1];
                   for(j=0;j<atom[id].ncnt;j++)
                     if(atom[id].cnt[j]==i) {
                        atom[id].bod[j]=atom[i].bod[b1];
                        break;
                     }
                   id=atom[i].cnt[b2];
                   for(j=0;j<atom[id].ncnt;j++)
                     if(atom[id].cnt[j]==i) {
                        atom[id].bod[j]=atom[i].bod[b2];
                        break;
                     }
                }
              } else {
                int a1,a2,a3;
                if(strcmp(atom[ atom[i].cnt[0] ].name,"C")==0) { a1=0;a2=1;a3=2; }
                else if(strcmp(atom[ atom[i].cnt[1] ].name,"C")==0) { a1=1;a2=0;a3=2; }
                else if(strcmp(atom[ atom[i].cnt[2] ].name,"C")==0) { a1=2;a2=0;a3=1; }
                else { a1=0;a2=1;a3=2; } //pick any non C atom to be double bond
                atom[i].bod[a1]=2;
                atom[i].bod[a2]=1;
                atom[i].bod[a3]=1;
                id=atom[i].cnt[a1];
                for(j=0;j<atom[id].ncnt;j++)
                  if(atom[id].cnt[j]==i) {
                     atom[id].bod[j]=atom[i].bod[a1];
                     break;
                  }
                id=atom[i].cnt[a2];
                for(j=0;j<atom[id].ncnt;j++)
                  if(atom[id].cnt[j]==i) {
                     atom[id].bod[j]=atom[i].bod[a2];
                     break;
                  }
                id=atom[i].cnt[a3];
                for(j=0;j<atom[id].ncnt;j++)
                  if(atom[id].cnt[j]==i) {
                     atom[id].bod[j]=atom[i].bod[a3];
                     break;
                  }
              }
          }
        } else if(atom[i].ncnt==4) {
          for(k=0;k<atom[i].ncnt;k++) {
             atom[i].bod[k]=1;
             id=atom[i].cnt[k];
             for(j=0;j<atom[id].ncnt;j++)
               if(atom[id].cnt[j]==i) {
                  atom[id].bod[j]=1;
                  break;
               }
          }
        } else printf("Warning CNT (%d) of atom %d element C may be wrong\n",atom[i].ncnt,i+1);
        ckf[i]=1;
      } 
    }

    delete [] ckf;
    return 0;
}

int MODEL::dst2bond()
{
    double tol=1.15;
    double dist,cov;
    int i,j,k;

    vector<int> tbond[2];
    double sft[3],sft1[3];

    nbond=0; 
    for(i=0;i<natom;i++) {
       for(j=i+1;j<natom;j++) {
          for(k=0;k<3;k++) sft[k]=atom[j].pv[k]-atom[i].pv[k];  //place atom i to the center of the unit cell
          cell.mindst(sft,sft1);
          dist=0;
          for(k=0;k<3;k++) dist+=(sft1[k]*sft1[k]);
          dist=sqrt(dist);
/*
          printf("%10.4f %10.4f %10.4f\n",sft[0],sft[1],sft[2]);
          printf("%10.4f %10.4f %10.4f\n",sft1[0],sft1[1],sft1[2]);
*/
          cov=(atom[i].crad+atom[j].crad)*tol;
          if(dist<=cov) { //this is a covalent bond
             tbond[0].push_back(i); tbond[1].push_back(j); nbond++; 
          }
       }
    }
    cout<<"nbond="<<nbond<<endl;
    init_bond();
    for(i=0;i<nbond;i++) {
       bond[i].id=i;
       bond[i].atm[0]=tbond[0][i];
       bond[i].atm[1]=tbond[1][j]; 
       bond[i].lbd=-1;
    }

    return 0;
}

int MODEL::mol2bond()
{

    double tol=1.15;
    double dist,cov;
    int i,j,k,x,xi,xj;

    vector<int> tbond[2];
    double sft[3],sft1[3];

    MOLECULE tmol[nmol];
    for(i=0;i<natom;i++) {
        j=atom[i].mol;
        tmol[j].natom++; 
    }
    for(j=0;j<nmol;j++) { tmol[j].init_atom(); tmol[j].natom=0; }

    for(i=0;i<natom;i++) {
        j=atom[i].mol;
        tmol[j].atm[ tmol[j].natom ]=i;
        tmol[j].natom++;
    }

    nbond=0; 
    for(x=0;x<nmol;x++) {
       for(i=0;i<tmol[x].natom;i++) {
          xi=tmol[x].atm[i];
          for(j=i+1;j<tmol[x].natom;j++) {
             xj=tmol[x].atm[j];
             for(k=0;k<3;k++) sft[k]=atom[xj].pv[k]-atom[xi].pv[k];  //place atom i to the center of the unit cell
             cell.mindst(sft,sft1);
             dist=0;
             for(k=0;k<3;k++) dist+=(sft1[k]*sft1[k]);
             dist=sqrt(dist);
/*
             printf("%10.4f %10.4f %10.4f\n",sft[0],sft[1],sft[2]);
             printf("%10.4f %10.4f %10.4f\n",sft1[0],sft1[1],sft1[2]);
*/
             cov=(atom[xi].crad+atom[xj].crad)*tol;
             if(dist<=cov) { //this is a covalent bond
                tbond[0].push_back(xi); tbond[1].push_back(xj); nbond++; 
             }
          }
       }
    }
    cout<<"nbond="<<nbond<<endl;
    init_bond();
    for(i=0;i<nbond;i++) {
       bond[i].id=i;
       bond[i].atm[0]=tbond[0][i];
       bond[i].atm[1]=tbond[1][i]; 
       bond[i].lbd=-1;
    }
    return 0;
}

int MODEL::findlongbond()
{
    int i,k,a1,a2;
    double odist,cov;
    double tol=1.15;
    for(i=0;i<nbond;i++) {
       a1=bond[i].atm[0];
       a2=bond[i].atm[1];
       odist=0;
       for(k=0;k<3;k++) odist+=(atom[a2].pv[k]-atom[a1].pv[k])*(atom[a2].pv[k]-atom[a1].pv[k]);
       odist=sqrt(odist);
       cov=(atom[a1].crad+atom[a2].crad)*tol;
       if(odist>cov) { //bond crossing the cell boundary
          bond[i].lbd=1;
       }
    }

    return 0;
}

int MODEL::fixlongbond()
{
    int i,j,k; 
    double sft[3],sft1[3];
    bond2mol();
    int a1,a2,m1,w;

    for(i=0;i<nbond;i++) {
       if(atom[bond[i].atm[0]].mass>atom[bond[i].atm[1]].mass) a1=bond[i].atm[0];
       else a1=bond[i].atm[1]; //position of heavy atom will not be changed
       m1=atom[a1].mol;
       for(j=0;j<mol[m1].natom;j++) {
          if(mol[m1].atm[j]==a1) continue;
          a2=mol[m1].atm[j];
          for(k=0;k<3;k++) sft[k]=atom[a2].pv[k]-atom[a1].pv[k]+cell.cb[k]; //placing H atom at the center of unit cell
          cell.Map2UnitCell(sft,sft1);
          for(k=0;k<3;k++) { atom[a2].pv[k]+=sft1[k]; atom[a2].shift[k]=0; }
       }
    }
    return 0;
}

int MODEL::farneigh()
{
    int *visited=new int [natom+1]; //flag for visited atoms
    int *cntvstd=new int [natom+1];  //ith connected atom visited
    int *path=new int [natom+1];     //current search path
    int (*tcnt)[maxcnt]=new int [natom][maxcnt];
    int curatm,nxtatm,faratm,fardst,curdst;
    int i,j,k;
    //store atom connectivity in tcnt

    for(i=0;i<natom;i++) {
       for(j=0;j<natom;j++) { 
          for(k=0;k<atom[j].ncnt;k++) tcnt[j][k]=atom[j].cnt[k];
          visited[j]=cntvstd[j]=path[j]=0;
       }  visited[j]=cntvstd[j]=path[j]=0;
       faratm=curatm=i;
       fardst=curdst=0;
       while(cntvstd[curatm]<atom[curatm].ncnt) {
          path[curdst]=curatm;
          visited[curatm]=1;

          //next atom to visit
          if ( cntvstd[curatm]<atom[curatm].ncnt ) { 
             //move to next connected atom
             nxtatm=tcnt[curatm][ cntvstd[curatm]++ ];
             if(curdst>0 && nxtatm==path[curdst-1]) ;
             else curdst++;
             if(visited[nxtatm]==0) {
               //move current atom to its last connectivity
               for(j=0;j<atom[nxtatm].ncnt;j++) {
                  if(tcnt[nxtatm][j]==curatm) {
                     SWAP(&tcnt[nxtatm][j],&tcnt[nxtatm][atom[nxtatm].ncnt-1]);
                     break;
                  }
               }
             } else {
               //move back to avoid visiting the same atom twice
               nxtatm=path[--curdst];
             }
          } else {
             //move back
             nxtatm=path[--curdst];
          }
          curatm=nxtatm;
          if(curdst>fardst) { fardst=curdst; faratm=curatm; }
          //printf("stlin i %d curatm %d curdst %d faratm %d fardst %d cntvstd %d ncnt %d\n",i,curatm,curdst,faratm,fardst,cntvstd[curatm],atom[curatm].ncnt);
       }
       atom[i].faratm=faratm;
       atom[i].fardst=fardst;
       //printf("atom %d faratm %d fardst %d\n",i+1,atom[i].faratm+1,atom[i].fardst);
    }
    delete [] visited;
    delete [] cntvstd;
    delete [] path;
    delete [] tcnt;
    return 0;
}

int MODEL::headtail()
{
    int i,j;
    int fardst;
    for(i=0;i<nmol;i++) {
       fardst=0;
       for(j=0;j<mol[i].natom;j++) {
          if(atom[mol[i].atm[j]].fardst>fardst) {
             fardst=atom[mol[i].atm[j]].fardst;
             mol[i].head=mol[i].atm[j];
             mol[i].tail=atom[mol[i].atm[j]].faratm;
          }
       }
       if(mol[i].tail<mol[i].head) SWAP(&mol[i].tail,&mol[i].head);
       //printf("stlin mol %d head %d tail %d fardst %d\n",i+1,mol[i].head+1,mol[i].tail+1,fardst);
    }
    return 0;
}

int MODEL::rd_grp(char *ingrpf)
{
    int i,j;
    char null[1024];
    if(strcmp(ingrpf,"none")==0) {
       i=0;
       ngrp=1;
       init_grp();
       grp[i].natom=natom;
       grp[i].init_atom();
       for(j=0;j<grp[i].natom;j++) grp[i].atm[j]=j;
       cal_grp_mass();
    } else {
       cout<<"Reading group file "<<ingrpf<<endl;
       ifstream inf(ingrpf,ios::in);
       if(!inf.is_open()) {
           cout<<"Error: group file "<<ingrpf<<" cannot be opened."<<endl;
           return filenotopen;
       }
       inf>>ngrp; ngrp++; //last group would be the whole system
       init_grp();
       for(i=0;i<ngrp-1;i++) {
          inf>>grp[i].natom;
          grp[i].init_atom();
          for(j=0;j<grp[i].natom;j++) { inf>>grp[i].atm[j]; grp[i].atm[j]--; }
       }
       i=ngrp-1;
       grp[i].natom=natom;
       grp[i].init_atom();
       for(j=0;j<grp[i].natom;j++) grp[i].atm[j]=j;
       if(inf.eof()) {
          cout<<"Error: group file format error"<<endl;
          return 1;
       }
       cal_grp_mass();
       unsigned long fl = inf.tellg();
       do{inf>>null; } while (strcmp(null,"Constraints")!=0&&!inf.eof());
       if(!inf.eof()) {
          for(i=0;i<ngrp;i++) inf>>grp[i].constraint;
          if(inf.eof()) {
             cout<<"Error: group file constraint format error"<<endl;
             return 1;
          }
       }
       inf.clear();
       inf.seekg(fl,ios::beg);
       do{inf>>null;} while (strcmp(null,"RotationalSymmetryNumber")!=0&&!inf.eof());
       if(!inf.eof()) {
          for(i=0;i<ngrp;i++) inf>>grp[i].rotsym;
          if(inf.eof()) {
             cout<<"Error: group file rotsym format error"<<endl;
             return 1;
          }
       }
       inf.clear();
       inf.seekg(fl,ios::beg);
       do{inf>>null;} while (strcmp(null,"LinearMoleculeFlag")!=0&&!inf.eof());
       if(!inf.eof()) {
          for(i=0;i<ngrp;i++) inf>>grp[i].linear;
          if(inf.eof()) {
             cout<<"Error: group file linear format error"<<endl;
             return 1;
          }
       }
       inf.clear();
       inf.seekg(fl,ios::beg);
       do{inf>>null;} while (strcmp(null,"PartialMolarVolumeRatio")!=0&&!inf.eof());
       if(!inf.eof()) {
          for(i=0;i<ngrp;i++) inf>>grp[i].pmvr;
          if(inf.eof()) {
             cout<<"Error: group file pmvr format error"<<endl;
             return 1;
          }
       }
       inf.clear();
       inf.seekg(fl,ios::beg);
       do{inf>>null;} while (strcmp(null,"AtomSize")!=0&&!inf.eof());
       if(!inf.eof()) {
          for(i=0;i<ngrp;i++) { 
              inf>>null;
              if(strcmp(null,"default")==0) grp[i].atomsize=-1;
              else grp[i].atomsize=atof(null);
          }
          if(inf.eof()) {
             cout<<"Error: group file AtomSize format error"<<endl;
             return 1;
          }
       }
       inf.clear();
       inf.seekg(fl,ios::beg);
       do{inf>>null;} while (strcmp(null,"BondSize")!=0&&!inf.eof());
       if(!inf.eof()) {
          for(i=0;i<ngrp;i++) {
              inf>>null;
              if(strcmp(null,"default")==0) grp[i].bondsize=-1;
              else grp[i].bondsize=atof(null);
          }
          if(inf.eof()) {
             cout<<"Error: group file BondSize format error"<<endl;
             return 1;
          }
       }
       inf.clear();
       inf.seekg(fl,ios::beg);
       do{inf>>null;} while (strcmp(null,"LabelSize")!=0&&!inf.eof());
       if(!inf.eof()) {
          for(i=0;i<ngrp;i++) {
              inf>>null;
              if(strcmp(null,"default")==0) grp[i].labelsize=-1;
              else grp[i].labelsize=atof(null);
          }
          if(inf.eof()) {
             cout<<"Error: group file LabelSize format error"<<endl;
             return 1;
          }
       }
       inf.clear();
       inf.seekg(fl,ios::beg);
       do{inf>>null;} while (strcmp(null,"Color")!=0&&!inf.eof());
       if(!inf.eof()) {
          for(i=0;i<ngrp;i++) {
              inf>>null;
              if(strcmp(null,"default")==0) grp[i].color=-1;
              else grp[i].color=atoi(null);
          }
          if(inf.eof()) {
             cout<<"Error: group file Color format error"<<endl;
             return 1;
          }
       }
       inf.close();
       cout<<"Finish reading group file "<<ingrpf<<endl;
    }

    return 0;
}

int MODEL::init_grp()
{
    if(ngrp==0) return 1;
    if(grp!=NULL) delete [] grp;
    grp=new GROUP [ngrp];
    return 0;
}

int MODEL::init_mol()
{
    if(nmol==0) return 1;
    if(mol!=NULL) delete [] mol;
    mol=new MOLECULE [nmol];
    return 0;
}

int MODEL::init_atom()
{

    if(natom==0) return 1;
    if(atom!=NULL) delete [] atom;

    atom=new ATOM [natom];

    return 0;
}

int MODEL::init_bond()
{
    if(nbond==0) return 1;
    if(bond!=NULL) delete [] bond;
    bond=new BOND [nbond];
    return 0;
}

int MODEL::init_angle()
{
    if(nangle==0) return 1;
    if(angle!=NULL) delete [] angle;
    angle=new ANGLE [nangle];
    return 0;
}

int MODEL::init_tor()
{
    if(ntor==0) return 1;
    if(tor!=NULL) delete [] tor;
    tor=new TORSION [ntor];
    return 0;
}

int MODEL::init_inv()
{
    if(ninv==0) return 1;
    if(inv!=NULL) delete [] inv;
    inv=new INVERSION [ninv];
    return 0;
}

int MODEL::init_imp()
{
    if(nimp==0) return 1;
    if(imp!=NULL) delete [] imp;
    imp=new IMPROPER [nimp];
    return 0;
}

int MODEL::mass2types()
{
    int i;
    double m;
    for(i=0;i<natom;i++) {
       if(strcmp(atom[i].name,"X")==0) { //find element name from its mass
          m=atom[i].mass;
          if(m>0.9&&m<1.1) strcpy(atom[i].name,"H"); 
          else if(m>1.9&&m<2.1) strcpy(atom[i].name,"He");
          else if(m>11.9&&m<12.1) strcpy(atom[i].name,"C");
          else if(m>13.9&&m<14.1) strcpy(atom[i].name,"N");
          else if(m>15.9&&m<16.1) strcpy(atom[i].name,"O");
          else if(m>18.9&&m<19.1) strcpy(atom[i].name,"F");
          else if(m>22.9&&m<23.1) strcpy(atom[i].name,"Na");
          else if(m>31.9&&m<32.1) strcpy(atom[i].name,"S");
          else if(m>35.4&&m<35.6) strcpy(atom[i].name,"Cl");
          else if(m>194.9&&m<195.1) strcpy(atom[i].name,"Pt");
          else if(m>39.9&&m<40.0) strcpy(atom[i].name,"Ar");
       }
       if(strcmp(atom[i].fftype,"X_")==0) { //assign force field type
          sprintf(atom[i].fftype,"%s_",atom[i].name);
       }
    }
    return 0;
}

int MODEL::element2prp()
{
    int i;
    double m;
    for(i=0;i<natom;i++) {
          if(     strcmp(atom[i].name,"H")==0)    {
              atom[i].color=WHITE; 
              atom[i].radius=1.20;
              atom[i].crad=0.37;
              atom[i].mass=1.008;
          } else if(strcmp(atom[i].name,"He")==0) { 
              atom[i].color=DEEPPINK; 
              atom[i].radius=1.40;
              atom[i].crad=0.00;
              atom[i].mass=4.003;
          } else if(strcmp(atom[i].name,"C")==0)  { 
              atom[i].color=GREY; 
              atom[i].radius=1.70;
              atom[i].crad=0.77;
              atom[i].mass=12.01;
          } else if(strcmp(atom[i].name,"N")==0)  { 
              atom[i].color=LIGHTBLUE;
              atom[i].radius=1.55;
              atom[i].crad=0.74;
              atom[i].mass=14.01;
          } else if(strcmp(atom[i].name,"O")==0)  { 
              atom[i].color=RED;
              atom[i].radius=1.52;
              atom[i].crad=0.74;
              atom[i].mass=16.00;
          } else if(strcmp(atom[i].name,"F")==0)  { 
              atom[i].color=GREEN;
              atom[i].radius=1.47;
              atom[i].crad=0.72;
              atom[i].mass=19.00;
          } else if(strcmp(atom[i].name,"Na")==0) { 
              atom[i].color=BLUE;
              atom[i].radius=2.27;
              atom[i].crad=1.66;
              atom[i].mass=22.99;
          } else if(strcmp(atom[i].name,"S")==0)  { 
              atom[i].color=YELLOW;
              atom[i].radius=1.80;
              atom[i].crad=1.04;
              atom[i].mass=32.07;
          } else if(strcmp(atom[i].name,"Cl")==0) { 
              atom[i].color=CYAN;
              atom[i].radius=1.75;
              atom[i].crad=0.99;
              atom[i].mass=35.45;
          } else if(strcmp(atom[i].name,"Pt")==0) { 
              atom[i].color=DEEPPINK;
              atom[i].radius=1.75;
              atom[i].crad=1.21;
              atom[i].mass=195.08;
          } else if(strcmp(atom[i].name,"Ar")==0) { 
              atom[i].color=PURPLE;
              atom[i].radius=1.75;
              atom[i].crad=1.21;
              atom[i].mass=39.948;
          } else { 
              atom[i].color=DEEPPINK;
              atom[i].radius=0.9;
              atom[i].crad=0.10;
              atom[i].mass=0.01;
          }
          if(strcmp(atom[i].fftype,"X_")==0) { 
              strcpy(atom[i].fftype,atom[i].name);
              strcat(atom[i].fftype,"_");
          }
    }
    return 0;
}

int MODEL::ff2mass()
{
    int i,j;
    double m;
    for(i=0;i<natom;i++) {
       for(j=0;j<ff.natmt;j++) {
          if(strcmp(atom[i].fftype,ff.atmt[j].fftype)==0) {
             atom[i].mass=ff.atmt[j].mass;
             break;
          }
       }
       if(j>=ff.natmt) cout<<" Warning: Mass assignment for atom "<<i+1<<" "<<atom[i].fftype<<" failed"<<endl;
    }
    return 0;
}

int MODEL::cal_mass()
{
    int i,j;
    prp.mass=0;
    for(i=0;i<natom;i++) {
       prp.mass += atom[i].mass;
    }
    for(i=0;i<nmol;i++) {
       mol[i].cal_mass();
    }
    return 0;
}

int MODEL::cal_grp_mass()
{
    int i,j;
    for(i=0;i<ngrp;i++) {
       grp[i].mass=0;
       for(j=0;j<grp[i].natom;j++) grp[i].mass += atom[grp[i].atm[j]].mass;
    }
    if(0) {
      for(i=0;i<ngrp;i++) printf("group %d mass %f\n",i+1,grp[i].mass);
    }
    return 0;
}

int MODEL::ck_range()
{
    if (periodic==0) {
       int i,k;
       double min[3],max[3];
       for(k=0;k<3;k++) { 
         cell.o[k]=0.0;
         max[k]=-1e99;
         min[k]=1e99;
       }
       for(i=0;i<natom;i++) {
         for(k=0;k<3;k++) {
            if(atom[i].pv[k]>max[k]) max[k]=atom[i].pv[k];
            else if(atom[i].pv[k]<min[k]) min[k]=atom[i].pv[k];
         }
       }
       cell.la=max[0]-min[0]; cell.alpha=90;
       cell.lb=max[1]-min[1]; cell.beta=90;
       cell.lc=max[2]-min[2]; cell.gamma=90;
       cell.labc2H();
       cell.H2others();
    }
    return 0;
}


int MODEL::cal_multipole()
{
    int i,j,k,aid;
    for(i=0;i<ngrp;i++) grp[i].chg=grp[i].mu[0]=grp[i].mu[1]=grp[i].mu[2]=0;
    prp.chg=prp.mu[0]=prp.mu[1]=prp.mu[2]=0;
    for(i=0;i<nmol;i++) {mol[i].chg=mol[i].mu[0]=mol[i].mu[1]=mol[i].mu[2]=0;}

    for(i=0;i<ngrp;i++) {
       for(j=0;j<grp[i].natom;j++) {
          aid=grp[i].atm[j];
          grp[i].chg += atom[aid].chg;
          for(k=0;k<3;k++) {
             grp[i].mu[k] += (atom[aid].chg*atom[aid].pv[k]*eA2debye);
          }
       }
    }

    for(i=0;i<nmol;i++) {
       mol[i].cal_chgmu();
       //for(j=0;j<mol[i].natom;j++) {
       //   aid=mol[i].atm[j];
       //   mol[i].chg += atom[aid].chg;
       //   for(k=0;k<3;k++) {
       //      mol[i].mu[k] += (atom[aid].chg*atom[aid].pv[k]*eA2debye);
       //   }
       //}
    }
 
if(0) {
    for(i=0;i<natom;i++) printf(" atom %d chg %f \n",i+1,atom[i].chg);
    for(i=0;i<ngrp;i++)  printf(" grp  %d chg %f mu %f %f %f |mu| %f\n",i+1,grp[i].chg,grp[i].mu[0],grp[i].mu[1],grp[i].mu[2],sqrt(grp[i].mu[0]*grp[i].mu[0]+grp[i].mu[1]*grp[i].mu[1]+grp[i].mu[2]*grp[i].mu[2]));
    for(i=0;i<nmol;i++)  printf(" mol  %d chg %f mu %f %f %f |mu| %f\n",i+1,mol[i].chg,mol[i].mu[0],mol[i].mu[1],mol[i].mu[2],sqrt(mol[i].mu[0]*mol[i].mu[0]+mol[i].mu[1]*mol[i].mu[1]+mol[i].mu[2]*mol[i].mu[2]));
}
    j=ngrp-1;
    prp.chg=grp[j].chg;
    for(i=0;i<3;i++) prp.mu[i]=grp[j].mu[i];
    prp.dipole=sqrt(prp.mu[0]*prp.mu[0]+prp.mu[1]*prp.mu[1]+prp.mu[2]*prp.mu[2]);
    return 0;
}

int MODEL::cal_vcomp()
{
    int i;
    for(i=0;i<nmol;i++) mol[i].cal_vc();
    return 0;
}

int MODEL::cal_T()
{
    int i,j;
    double v2,df;
    prp.T=prp.Ek=0;
    for(i=0;i<natom;i++) {
       v2=0;
       for(j=0;j<3;j++) v2+= (atom[i].vel[j]*atom[i].vel[j]);
       prp.Ek += (atom[i].mass*v2);
    }
    df=3.0*natom-(periodic>0?3:0);

    prp.T = prp.Ek/(df*R*0.1);  //0.1 is due to unit conversion
    prp.Ek /= (2*caltoj*100);  //in kcal/mol
    return 0;
}

int MODEL::cal_P()
{
    int i;
    double df;
    double unitc=1000*caltoj*1e21/Na; //conver from kcal/mol-A3 to GPa
    for(i=0;i<6;i++) prp.strs[i]*=(unitc/prp.V); //in GPa
    df=3.0*natom-(periodic>0?3:0);
    prp.P=(df*kb*prp.T)/(3*prp.V*1e-21)+(prp.strs[0]+prp.strs[1]+prp.strs[2])/3; //in GPa

    //printf("%10.4f %10.4f %10.4f\n",(df*kb*prp.T)/(3*prp.V*1e-21),(prp.strs[0]+prp.strs[1]+prp.strs[2])/3,prp.P);
    return 0;
}
 
int MODEL::find_cm()
{
    int i,k;
    prp.cmpv[0]=prp.cmpv[1]=prp.cmpv[2]=0;
    for(i=0;i<natom;i++) {
       for(k=0;k<3;k++) prp.cmpv[k] += atom[i].mass*atom[i].pv[k];
    }
    for(k=0;k<3;k++) prp.cmpv[k]/=prp.mass;
    return 0;
}

int MODEL::find_molcm()
{
    int i,j,k,id;
    for(i=0;i<nmol;i++) {
       mol[i].find_cm();
       //mol[i].pv[0]=mol[i].pv[1]=mol[i].pv[2]=0;
       //for(j=0;j<mol[i].natom;j++) {
       //   id=mol[i].atm[j];
       //   for(k=0;k<3;k++) mol[i].pv[k] += atom[id].mass*atom[id].pv[k];
       //}
       //for(k=0;k<3;k++) mol[i].pv[k]/=mol[i].mass;
       //printf("mol %d cm %10.4f %10.4f %10.4f\n",i+1,mol[i].pv[0],mol[i].pv[1],mol[i].pv[2]);
    }
    return 0;
}

int MODEL::find_molingrp()
{
    int i,j,k,l,match;
    for(i=0;i<ngrp;i++) {
        grp[i].nmol=0;
        for(j=0;j<nmol;j++) {
           match=0;
           for(k=0;k<mol[j].natom;k++) {
              for(l=0;l<grp[i].natom;l++) if(grp[i].atm[l]==mol[j].atm[k]) {match++;break;}
           }
           if(match==mol[j].natom) grp[i].nmol++;
        }
        grp[i].mol=new int [grp[i].nmol];
        grp[i].nmol=0;
        for(j=0;j<nmol;j++) {
           match=0;
           for(k=0;k<mol[j].natom;k++) {
              for(l=0;l<grp[i].natom;l++) if(grp[i].atm[l]==mol[j].atm[k]) {match++;break;}
           }
           if(match==mol[j].natom) grp[i].mol[ grp[i].nmol++ ]=j;
        }
    }

    return 0;
}

int MODEL::superlattice()
{
    int i,k;
    for(i=0;i<natom;i++) for(k=0;k<3;k++) { atom[i].pv[k]+=atom[i].shift[k]; atom[i].shift[k]=0; }
    return 0;
}

int MODEL::mapatom(int type)
{
    int i,j,id;
    double *inpv,*shift;
    switch(type) {
        case MAPDEFAULT: //ensures the cm of each molecule in the unit cell
             find_molcm();
             for(i=0;i<nmol;i++) {
                inpv=&(mol[i].pv[0]);
                for(j=0;j<mol[i].natom;j++) {
                   id=mol[i].atm[j];
                   shift=&(atom[id].shift[0]);
                   cell.Map2UnitCell(inpv,shift);
                }
             }
             break;
        case MAPONECELL:
             for(i=0;i<natom;i++) {
                inpv=&(atom[i].pv[0]);
                shift=&(atom[i].shift[0]);
                cell.Map2UnitCell(inpv,shift);
                //printf("x %f %f y %f %f z %f %f\n",shift[0],atom[i].shift[0],shift[1],atom[i].shift[1],shift[2],atom[i].shift[2]);
             }
             break;
        case MAPORIGINAL:
        default:
             for(i=0;i<natom;i++) {
                atom[i].shift[0]=atom[i].shift[1]=atom[i].shift[2]=0;
             }
             break;
    }
    return 0;
}

int MODEL::init_trj()
{
    ctl.in_trj_flag=trj.init_trj(
      ctl.in_trj_n,ctl.in_trj,ctl.in_trj_flag,ctl.out_trj_flag,
      &cell,atom,&prp
      );
    //check if the model and trj agree
    int i;
    for(i=0;i<ctl.in_trj_n;i++) {
       if(natom!=trj.strj[i].header.nmovatm) {
          cout<<"Error Mismatch of atoms in model and trj # "<<i+1<<": "<<natom<<" "<<trj.strj[i].header.nmovatm<<endl;
          return 1;
       }
    }
    if(ctl.ana_flag) {
      if(ctl.ana_fframe==0||ctl.ana_fframe>trj.totframe) ctl.ana_fframe=trj.totframe;
      if(ctl.ana_iframe< 0||ctl.ana_iframe>trj.totframe) ctl.ana_iframe=1;
      if(ctl.ana_sframe<=0) ctl.ana_sframe=1;
    }

    if(ctl.visu_flag) {
      if(ctl.visu_fframe==0||ctl.visu_fframe>trj.totframe) ctl.visu_fframe=trj.totframe;
      if(ctl.visu_iframe< 0||ctl.visu_iframe>trj.totframe) ctl.visu_iframe=1;
      if(ctl.visu_sframe==0) ctl.visu_sframe=1;
    }

    if(ctl.in_trj_flag==cpmdtrj) {
       for(i=0;i<ctl.in_trj_n;i++) {
          trj.strj[i].header.timestep=timestep;
       }
    }
    return 0;
}

int MODEL::ck_exocyclic(int iatom1,int iatom2)
{
    //determines whet;her bond between iatom1 and iatom2 is in or near an "aromatic" ring
    //returns 1 for an endocyclic bond (modified from arringbond.f of Biograf)
    //        2 for a bond not involved with a ring
    //        3 for an exocyclic bond

    int MAXDEPTH=10; //maximum distance to search
    int i,k,ck1_sp2,ck2_sp2;
    int ring;
    int ringsize,depth,atm,iatm1,iatm2,path,next;

    int atomsave[MAXDEPTH+1];
    int pathsave[MAXDEPTH+1];
    int used[natom];
    int retry;

    iatm1 = iatom1;
    iatm2 = iatom2;
    retry = true;

    k=iatm1;
    ck1_sp2=(strncmp(atom[k].fftype,"C_R",3)==0 || strncmp(atom[k].fftype,"O_2",3)==0 || strncmp(atom[k].fftype,"N_R",3)==0);
    k=iatm2;
    ck2_sp2=(strncmp(atom[k].fftype,"C_R",3)==0 || strncmp(atom[k].fftype,"O_2",3)==0 || strncmp(atom[k].fftype,"N_R",3)==0);

    if(!ck1_sp2 || !ck2_sp2 ) return 2; //bond not involved with an aromatic ring
 
    //USED will keep track of atoms we've been to 

    while (true) { 
      for (i=0;i<natom;i++) used[i]=false;
 
      ring=false;
      ringsize=999999;
      depth=0;
      atm=iatm1;
      used[atm]=true;
      used[iatm2]=true;
      path=0;

      while (true) {

         if (path>=atom[atm].ncnt) next=-1;
         else next=atom[atm].cnt[path];
         path=path+1;
 
         if (next==-1) { // exhausted connections (paths) at this level,back up 
  
            if (depth==0) {
               if(ring) return 3;   //exocyclic bond
               else if(retry) {
                  retry =false;
                  iatm1 =iatom2;
                  iatm2 =iatom1;
                  break;            //get out of the inner while loop
               } else
                  return 2;         //bond not involved with a ring
            }

            used[atm]=false;
            depth=depth-1;
            atm=atomsave[depth];
            path=pathsave[depth];
            continue;
         }
 
         //don't allow us to go straight back to previous atom 
 
         if (depth>0) 
            if (next==atomsave[depth-1]) continue;
 
         //see if we are back to the starting atom
 
         if (next==iatm1) {
            ring=true;
            if(ringsize>(depth+1)) ringsize=depth+1;
            //ringsize=min(ringsize,depth+1);
            continue;
         }
 
         //check if we are at iatm2
         if( next==iatm2 && depth>0 )  return 1; //endocyclic bond
 
         //don't allow us to reuse an atom (this stops side chains from being added
         //also skip things that are not "aromatic"
         k=next;
         ck1_sp2=(strncmp(atom[k].fftype,"C_R",3)==0 || strncmp(atom[k].fftype,"O_2",3)==0 || strncmp(atom[k].fftype,"N_R",3)==0);
         if (used[next] || !ck1_sp2) continue;

         //check we are not too deep
         if (depth>MAXDEPTH) continue;
 
         //push everything on the stack and continue
 
         used[next]=true;
         atomsave[depth]=atm;
         pathsave[depth]=path;
         depth=depth+1;
         atm=next;
         path=0;
      }
    }
    return 0;
}

int MODEL::out_structure(ostream *outf)
{
    char null[1024];
    int i;
    *outf<<endl<<"Atom positions"<<endl;
    for(i=0;i<natom;i++) {
       sprintf(null,"%5d %-5s %5d %10.4f %10.4f %10.4f",i+1,atom[i].fftype,atom[i].id+1,atom[i].pv[0],atom[i].pv[1],atom[i].pv[2]);
       *outf<<null<<endl;
    }

    *outf<<endl<<"Bond lengths"<<endl;
    for(i=0;i<nbond;i++) {
       bond[i].cal_len();
       sprintf(null,"%5d %-5s %-5s %5d %5d %10.4f",i+1,atom[bond[i].atm[0]].fftype,atom[bond[i].atm[1]].fftype,atom[bond[i].atm[0]].id+1,atom[bond[i].atm[1]].id+1,bond[i].len);
       *outf<<null<<endl;
    }

    *outf<<endl<<"Angles"<<endl;
    for(i=0;i<nangle;i++) {
       angle[i].cal_theta();
       sprintf(null,"%5d %-5s %-5s %-5s %5d %5d %5d %10.4f",i+1,atom[angle[i].atm[0]].fftype,atom[angle[i].atm[1]].fftype,atom[angle[i].atm[2]].fftype,atom[angle[i].atm[0]].id+1,atom[angle[i].atm[1]].id+1,atom[angle[i].atm[2]].id+1,angle[i].theta/PI*180.0);
       *outf<<null<<endl;
    }

    *outf<<endl<<"Dihedral Torsions"<<endl;
    for(i=0;i<ntor;i++) {
       tor[i].cal_theta();
       sprintf(null,"%5d %-5s %-5s %-5s %-5s %5d %5d %5d %5d %10.4f",i+1,atom[tor[i].atm[0]].fftype,atom[tor[i].atm[1]].fftype,atom[tor[i].atm[2]].fftype,atom[tor[i].atm[3]].fftype,atom[tor[i].atm[0]].id+1,atom[tor[i].atm[1]].id+1,atom[tor[i].atm[2]].id+1,atom[tor[i].atm[3]].id+1,tor[i].theta/PI*180.0);
       *outf<<null<<endl;
    }

    *outf<<endl<<"Inversions"<<endl;
    for(i=0;i<ninv;i++) {
       inv[i].cal_theta();
       sprintf(null,"%5d %-5s %-5s %-5s %-5s %5d %5d %5d %5d %10.4f",i+1,atom[inv[i].atm[0]].fftype,atom[inv[i].atm[1]].fftype,atom[inv[i].atm[2]].fftype,atom[inv[i].atm[3]].fftype,atom[inv[i].atm[0]].id+1,atom[inv[i].atm[1]].id+1,atom[inv[i].atm[2]].id+1,atom[inv[i].atm[3]].id+1,inv[i].theta/PI*180.0);
       *outf<<null<<endl;
    }

    return 0;
}

