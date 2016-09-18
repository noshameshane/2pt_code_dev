#include "trajectory.h"

int TRJheader::init_header()
{
    if(nmovatm==0) return 1;
    if(mvatmofst!=NULL) delete [] mvatmofst;
    mvatmofst=new int [nmovatm];
    int i;
    for(i=0;i<nmovatm;i++) { mvatmofst[i]=i+1; }
    mvatmpfu[0]=nmovatm;   //number of movable atoms
    natmpfu[0]=nmovatm;    //total number of atoms

    return 0;
}

int TRJheader::CleanTRJheader()
{
	int i,j;
	for(i=0;i<5;i++) hdr[i]='\0';  
	for(i=0;i<20;i++) icntrl[i]=0; 
	ntrjti=neexti=0;
	for(i=0;i<10;i++) 
	{ 
		for(j=0;j<80;j++) trjtic[i][j]=eextic[i][j]='\0'; 
	}
	period=molxtl=lcanon=defcel=prtthrm=lnose=lnpecan=ltmpdamp=nflusd=0;
	for(i=0;i<MAXFILE;i++) { mvatmpfu[i]=natmpfu[i]=0; for(j=0;j<8;j++) decusd[i][j]='\0'; }
	totmov=0;
	//for(i=0;i<nmovatm;i++) { mvatmofst[i]=0; }
	leexti=lparti=natom=nmovatm=movatm1=movatmn=version=0;
	for(i=0;i<80;i++) eextit[i]=partit[i]='\0';
        xyzfreq=velfreq=engfreq=strfreq=0;
        nxyz=nvel=neng=nstr=0;
        timestep=0;
        totaccustep=0;

        //stlin 12/07/2007
        strcpy(hdr,"MDTR");
        version=icntrl[0]=2010;
        ntrjti=1; strcpy(trjtic[0],"COMET trajectory");
        nflusd=1;
        strcpy(decusd[0],"TEST");
        leexti=7;
        strcpy(eextit,"NOTITLE");
        lparti=5;
        strcpy(partit,"NOPAR");
     
        return 0;
}

int TRJheader::ReadBinHeader(ifstream *intrj)
{
        int i;
        float null;  //dummy variable, one carriage return
        double crtn;

        CleanTRJheader();
	
//     ----- Read header information -----
        intrj->read((char *)&null,sizeof(null)); 
        intrj->read((char *)&hdr,sizeof(hdr)-1);

	for(i=0;i<20;i++) intrj->read((char *)&icntrl[i],sizeof(icntrl[i])); 
        intrj->read((char *)&crtn,sizeof(crtn)); //read carriage return produced in unformatted Fortran file

	version=icntrl[0];

	if (version>=311) {
         intrj->read((char *)&ntrjti,sizeof(ntrjti));
         for (i=0;i<ntrjti;i++)	intrj->read((char *)&trjtic[i],sizeof(trjtic[i]));
	}
        intrj->read((char *)&crtn,sizeof(crtn));

	intrj->read((char *)&neexti,sizeof(neexti));
	for (i=0;i<neexti;i++) intrj->read((char *)&eextic[i],sizeof(eextic[i]));
        intrj->read((char *)&crtn,sizeof(crtn));

        if (version<=150) {
         period=false;
         molxtl=false;
         lcanon=false;
         defcel=false;
         prtthrm=false;
	} else if (version<300) {
         intrj->read((char *)&period,sizeof(period));
	 intrj->read((char *)&molxtl,sizeof(molxtl));
	 intrj->read((char *)&lcanon,sizeof(lcanon));
	 intrj->read((char *)&defcel,sizeof(defcel));
	 intrj->read((char *)&prtthrm,sizeof(prtthrm));
         intrj->read((char *)&crtn,sizeof(crtn));
         lnose =false;
         lnpecan=false;
         ltmpdamp=false;
	} else {
         intrj->read((char *)&period,sizeof(period)); //Periodicity
	 intrj->read((char *)&molxtl,sizeof(molxtl)); //MolXtl
	 intrj->read((char *)&lcanon,sizeof(lcanon)); //Canonical
	 intrj->read((char *)&defcel,sizeof(defcel)); //DefCel
	 intrj->read((char *)&prtthrm,sizeof(prtthrm)); //PertTheory
	 intrj->read((char *)&lnose,sizeof(lnose)); //NoseOrHoover
	 intrj->read((char *)&lnpecan,sizeof(lnpecan)); //NpTCanon
	 intrj->read((char *)&ltmpdamp,sizeof(ltmpdamp)); //TempDamping
         intrj->read((char *)&crtn,sizeof(crtn));
	}

	if (version>=200)
	{
        intrj->read((char *)&nflusd,sizeof(nflusd));
	for(i=0;i<nflusd;i++) intrj->read((char *)&mvatmpfu[i],sizeof(mvatmpfu[i]));
	for(i=0;i<nflusd;i++) intrj->read((char *)&natmpfu[i],sizeof(natmpfu[i]));
	for(i=0;i<nflusd;i++) intrj->read((char *)&decusd[i],sizeof(decusd[i]));
        intrj->read((char *)&crtn,sizeof(crtn));
       
        intrj->read((char *)&totmov,sizeof(totmov));
        nmovatm = totmov;
        init_header();

	for(i=0;i<totmov;i++) intrj->read((char *)&mvatmofst[i],sizeof(mvatmofst[i]));
        intrj->read((char *)&crtn,sizeof(crtn));

	}   else
	{
        intrj->read((char *)&natom,sizeof(natom));
 	intrj->read((char *)&nmovatm,sizeof(nmovatm));
	intrj->read((char *)&movatm1,sizeof(movatm1));
	intrj->read((char *)&movatmn,sizeof(movatmn));
        intrj->read((char *)&crtn,sizeof(crtn));
	}

	intrj->read((char *)&leexti,sizeof(leexti));
	for(i=0;i<leexti;i++) intrj->read((char *)&eextit[i],sizeof(eextit[i]));
        intrj->read((char *)&crtn,sizeof(crtn));

        intrj->read((char *)&lparti,sizeof(lparti));
	for(i=0;i<lparti;i++) intrj->read((char *)&partit[i],sizeof(partit[i]));
        intrj->read((char *)&crtn,sizeof(crtn));

//     ----- end of header information -----
//	intrj->read((char *)&null,sizeof(null));
        return 0;
}

int TRJheader::ReadAscHeader(ifstream *intrj)
{
        int i; 
        char check[1024];
        double null;  //dummy variable

        CleanTRJheader();

//     ----- Read header information -----
        *intrj>>check>>hdr;
        if(strcmp(check,"MPSim-trj_ASCII_file_for_BINARY_conversion")!=0) {
         cout<<" ASCII file does not have the right format. Reading fails."<<endl;
         return 1;
        }
        for(i=0;i<20;i++) *intrj>>icntrl[i];
        version=icntrl[0];

        if (version>=311) {
                *intrj>>ntrjti;
                for (i=0;i<ntrjti;i++) { 
                  intrj->getline(trjtic[i],sizeof(trjtic[i]));
                  intrj->clear();
                }
        }
        *intrj>>neexti;
        for (i=0;i<neexti;i++) {
              intrj->getline(eextic[i],sizeof(eextic[i]));
              intrj->clear();
        }
        if(version<=150) ;
        else if (version<300) {
                *intrj>>period>>molxtl>>lcanon>>defcel>>prtthrm;
        } else {
                *intrj>>period>>molxtl>>lcanon>>defcel>>prtthrm>>lnose>>lnpecan>>ltmpdamp;
        }
        //printf("stlin period %d molxtl %d lcanon %d defcel %d prtthrm %d lnose %d lnpecan %d ltmpdamp %d \n",period,molxtl,lcanon,defcel,prtthrm,lnose,lnpecan,ltmpdamp);
        if (version>=200) {
            *intrj>>nflusd;
            for(i=0;i<nflusd;i++) *intrj>>mvatmpfu[i];
            for(i=0;i<nflusd;i++) *intrj>>natmpfu[i];
            for(i=0;i<nflusd;i++) *intrj>>decusd[i];
            *intrj>>totmov;
            nmovatm = totmov;
            init_header();
            for(i=0;i<totmov;i++) *intrj>>mvatmofst[i];
        } else {
                *intrj>>natom>>nmovatm>>movatm1>>movatmn;
        }

        *intrj>>leexti;
        for(i=0;i<leexti;i++) *intrj>>eextit[i];
        *intrj>>lparti;
        for(i=0;i<lparti;i++) *intrj>>partit[i];

//     ----- end of header information -----
    return 0;
}

int TRJheader::ReadUSRHeader(ifstream *intrj)
{
    CleanTRJheader();
    intrj->read((char *)&totmov,sizeof(totmov));
    nmovatm = totmov;
    init_header();
    cout<<"stlin usrtrj natom "<<nmovatm<<endl;

    return 0;
}

int TRJheader::ReadCPMDHeader(ifstream *intrj,ifstream *ineng,ifstream *instr)
{
    int i,fstp,k;
    char null[1024];

    CleanTRJheader();

    unsigned long tmp;
    tmp=intrj->tellg();
    intrj->clear();
    intrj->seekg(0,ios::beg);
    //find number of atoms
    *intrj>>fstp; i=fstp;
    totmov=0;
    while(i==fstp) {
      totmov++;
      intrj->getline(null,1024);
      *intrj>>i;
    }
    nmovatm = totmov;  //number of movable atoms
    init_header();

    //stlin 12/07/2007
    ntrjti=1; strcpy(trjtic[0],"CPMD trajectory");
    period=3;

    //determine xyzfreq, velfreq
    xyzfreq=velfreq=i-fstp;
    //cout<<xyzfreq<<endl;
    intrj->seekg(tmp, ios::beg);

    //determine engfreq
    tmp=ineng->tellg();
    ineng->clear();
    ineng->seekg(0,ios::beg);
    *ineng>>fstp; ineng->getline(null,1024); *ineng>>i;
    engfreq=i-fstp;
    ineng->seekg(tmp, ios::beg);

    //determine strfreq
    tmp=instr->tellg();
    instr->clear();
    instr->seekg(0,ios::beg);
    do{*instr>>null;} while(strcmp(null,"Step:")!=0);
    *instr>>fstp;
    do{*instr>>null;} while(strcmp(null,"Step:")!=0);
    *instr>>i;
    strfreq=i-fstp;
    instr->seekg(tmp, ios::beg);

    printf("xyzfreq %d velfreq %d engfreq %d strfreq %d\n",xyzfreq,velfreq,engfreq,strfreq);
    //determine the skips of each trj, not fully optimized yet
    if(xyzfreq%strfreq) {
      nxyz=strfreq;
      nstr=xyzfreq;
      neng=xyzfreq*strfreq;
    } else {
      nxyz=nvel=1;
      nstr=xyzfreq/strfreq;
      neng=xyzfreq/engfreq;
    }
    printf("nxyz %d nvel %d neng %d nstr %d\n",nxyz,nvel,neng,nstr);
    return 0;
}

int TRJheader::ReadLMPHeader(ifstream *intrj,ifstream *ineng)
{
    int i,fstp,k;
    char null[1024];

    CleanTRJheader();

    unsigned long tmp;
    tmp=intrj->tellg();
    intrj->clear();
    intrj->seekg(0,ios::beg);
    //find number of atoms
    *intrj>>fstp; i=fstp;
    totmov=0;
    while(i==fstp) {
      totmov++;
      intrj->getline(null,1024);
      *intrj>>i;
    }
    nmovatm = totmov;  //number of movable atoms
    init_header();

    //stlin 12/07/2007
    ntrjti=1; strcpy(trjtic[0],"CPMD trajectory");
    period=3;

    //determine xyzfreq, velfreq
    xyzfreq=velfreq=i-fstp;
    //cout<<xyzfreq<<endl;
    intrj->seekg(tmp, ios::beg);

    //determine engfreq
    tmp=ineng->tellg();
    ineng->clear();
    ineng->seekg(0,ios::beg);
    *ineng>>fstp; ineng->getline(null,1024); *ineng>>i;
    engfreq=i-fstp;
    ineng->seekg(tmp, ios::beg);

    printf("xyzfreq %d velfreq %d engfreq %d strfreq %d\n",xyzfreq,velfreq,engfreq,strfreq);
    //determine the skips of each trj, not fully optimized yet
    if(xyzfreq%strfreq) {
      nxyz=strfreq;
      nstr=xyzfreq;
      neng=xyzfreq*strfreq;
    } else {
      nxyz=nvel=1;
      nstr=xyzfreq/strfreq;
      neng=xyzfreq/engfreq;
    }
    printf("nxyz %d nvel %d neng %d nstr %d\n",nxyz,nvel,neng,nstr);
    return 0;
}

int TRJheader::WriteASCIIheader(ostream *outf, int oformat, int ndel, int iddel[])
{

	int i;
	switch(oformat)
	{
	case 0:
		*outf<<"MPSim-trj_ASCII_file_for_BINARY_conversion"<<endl;
		*outf<<hdr<<" ";
		for(i=0;i<20;i++) *outf<<icntrl[i]<<" ";

		if (version>=311) 
		{
			*outf<<endl<<ntrjti<<" ";
			for (i=0;i<ntrjti;i++) *outf<<trjtic[i]<<endl;
		}

		*outf<<neexti<<" ";
		for (i=0;i<neexti;i++) *outf<<eextic[i]<<endl;

		if(version<=150) ;
		else if (version<300) 
		{
				*outf<<endl<<period<<" "<<molxtl<<" "<<lcanon<<" "
					<<defcel<<" "<<prtthrm<<endl;
		} else 
		{
			*outf<<endl<<period<<" "<<molxtl<<" "<<lcanon<<" "<<defcel<<" "
				<<prtthrm<<" "<<lnose<<" "<<lnpecan<<" "<<ltmpdamp<<endl;
		}
		if (version>=200)
		{
			*outf<<nflusd<<" ";
			for(i=0;i<nflusd;i++) *outf<<mvatmpfu[i]<<" ";
			for(i=0;i<nflusd;i++) *outf<<natmpfu[i]<<" ";
			for(i=0;i<nflusd;i++) *outf<<decusd[i]<<" ";
			*outf<<endl<<totmov<<" ";
			for(i=0;i<totmov;i++) *outf<<mvatmofst[i]<<" ";
		} else 
		{
			*outf<<natom<<" "<<nmovatm<<" "<<movatm1<<" "<<movatmn<<" ";
		}
		*outf<<endl<<leexti<<" ";
		for(i=0;i<leexti;i++) *outf<<eextit[i]<<" ";
		*outf<<endl<<lparti<<" ";
		for(i=0;i<lparti;i++) *outf<<partit[i]<<" ";
		*outf<<endl; 
		break;
	case 1:
        *outf<<" Header: "<<hdr<<endl;
        *outf<<" Version: "<<version<<endl;
        *outf<<" NRemarks: "<<ntrjti<<endl;
        *outf<<" Neexti: "<<neexti<<endl;
        *outf<<" Flags: "<<period<<" "<<molxtl<<" "<<lcanon<<" "<<defcel
            <<" "<<prtthrm<<" "<<lnose<<" "<<lnpecan<<" "<<ltmpdamp<<endl;
        *outf<<" Nfileused: "<<nflusd<<" "<<mvatmpfu[0]<<" "<<natmpfu[0]
            <<" "<<decusd[0]<<endl;
        *outf<<" Totmovatm: "<<totmov<<" "<<nmovatm<<endl;
        *outf<<" Leexti: "<<leexti<<" "<<eextit<<endl;
        *outf<<" Parti: "<<lparti<<" "<<partit<<endl<<endl;
		*outf<<" Time(ps)  Temp(K) Etot(kcal/mol) Epot(kcal/mol) Eele(kcal/mol) Evdw(kcal/mol) Eval(kcal/mol) Ehb(kcal/mol)"<<endl;

		*outf<<"-----------------------------------------------------------------------------------------------------------"<<endl;
		break;
	case 2:
        	*outf<<"totpts= "<<100<<endl;
        	*outf<<"total area= "<<0<<endl;
        	*outf<<"   ID    X        Y          Z         AREA    TYPE CHARGE    POTENTIAL A1   A2  A3    Nx        Ny        Nz    "<<endl;
		break;
	case 3: /* removing certain number of atoms */
		*outf<<"MPSim-trj_ASCII_file_for_BINARY_conversion"<<endl;
		*outf<<hdr<<" ";
		for(i=0;i<20;i++) *outf<<icntrl[i]<<" ";

		if (version>=311) 
		{
			*outf<<endl<<ntrjti<<" ";
			for (i=0;i<ntrjti;i++) *outf<<trjtic[i]<<endl;
		}

		*outf<<neexti<<" ";
		for (i=0;i<neexti;i++) *outf<<eextic[i]<<endl;

		if(version<=150) ;
		else if (version<300) 
		{
				*outf<<endl<<period<<" "<<molxtl<<" "<<lcanon<<" "
					<<defcel<<" "<<prtthrm<<endl;
		} else 
		{
			*outf<<endl<<period<<" "<<molxtl<<" "<<lcanon<<" "<<defcel<<" "
				<<prtthrm<<" "<<lnose<<" "<<lnpecan<<" "<<ltmpdamp<<endl;
		}
		if (version>=200)
		{
			*outf<<nflusd<<" ";
			for(i=0;i<nflusd;i++) *outf<<mvatmpfu[i]-ndel<<" ";
			for(i=0;i<nflusd;i++) *outf<<natmpfu[i]-ndel<<" ";
			for(i=0;i<nflusd;i++) *outf<<decusd[i]<<" ";
			*outf<<endl<<totmov-ndel<<" ";
			for(i=0;i<totmov;i++) {
				if(iddel[i]) *outf<<mvatmofst[i]<<" ";
			}
		} else 
		{
			*outf<<natom<<" "<<nmovatm<<" "<<movatm1<<" "<<movatmn<<" ";
		}
		*outf<<endl<<leexti<<" ";
		for(i=0;i<leexti;i++) *outf<<eextit[i]<<" ";
		*outf<<endl<<lparti<<" ";
		for(i=0;i<lparti;i++) *outf<<partit[i]<<" ";
		*outf<<endl; 
		break;
	case 4:
		break;
	case 5:
		break;
	default:
		break;
	}
        return 0;
}

int TRJheader::WriteBinaryheader(ostream *outf, int oformat)
{
	int i;
    int size;	//size of the memory between one line (Fortran write format)

    size=sizeof(hdr)-1+20*sizeof(icntrl[0]); 
	outf->write((char *)&size,sizeof(size));
    outf->write((char *)&hdr,sizeof(hdr)-1);
	for(i=0;i<20;i++) outf->write((char *)&icntrl[i],sizeof(icntrl[i])); 
	outf->write((char *)&size,sizeof(size));


	if (version>=311) 
	{
		size=sizeof(ntrjti)+ntrjti*sizeof(trjtic[0]);
		outf->write((char *)&size,sizeof(size));
        outf->write((char *)&ntrjti,sizeof(ntrjti));
        for (i=0;i<ntrjti;i++)	outf->write((char *)&trjtic[i],sizeof(trjtic[i]));
		outf->write((char *)&size,sizeof(size));
	}

	size=sizeof(neexti)+neexti*sizeof(eextic[0]);
	outf->write((char *)&size,sizeof(size));
	outf->write((char *)&neexti,sizeof(neexti));
	for (i=0;i<neexti;i++) outf->write((char *)&eextic[i],sizeof(eextic[i]));
	outf->write((char *)&size,sizeof(size));

    if (version<=150) ;	
	else if (version<300) 
	{
 	 size=sizeof(period)*5;
 	 outf->write((char *)&size,sizeof(size));
     outf->write((char *)&period,sizeof(period));
	 outf->write((char *)&molxtl,sizeof(molxtl));
	 outf->write((char *)&lcanon,sizeof(lcanon));
	 outf->write((char *)&defcel,sizeof(defcel));
	 outf->write((char *)&prtthrm,sizeof(prtthrm));
 	 outf->write((char *)&size,sizeof(size));
	}   else
	{
 	 size=sizeof(period)*8;
 	 outf->write((char *)&size,sizeof(size));
         outf->write((char *)&period,sizeof(period));
	 outf->write((char *)&molxtl,sizeof(molxtl));
	 outf->write((char *)&lcanon,sizeof(lcanon));
	 outf->write((char *)&defcel,sizeof(defcel));
	 outf->write((char *)&prtthrm,sizeof(prtthrm));
	 outf->write((char *)&lnose,sizeof(lnose));
	 outf->write((char *)&lnpecan,sizeof(lnpecan));
	 outf->write((char *)&ltmpdamp,sizeof(ltmpdamp));
 	 outf->write((char *)&size,sizeof(size));
	}

	if (version>=200)
	{
 	 size=sizeof(nflusd)+nflusd*(sizeof(mvatmpfu[0])+sizeof(natmpfu[0])+sizeof(decusd[0]));
 	 outf->write((char *)&size,sizeof(size));
        outf->write((char *)&nflusd,sizeof(nflusd));
	for(i=0;i<nflusd;i++) outf->write((char *)&mvatmpfu[i],sizeof(mvatmpfu[i]));
	for(i=0;i<nflusd;i++) outf->write((char *)&natmpfu[i],sizeof(natmpfu[i]));
	for(i=0;i<nflusd;i++) outf->write((char *)&decusd[i],sizeof(decusd[i]));
 	 outf->write((char *)&size,sizeof(size));
       
 	 size=sizeof(totmov)+totmov*(sizeof(mvatmofst[0]));
 	 outf->write((char *)&size,sizeof(size));
       outf->write((char *)&totmov,sizeof(totmov));
	for(i=0;i<totmov;i++) outf->write((char *)&mvatmofst[i],sizeof(mvatmofst[i]));
 	 outf->write((char *)&size,sizeof(size));
	}   else
	{
 	 size=sizeof(natom)+sizeof(nmovatm)+sizeof(movatm1)+sizeof(movatmn);
 	 outf->write((char *)&size,sizeof(size));
        outf->write((char *)&natom,sizeof(natom));
 	outf->write((char *)&nmovatm,sizeof(nmovatm));
	outf->write((char *)&movatm1,sizeof(movatm1));
	outf->write((char *)&movatmn,sizeof(movatmn));
 	 outf->write((char *)&size,sizeof(size));
	}

 	 size=sizeof(leexti)+leexti*sizeof(eextit[0]);
 	 outf->write((char *)&size,sizeof(size));
	outf->write((char *)&leexti,sizeof(leexti));
	for(i=0;i<leexti;i++) outf->write((char *)&eextit[i],sizeof(eextit[i]));
 	 outf->write((char *)&size,sizeof(size));

 	 size=sizeof(lparti)+lparti*sizeof(partit[0]);
 	 outf->write((char *)&size,sizeof(size));
	outf->write((char *)&lparti,sizeof(lparti));
	for(i=0;i<lparti;i++) outf->write((char *)&partit[i],sizeof(partit[i]));
 	 outf->write((char *)&size,sizeof(size));
        return 0;
}

int TRJ::init_contentDouble()
{
    
    if(header.nmovatm==0) return 1;
    //if(x!=NULL) delete [] x;
    //if(y!=NULL) delete [] y;
    //if(z!=NULL) delete [] z;
    //if(velx!=NULL) delete [] velx;
    //if(vely!=NULL) delete [] vely;
    //if(velz!=NULL) delete [] velz;

    //x=new double [header.nmovatm];
    //y=new double [header.nmovatm];
    //z=new double [header.nmovatm];
    //velx=new double [header.nmovatm];
    //vely=new double [header.nmovatm];
    //velz=new double [header.nmovatm];

    return 0;
}

int TRJ::CleanTRJcontentDouble()
{
	int i,j;
	prp->time=0.0;
	prp->step=0;
	prp->T=avetem=dstep=firstt=finalt=0.0;
        prp->Ep=prp->Ebond=prp->Eangle=prp->Etorsion=prp->Einversion
               =prp->Evdw=prp->Eel=prp->Ehb=prp->Et=prp->Ek;
	//e=eb=et=ep=ei=enb=eel=ehb
	ec=eu=tint=tnb=ea=eba
		=eta=epa=eia=enba=eela=ehba
		=eca=eua=tinta=tnba   //=tote=totke
		=totea=tkea=0.0;
	for(i=0;i<12;i++) dum[i]=duma[i]=0.0;
	iconmp=imstep=iconfs=icstep=lvelwr=lfrcwr=lengwr=0;

	prp->P=prp->V=pvtota=pvkina=pvpota=radgyra=0.0;
	pressur =vol =pvtot =pvkin =pvpot =radgyr =0.0;

	signose=zfrict=snose=snoseh=ssdot=sigdyn[0]=sigdyn[1]=0.0;
	zprfrict=qcanon=gamtmp=0.0;
	tcel=ucel=tcela=ucela=0.0;
	//for(i=0;i<6;i++) s2r[i]=s2ra[i]=s2rdot[i]=0.0;
	for(i=0;i<6;i++) s2ra[i]=s2rdot[i]=0.0;
	//for(i=0;i<3;i++) for(j=0;j<3;j++) cell->H[i][j]=0;

	natmcel=0;
	for(i=0;i<6;i++) strsa[i]=prp->strs[i]=0.0;
	extstrs=extstrsa=eabtota=eabvala=eabelha=eabnba=eabmisa=dltfaba=expprta=0.0;
        for(i=0;i<header.nmovatm;i++) //x[i]=y[i]=z[i]=velx[i]=vely[i]=velz[i]=0.0;
          for(j=0;j<3;j++) atom[i].pv[j]=atom[i].vel[j]=0.0;

        return 0;
}

int TRJ::ReadBinEng()
{
        int i;
        double dnull;
        double crtn; //two carriage returns (each carriage return in fortran is 4 bytes)

	intrj.read((char *)&dnull,sizeof(dnull));
			prp->time=dnull; //current time (ps)
			intrj.read((char *)&(prp->step),sizeof(prp->step)); //number of steps
			intrj.read((char *)&(prp->T),sizeof(prp->T)); //instantaneous temperature
			intrj.read((char *)&avetem,sizeof(avetem)); //running averaged temperature
			intrj.read((char *)&dstep,sizeof(dstep)); //time step (ps)
			intrj.read((char *)&firstt,sizeof(firstt)); //initial temperature setting
			intrj.read((char *)&finalt,sizeof(finalt)); //final temperature setting
			intrj.read((char *)&(prp->Ep),sizeof(prp->Ep)); //potential energy
			intrj.read((char *)&(prp->Ebond),sizeof(prp->Ebond)); //bond
			intrj.read((char *)&(prp->Eangle),sizeof(prp->Eangle)); //angle
			intrj.read((char *)&(prp->Etorsion),sizeof(prp->Etorsion)); //torsion
			intrj.read((char *)&(prp->Einversion),sizeof(prp->Einversion)); //inversion
			intrj.read((char *)&(prp->Evdw),sizeof(prp->Evdw)); //vdw
			intrj.read((char *)&(prp->Eel),sizeof(prp->Eel)); //coulobm
			intrj.read((char *)&(prp->Ehb),sizeof(prp->Ehb)); //h-bond
			intrj.read((char *)&ec,sizeof(ec)); //restraint
			for(i=0;i<12;i++) intrj.read((char *)&dum[i],sizeof(dum[i])); //12 cross terms
			intrj.read((char *)&eu,sizeof(eu)); //user
			intrj.read((char *)&tint,sizeof(tint)); //valence 
			intrj.read((char *)&tnb,sizeof(tnb)); //nonbond
			intrj.read((char *)&ea,sizeof(ea));
			intrj.read((char *)&eba,sizeof(eba));
			intrj.read((char *)&eta,sizeof(eta));
			intrj.read((char *)&epa,sizeof(epa));
			intrj.read((char *)&eia,sizeof(eia));
			intrj.read((char *)&enba,sizeof(enba));
			intrj.read((char *)&eela,sizeof(eela));
			intrj.read((char *)&ehba,sizeof(ehba));
			intrj.read((char *)&eca,sizeof(eca));
			intrj.read((char *)&eua,sizeof(eua));
			intrj.read((char *)&tinta,sizeof(tinta));
			for(i=0;i<12;i++) intrj.read((char *)&duma[i],sizeof(duma[i]));
			intrj.read((char *)&tnba,sizeof(tnba));
			intrj.read((char *)&(prp->Et),sizeof(prp->Et));
			intrj.read((char *)&(prp->Ek),sizeof(prp->Ek));
			intrj.read((char *)&totea,sizeof(totea));
			intrj.read((char *)&tkea,sizeof(tkea));
			intrj.read((char *)&iconmp,sizeof(iconmp)); //flag
			intrj.read((char *)&imstep,sizeof(imstep)); //flag
			intrj.read((char *)&lvelwr,sizeof(lvelwr)); //write velocity
			intrj.read((char *)&lfrcwr,sizeof(lfrcwr)); //write force
                        if (header.version>=2012) intrj.read((char *)&lengwr,sizeof(lengwr)); //atomic energy
			intrj.read((char *)&iconfs,sizeof(iconfs)); //flag
			intrj.read((char *)&icstep,sizeof(icstep)); //flag
                        intrj.read((char *)&crtn,sizeof(crtn));

			if (header.version>=155) 
			{
				intrj.read((char *)&pressur,sizeof(pressur));
				intrj.read((char *)&vol,sizeof(vol));
				intrj.read((char *)&pvtot,sizeof(pvtot));
				intrj.read((char *)&pvkin,sizeof(pvkin));
				intrj.read((char *)&pvpot,sizeof(pvpot));
				intrj.read((char *)&radgyr,sizeof(radgyr));

				intrj.read((char *)&(prp->P),sizeof(prp->P));
				intrj.read((char *)&(prp->V),sizeof(prp->V));
				intrj.read((char *)&pvtota,sizeof(pvtota));
				intrj.read((char *)&pvkina,sizeof(pvkina));
				intrj.read((char *)&pvpota,sizeof(pvpota));
				intrj.read((char *)&radgyra,sizeof(radgyra));
                                intrj.read((char *)&crtn,sizeof(crtn));
			}

			if (header.lcanon) 
			{
				  if (header.version<300) 
				  {
					intrj.read((char *)&signose,sizeof(signose));
					intrj.read((char *)&zfrict,sizeof(zfrict));
					intrj.read((char *)&zprfrict,sizeof(zprfrict));
                                        intrj.read((char *)&crtn,sizeof(crtn));
				  } else
				  {
					  if(header.lnose) 
					  {
						  intrj.read((char *)&snose,sizeof(snose));
						  intrj.read((char *)&snoseh,sizeof(snoseh));
						  intrj.read((char *)&ssdot,sizeof(ssdot));
						  intrj.read((char *)&qcanon,sizeof(qcanon));
                                                  intrj.read((char *)&crtn,sizeof(crtn));
					  }	else 
					  {
						  intrj.read((char *)&signose,sizeof(signose));
						  intrj.read((char *)&zfrict,sizeof(zfrict));
						  intrj.read((char *)&zprfrict,sizeof(zprfrict));
						  intrj.read((char *)&qcanon,sizeof(qcanon));
                                                  intrj.read((char *)&crtn,sizeof(crtn));
					  }
				  }
			}
			
			if(header.version>=220)
			{
				if(header.period) 
				{
					intrj.read((char *)&tcel,sizeof(tcel));
					intrj.read((char *)&tcela,sizeof(tcela));
					//for(i=0;i<6;i++) intrj.read((char *)&s2r[i],sizeof(s2r[i]));
					intrj.read((char *)&(cell->H[0][0]),sizeof(cell->H[0][0]));
					intrj.read((char *)&(cell->H[1][1]),sizeof(cell->H[1][1]));
					intrj.read((char *)&(cell->H[2][2]),sizeof(cell->H[2][2]));
					intrj.read((char *)&(cell->H[2][1]),sizeof(cell->H[2][1]));
					intrj.read((char *)&(cell->H[2][0]),sizeof(cell->H[2][0]));
					intrj.read((char *)&(cell->H[1][0]),sizeof(cell->H[1][0]));
					for(i=0;i<6;i++) intrj.read((char *)&s2rdot[i],sizeof(s2rdot[i]));
					intrj.read((char *)&ucel,sizeof(ucel));
					intrj.read((char *)&ucela,sizeof(ucela));
					for(i=0;i<6;i++) intrj.read((char *)&s2ra[i],sizeof(s2ra[i]));
                                        intrj.read((char *)&crtn,sizeof(crtn));
				}
/*stlin 8/27/2002*/
			}  else if (header.defcel) 
				{
					intrj.read((char *)&tcela,sizeof(tcela));
					//for(i=0;i<6;i++) intrj.read((char *)&s2r[i],sizeof(s2r[i]));
                                        intrj.read((char *)&(cell->H[0][0]),sizeof(cell->H[0][0]));
                                        intrj.read((char *)&(cell->H[1][1]),sizeof(cell->H[1][1]));
                                        intrj.read((char *)&(cell->H[2][2]),sizeof(cell->H[2][2]));
                                        intrj.read((char *)&(cell->H[2][1]),sizeof(cell->H[2][1]));
                                        intrj.read((char *)&(cell->H[2][0]),sizeof(cell->H[2][0]));
                                        intrj.read((char *)&(cell->H[1][0]),sizeof(cell->H[1][0]));
                                        intrj.read((char *)&crtn,sizeof(crtn));
			}

			if (header.period) 
			{
				  if ( (header.version==210) || (header.version>=300) ) 
				  {
					  intrj.read((char *)&natmcel,sizeof(natmcel));
					  for(i=0;i<6;i++) intrj.read((char *)&(prp->strs[i]),sizeof(prp->strs[i]));
					  intrj.read((char *)&extstrs,sizeof(extstrs));
					  for(i=0;i<6;i++) intrj.read((char *)&strsa[i],sizeof(strsa[i]));
					  intrj.read((char *)&extstrsa,sizeof(extstrsa));
                                          intrj.read((char *)&crtn,sizeof(crtn));
				  } else
				  {
					intrj.read((char *)&natmcel,sizeof(natmcel));
					for(i=0;i<6;i++) intrj.read((char *)&(prp->strs[i]),sizeof(prp->strs[i]));
                                        intrj.read((char *)&crtn,sizeof(crtn));
				  }
			}

			if ( header.version>=300 ) 
			{
				  if (header.period && header.lnpecan ) 
				  {
					  intrj.read((char *)&sigdyn[0],sizeof(sigdyn[0]));
					  intrj.read((char *)&sigdyn[1],sizeof(sigdyn[1]));
					  intrj.read((char *)&qcanon,sizeof(qcanon));
                                          intrj.read((char *)&crtn,sizeof(crtn));
				  }
				  if (header.ltmpdamp) 
                                  {
                                          intrj.read((char *)&gamtmp,sizeof(gamtmp));
                                          intrj.read((char *)&crtn,sizeof(crtn));
                                  }
			}

			if (header.prtthrm)
			{
				intrj.read((char *)&eabtota,sizeof(eabtota));
				intrj.read((char *)&eabvala,sizeof(eabvala));
				intrj.read((char *)&eabelha,sizeof(eabelha));
				intrj.read((char *)&eabnba,sizeof(eabnba));
				intrj.read((char *)&eabmisa,sizeof(eabmisa));
				intrj.read((char *)&dltfaba,sizeof(dltfaba));
				intrj.read((char *)&expprta,sizeof(expprta));
                                intrj.read((char *)&crtn,sizeof(crtn));
			}

           return 0;
}

int TRJ::ReadBinFrame()
{
        int i;
        double crtn; //two carriage returns (each carriage return in fortran is 4 bytes)
        CleanTRJcontentDouble();
        ReadBinEng();

	for(i=0;i<header.nmovatm;i++) intrj.read((char *)&(atom[i].pv[0]),sizeof((atom[i].pv[0])));
        intrj.read((char *)&crtn,sizeof(crtn));
	for(i=0;i<header.nmovatm;i++) intrj.read((char *)&(atom[i].pv[1]),sizeof((atom[i].pv[1])));
        intrj.read((char *)&crtn,sizeof(crtn));
	for(i=0;i<header.nmovatm;i++) intrj.read((char *)&(atom[i].pv[2]),sizeof((atom[i].pv[2])));
        intrj.read((char *)&crtn,sizeof(crtn));

	//     ----- velocities if needed -----

	if ( lvelwr ) 
	{
		for(i=0;i<header.nmovatm;i++) intrj.read((char *)&(atom[i].vel[0]),sizeof(atom[i].vel[0]));
                intrj.read((char *)&crtn,sizeof(crtn));
		for(i=0;i<header.nmovatm;i++) intrj.read((char *)&(atom[i].vel[1]),sizeof((atom[i].vel[1])));
                intrj.read((char *)&crtn,sizeof(crtn));
		for(i=0;i<header.nmovatm;i++) intrj.read((char *)&(atom[i].vel[2]),sizeof((atom[i].vel[2])));
                intrj.read((char *)&crtn,sizeof(crtn));
	}
        //     ----- atomic energy if needed -----
        if ( lengwr )
        {
                for(i=0;i<header.nmovatm;i++) intrj.read((char *)&(atom[i].eng),sizeof(atom[i].eng));
                intrj.read((char *)&crtn,sizeof(crtn));
        }
        return 0;
}

int TRJ::ReadAscFrame()
{
        int i;
        double null;  //dummy variable
        CleanTRJcontentDouble();

        intrj>>prp->time>>prp->step>>prp->T>>avetem>>dstep>>firstt
                >>finalt>>prp->Ep>>prp->Ebond>>prp->Eangle>>prp->Etorsion>>prp->Einversion
                >>prp->Evdw>>prp->Eel>>prp->Ehb>>ec>>eu>>tint
                >>tnb>>ea>>eba>>eta>>epa>>eia
                >>enba>>eela>>ehba>>eca>>eua>>tinta
                >>tnba>>prp->Et>>prp->Ek>>totea>>tkea;
//cout<<"stlin "<<prp->time<<" "<<intrj.tellg()<<endl;;
        for(i=0;i<12;i++) intrj>>dum[i];
        for(i=0;i<12;i++) intrj>>duma[i];
        intrj>>iconmp>>imstep>>lvelwr>>lfrcwr;
        if (header.version>=2012) intrj>>lengwr;
        intrj>>iconfs>>icstep;

        if (header.version>=155) {
                intrj>>pressur>>vol>>pvtot
                >>pvkin>>pvpot>>radgyr;
                intrj>>prp->P>>prp->V>>pvtota
                >>pvkina>>pvpota>>radgyra;
        }

        if (header.lcanon) {
           if (header.version<300) {
                 intrj>>signose>>zfrict>>zprfrict;
           } else {
                if(header.lnose) {
                     intrj>>snose>>snoseh>>ssdot>>qcanon;
                } else {
                     intrj>>signose>>zfrict>>zprfrict>>qcanon;
                }
           }
        }

        if(header.version>=220) {
           if(header.period) {
                   intrj>>tcel>>tcela;
                   //for(i=0;i<6;i++) intrj>>s2r[i];
                   intrj>>cell->H[0][0]>>cell->H[1][1]>>cell->H[2][2]>>cell->H[2][1]>>cell->H[2][0]>>cell->H[1][0];
                   for(i=0;i<6;i++) intrj>>s2rdot[i];
                   intrj>>ucel>>ucela;
                   for(i=0;i<6;i++) intrj>>s2ra[i];
           }
/*stlin 8/27/2002*/
        }  else if (header.defcel) {
            intrj>>tcela;
            //for(i=0;i<6;i++) intrj>>s2r[i];
            intrj>>cell->H[0][0]>>cell->H[1][1]>>cell->H[2][2]>>cell->H[2][1]>>cell->H[2][0]>>cell->H[1][0];
        }

        if (header.period) {
            if ( (header.version==210) || (header.version>=300) ) {
                    intrj>>natmcel;
                    for(i=0;i<6;i++) intrj>>prp->strs[i];
                    intrj>>extstrs;
                    for(i=0;i<6;i++) intrj>>strsa[i];
                    intrj>>extstrsa;
            } else {
                  intrj>>natmcel;
                  for(i=0;i<6;i++) intrj>>prp->strs[i];
            }
        }

        if ( header.version>=300 ) {
             if (header.period && header.lnpecan ) {
                intrj>>sigdyn[0]>>sigdyn[1]>>qcanon;
             }
             if (header.ltmpdamp)  intrj>>gamtmp;
        }

        if (header.prtthrm) {
            intrj>>eabtota>>eabvala>>eabelha>>eabnba
                     >>eabmisa>>dltfaba>>expprta;
        }

        for(i=0;i<header.nmovatm;i++) intrj>>atom[i].pv[0];
        for(i=0;i<header.nmovatm;i++) intrj>>atom[i].pv[1];
        for(i=0;i<header.nmovatm;i++) intrj>>atom[i].pv[2];

//     ----- velocities if needed -----

        if ( lvelwr ) {
                for(i=0;i<header.nmovatm;i++) intrj>>atom[i].vel[0];
                for(i=0;i<header.nmovatm;i++) intrj>>atom[i].vel[1];
                for(i=0;i<header.nmovatm;i++) intrj>>atom[i].vel[2];
        }

//     ----- atomic energy if needed -----
        if ( lengwr ) {
                for(i=0;i<header.nmovatm;i++) intrj>>atom[i].eng;
        }
}

int TRJ::ReadUSRFrame()
{
        int i,j;
        double la,lb,lc;
        CleanTRJcontentDouble();
 
        intrj.read((char *)&(prp->step),sizeof(prp->step)); //number of steps
        intrj.read((char *)&(prp->time),sizeof(prp->time)); //time step (ps)
        intrj.read((char *)&(prp->T),sizeof(prp->T)); //instantaneous temperature
        intrj.read((char *)&la,sizeof(double)); //la
        intrj.read((char *)&lb,sizeof(double)); //lb
        intrj.read((char *)&lc,sizeof(double)); //lc
        prp->V=cell->volume=la*lb*lc;

	for(i=0;i<header.nmovatm;i++) {
           for(j=0;j<3;j++){
               intrj.read((char *)&(atom[i].pv[j]),sizeof(atom[i].pv[j]));
               intrj.read((char *)&(atom[i].vel[j]),sizeof(atom[i].vel[j]));
               atom[i].vel[j]*=sqrt(418.4032);
           }
        }

//for(i=0;i<header.nmovatm;i++) printf("atom %d v %8.6f %8.6f %8.6f c %8.6f %8.6f %8.6f\n",i,atom[i].vel[0],atom[i].vel[1],atom[i].vel[2],atom[i].pv[0],atom[i].pv[1],atom[i].pv[2]);

}

int TRJ::ReadCPMDFrame()
{
    int i,j,skip;
    double tmp;  //dummy variable
    char null[1024];
    CleanTRJcontentDouble();
    //read xyz and vel
    dstep=header.timestep;
    for(skip=0;skip<header.nxyz;skip++) {
        for(i=0;i<header.nmovatm;i++) {
          intrj>>null;
          if(strcmp(null,"<<<<<<")==0) {
             intrj.getline(null,1024); intrj>>null;
          }
          prp->step=(int)atof(null);
          intrj>>atom[i].pv[0]>>atom[i].pv[1]>>atom[i].pv[2]>>atom[i].vel[0]>>atom[i].vel[1]>>atom[i].vel[2];
          prp->step+=header.totaccustep;
          for(j=0;j<3;j++) {
             atom[i].pv[j]*=Bohr2A; //in A
             atom[i].vel[j]*=2188491.52E-2; //in A/ps
          }
          prp->time=(prp->step-1.0)*(header.timestep);
          //if(prp->step>99999&&prp->step<100004) printf("%d %s %f %f %f %f\n",prp->step,null,prp->time,atom[i].pv[0],atom[i].pv[1],atom[i].pv[2]);
        }
    }
    //read energy
    ineng>>i>>prp->Eke>>prp->T>>prp->Ep>>prp->Et>>prp->Ete>>tmp>>tmp; 
    if(i!=prp->step-header.totaccustep) cout<<"WARNING: stept "<<prp->step<<" stepe "<<i<<endl;
    (prp->Eke) *= Hartree2kcalmol;
    (prp->Ep)  *= Hartree2kcalmol;
    (prp->Et)  *= Hartree2kcalmol;
    (prp->Ete) *= Hartree2kcalmol;
    prp->Ek=(prp->Et)-(prp->Ep);
    prp->V=cell->volume;
    for(skip=0;skip<header.neng-1;skip++) ineng>>null>>null>>null>>null>>null>>null>>null>>null;
     
    //read stress
    //instr.getline(null,1024);
    instr>>null;
    if(strcmp(null,"<<<<<<")==0) instr>>null>>null>>null>>null;
    instr>>null>>null>>null>>null>>i;
    if(i!=prp->step-header.totaccustep) cout<<"WARNING: stept "<<prp->step<<" steps "<<i<<endl;
    instr>>prp->strs[0]>>prp->strs[3]>>prp->strs[4]>>tmp>>prp->strs[1]>>prp->strs[5]>>tmp>>tmp>>prp->strs[2]; 
    for(i=0;i<6;i++) prp->strs[i]*=-0.1; //in GPa
    prp->P=(prp->strs[0]+prp->strs[1]+prp->strs[2])/3;
    for(skip=0;skip<header.nstr-1;skip++) {
       instr>>null;
       if(strcmp(null,"<<<<<<")==0) instr.getline(null,1024);
       instr.getline(null,1024); 
       instr.getline(null,1024); 
       instr.getline(null,1024); 
       instr.getline(null,1024); 
    }
    return 0;
}

int TRJ::ReadLMPFrame()
{
    int i,j,skip;
    double tmp;  //dummy variable
    char null[1024];
    CleanTRJcontentDouble();
    //read xyz and vel
    dstep=header.timestep;
    for(skip=0;skip<header.nxyz;skip++) {
        for(i=0;i<header.nmovatm;i++) {
          intrj>>null;
          if(strcmp(null,"<<<<<<")==0) {
             intrj.getline(null,1024); intrj>>null;
          }
          prp->step=(int)atof(null);
          intrj>>atom[i].pv[0]>>atom[i].pv[1]>>atom[i].pv[2]>>atom[i].vel[0]>>atom[i].vel[1]>>atom[i].vel[2];
          prp->step+=header.totaccustep;
          for(j=0;j<3;j++) {
             atom[i].pv[j]*=Bohr2A; //in A
             atom[i].vel[j]*=2188491.52E-2; //in A/ps
          }
          prp->time=(prp->step-1.0)*(header.timestep);
          //if(prp->step>99999&&prp->step<100004) printf("%d %s %f %f %f %f\n",prp->step,null,prp->time,atom[i].pv[0],atom[i].pv[1],atom[i].pv[2]);
        }
    }
    //read energy
    ineng>>i>>prp->Eke>>prp->T>>prp->Ep>>prp->Et>>prp->Ete>>tmp>>tmp; 
    if(i!=prp->step-header.totaccustep) cout<<"WARNING: stept "<<prp->step<<" stepe "<<i<<endl;
    (prp->Eke) *= Hartree2kcalmol;
    (prp->Ep)  *= Hartree2kcalmol;
    (prp->Et)  *= Hartree2kcalmol;
    (prp->Ete) *= Hartree2kcalmol;
    prp->Ek=(prp->Et)-(prp->Ep);
    prp->V=cell->volume;
    for(skip=0;skip<header.neng-1;skip++) ineng>>null>>null>>null>>null>>null>>null>>null>>null;
     
    //read stress
    //instr.getline(null,1024);
    instr>>null;
    if(strcmp(null,"<<<<<<")==0) instr>>null>>null>>null>>null;
    instr>>null>>null>>null>>null>>i;
    if(i!=prp->step-header.totaccustep) cout<<"WARNING: stept "<<prp->step<<" steps "<<i<<endl;
    instr>>prp->strs[0]>>prp->strs[3]>>prp->strs[4]>>tmp>>prp->strs[1]>>prp->strs[5]>>tmp>>tmp>>prp->strs[2]; 
    for(i=0;i<6;i++) prp->strs[i]*=-0.1; //in GPa
    prp->P=(prp->strs[0]+prp->strs[1]+prp->strs[2])/3;
    for(skip=0;skip<header.nstr-1;skip++) {
       instr>>null;
       if(strcmp(null,"<<<<<<")==0) instr.getline(null,1024);
       instr.getline(null,1024); 
       instr.getline(null,1024); 
       instr.getline(null,1024); 
       instr.getline(null,1024); 
    }
    return 0;
}

int TRJ::WriteASCIIcontent(ostream *outf, int oformat, int ndel, int iddel[])
{
	int i,j;
	switch(oformat)
	{
	case 0:  //convert all
        	*outf<<setprecision(20)<<endl;
		*outf<<prp->time<<" "<<prp->step<<" "<<prp->T
			<<" "<<avetem<<" "<<dstep<<" "<<firstt
			<<" "<<finalt<<" "<<prp->Ep<<" "<<prp->Ebond<<" "<<prp->Eangle
			<<" "<<prp->Etorsion<<" "<<prp->Einversion<<" "<<prp->Evdw
                        <<" "<<prp->Eel<<" "<<prp->Ehb<<" "<<ec<<" "<<eu<<" "<<tint
			<<" "<<tnb<<" "<<ea<<" "<<eba<<" "<<eta
			<<" "<<epa<<" "<<eia<<" "<<enba<<" "<<eela
			<<" "<<ehba<<" "<<eca<<" "<<eua<<" "<<tinta
			<<" "<<tnba<<" "<<prp->Et<<" "<<prp->Ek<<" "<<totea
			<<" "<<tkea<<" ";
		for(i=0;i<12;i++) *outf<<dum[i]<<" ";
		for(i=0;i<12;i++) *outf<<duma[i]<<" ";
                *outf<<iconmp<<" "<<imstep<<" "<<lvelwr<<" "<<lfrcwr;
                if (header.version>=2012) *outf<<" "<<lengwr;
		*outf<<" "<<iconfs<<" "<<icstep<<endl;

		if (header.version>=155) 
		{
			*outf<<pressur<<" "<<vol<<" "<<pvtot<<" "
				<<pvkin<<" "<<pvpot<<" "<<radgyr<<" ";
			*outf<<prp->P<<" "<<prp->V<<" "<<pvtota<<" "
				<<pvkina<<" "<<pvpota<<" "<<radgyra<<" ";
		}

		if (header.lcanon) 
		{
			  if (header.version<300) 
			  {
				*outf<<signose<<" "<<zfrict<<" "<<zprfrict<<endl;
           	  } else
			  {
				  if(header.lnose) 
				  {
					  *outf<<snose<<" "<<snoseh<<" "<<ssdot
						  <<" "<<qcanon<<endl;
				  }	else 
				  {
					  *outf<<signose<<" "<<zfrict<<" "<<zprfrict
						  <<" "<<qcanon<<endl;
				  }
			  }
		}
		
		if(header.version>=220)
		{
			if(header.period) 
			{
				*outf<<tcel<<" "<<tcela<<" ";
				//for(i=0;i<6;i++) *outf<<s2r[i]<<" ";
                                *outf<<cell->H[0][0]<<" "<<cell->H[1][1]<<" "<<cell->H[2][2]<<" "<<cell->H[2][1]<<" "<<cell->H[2][0]<<" "<<cell->H[1][0]<<" ";
				for(i=0;i<6;i++) *outf<<s2rdot[i]<<" ";
				*outf<<ucel<<" "<<ucela<<" ";
				for(i=0;i<6;i++) *outf<<s2ra[i]<<" ";
			}
/* stlin 8/27/2002*/
		}  else if (header.defcel) 
			{
				*outf<<endl<<tcela<<" ";
				//for(i=0;i<6;i++) *outf<<s2r[i]<<" ";
                                *outf<<cell->H[0][0]<<" "<<cell->H[1][1]<<" "<<cell->H[2][2]<<" "<<cell->H[2][1]<<" "<<cell->H[2][0]<<" "<<cell->H[1][0]<<" ";
		}

		if (header.period) 
		{
			  if ( (header.version==210) || (header.version>=300) ) 
			  {
				  *outf<<natmcel<<" ";
				  for(i=0;i<6;i++) *outf<<prp->strs[i]<<" ";
				  *outf<<extstrs<<endl;
				  for(i=0;i<6;i++) *outf<<strsa[i]<<" ";
				  *outf<<extstrsa<<endl;
			  } else
			  {
				*outf<<natmcel<<" ";
				for(i=0;i<6;i++) *outf<<prp->strs[i]<<" ";
			  }
		}

		if ( header.version>=300 ) 
		{
			  if (header.period && header.lnpecan ) 
			  {
				  *outf<<sigdyn[0]<<" "<<sigdyn[1]<<" "<<qcanon<<endl;
			  }
			  if (header.ltmpdamp) 
                              {
                                      *outf<<gamtmp<<endl;
                              }
		}

		if (header.prtthrm)
		{
			*outf<<eabtota<<" "<<eabvala<<" "<<eabelha<<" "<<eabnba
				<<" "<<eabmisa<<" "<<dltfaba<<" "<<expprta<<endl;
		}

		for(i=0;i<header.nmovatm;i++) *outf<<atom[i].pv[0]<<" ";
		*outf<<endl;
		for(i=0;i<header.nmovatm;i++) *outf<<atom[i].pv[1]<<" ";
		*outf<<endl;
		for(i=0;i<header.nmovatm;i++) *outf<<atom[i].pv[2]<<" ";
		*outf<<endl;

	//     ----- velocities if needed -----

		if ( lvelwr ) 
		{
			for(i=0;i<header.nmovatm;i++) *outf<<atom[i].vel[0]<<" ";
			*outf<<endl;
			for(i=0;i<header.nmovatm;i++) *outf<<atom[i].vel[1]<<" ";
			*outf<<endl;
			for(i=0;i<header.nmovatm;i++) *outf<<atom[i].vel[2]<<" ";
			*outf<<endl;
		}
        //     ----- atomic energy if needed -----		
                if ( lengwr )
                {
                        for(i=0;i<header.nmovatm;i++) *outf<<atom[i].eng<<" ";
                        *outf<<endl;
                }
		break;
	case 1: // just the energies
		*outf<<setiosflags(ios::fixed|ios::right)
			<<setw(9)<<setprecision(3)<<prp->time  //time in ps
			<<setw(9)<<setprecision(3)<<prp->T  //temperature
			<<setw(14)<<setprecision(4)<<prp->Et   //total energy
		      //<<setw(15)<<setprecision(4)<<prp->Ek  //kinetic energy
			<<setw(15)<<setprecision(4)<<prp->Ep  //potential energy
			<<setw(15)<<setprecision(4)<<prp->Eel    //electrostatic
			<<setw(15)<<setprecision(4)<<prp->Evdw    //van der waals
			<<setw(15)<<setprecision(4)<<tint   //valence energy
			<<setw(15)<<setprecision(4)<<prp->Ehb    //h-bond energy
			<<endl;
		break;
	case 2: /*temporary format for generating trace of certain atoms */
		j=0;
		for (i=1008;i<header.nmovatm;i+=5) {
			*outf<<setw(5)<<i+1<<setprecision(4)<<setiosflags(ios::fixed)
			     <<setw(10)<<atom[i].pv[0]<<setw(10)<<atom[i].pv[1]<<setw(10)<<atom[i].pv[2]
			     <<setw(10)<<1<<setw(3)<<0<<setw(10)<<(i+1-j-1000.0-30)/1000.0<<setw(10)<<0
			     <<setw(5)<<1<<setw(5)<<1<<setw(5)<<1<<setw(10)<<0<<setw(10)<<0<<setw(10)<<0
			     <<endl;
			//j++;
			//if(j==5) j=0;
		}
		break;
	case 3: /* removing certain number of atoms */
        	*outf<<setprecision(20)<<endl;
		*outf<<prp->time<<" "<<prp->step<<" "<<prp->T
			<<" "<<avetem<<" "<<dstep<<" "<<firstt
			<<" "<<finalt<<" "<<prp->Ep<<" "<<prp->Ebond<<" "<<prp->Eangle
			<<" "<<prp->Etorsion<<" "<<prp->Einversion<<" "<<prp->Evdw
                        <<" "<<prp->Eel<<" "<<prp->Ehb<<" "<<ec<<" "<<eu<<" "<<tint
			<<" "<<tnb<<" "<<ea<<" "<<eba<<" "<<eta
			<<" "<<epa<<" "<<eia<<" "<<enba<<" "<<eela
			<<" "<<ehba<<" "<<eca<<" "<<eua<<" "<<tinta
			<<" "<<tnba<<" "<<prp->Et<<" "<<prp->Ek<<" "<<totea
			<<" "<<tkea<<" ";
		for(i=0;i<12;i++) *outf<<dum[i]<<" ";
		for(i=0;i<12;i++) *outf<<duma[i]<<" ";
                *outf<<iconmp<<" "<<imstep<<" "<<lvelwr<<" "<<lfrcwr;
                if (header.version>=2012) *outf<<" "<<lengwr;
		*outf<<" "<<iconfs<<" "<<icstep<<endl;

		if (header.version>=155) 
		{
			*outf<<pressur<<" "<<vol<<" "<<pvtot<<" "
				<<pvkin<<" "<<pvpot<<" "<<radgyr<<" ";
			*outf<<prp->P<<" "<<prp->V<<" "<<pvtota<<" "
				<<pvkina<<" "<<pvpota<<" "<<radgyra<<" ";
		}

		if (header.lcanon) 
		{
			  if (header.version<300) 
			  {
				*outf<<signose<<" "<<zfrict<<" "<<zprfrict<<endl;
           	  } else
			  {
				  if(header.lnose) 
				  {
					  *outf<<snose<<" "<<snoseh<<" "<<ssdot
						  <<" "<<qcanon<<endl;
				  }	else 
				  {
					  *outf<<signose<<" "<<zfrict<<" "<<zprfrict
						  <<" "<<qcanon<<endl;
				  }
			  }
		}
		
		if(header.version>=220)
		{
			if(header.period) 
			{
				*outf<<tcel<<" "<<tcela<<" ";
				//for(i=0;i<6;i++) *outf<<s2r[i]<<" ";
                                *outf<<cell->H[0][0]<<" "<<cell->H[1][1]<<" "<<cell->H[2][2]<<" "<<cell->H[2][1]<<" "<<cell->H[2][0]<<" "<<cell->H[1][0]<<" ";
				for(i=0;i<6;i++) *outf<<s2rdot[i]<<" ";
				*outf<<ucel<<" "<<ucela<<" ";
				for(i=0;i<6;i++) *outf<<s2ra[i]<<" ";
			}
/*stlin 8/27/2002*/
		}  else if (header.defcel) 
			{
				*outf<<endl<<tcela<<" ";
				//for(i=0;i<6;i++) *outf<<s2r[i]<<" ";
                                *outf<<cell->H[0][0]<<" "<<cell->H[1][1]<<" "<<cell->H[2][2]<<" "<<cell->H[2][1]<<" "<<cell->H[2][0]<<" "<<cell->H[1][0]<<" ";
		}

		if (header.period) 
		{
			  if ( (header.version==210) || (header.version>=300) ) 
			  {
				  *outf<<natmcel-ndel<<" ";
				  for(i=0;i<6;i++) *outf<<prp->strs[i]<<" ";
				  *outf<<extstrs<<endl;
				  for(i=0;i<6;i++) *outf<<strsa[i]<<" ";
				  *outf<<extstrsa<<endl;
			  } else
			  {
				*outf<<natmcel<<" ";
				for(i=0;i<6;i++) *outf<<prp->strs[i]<<" ";
			  }
		}

		if ( header.version>=300 ) 
		{
			  if (header.period && header.lnpecan ) 
			  {
				  *outf<<sigdyn[0]<<" "<<sigdyn[1]<<" "<<qcanon<<endl;
			  }
			  if (header.ltmpdamp) 
                              {
                                      *outf<<gamtmp<<endl;
                              }
		}

		if (header.prtthrm)
		{
			*outf<<eabtota<<" "<<eabvala<<" "<<eabelha<<" "<<eabnba
				<<" "<<eabmisa<<" "<<dltfaba<<" "<<expprta<<endl;
		}

		for(i=0;i<header.nmovatm;i++) {
			if(iddel[i]) *outf<<atom[i].pv[0]<<" ";
		}
		*outf<<endl;
		for(i=0;i<header.nmovatm;i++) {
			if(iddel[i]) *outf<<atom[i].pv[1]<<" ";
		}
		*outf<<endl;
		for(i=0;i<header.nmovatm;i++) { 
			if(iddel[i]) *outf<<atom[i].pv[2]<<" ";
		}
		*outf<<endl;

	//     ----- velocities if needed -----

		if ( lvelwr ) 
		{
			for(i=0;i<header.nmovatm;i++) {
				if(iddel[i]) {
					 *outf<<atom[i].vel[0]<<" ";
				}
			}
			*outf<<endl;
			for(i=0;i<header.nmovatm;i++) {
				if(iddel[i]) {
                                         *outf<<atom[i].vel[1]<<" ";
				}
			}
			*outf<<endl;
			for(i=0;i<header.nmovatm;i++) {
				if(iddel[i]) {
                                         *outf<<atom[i].vel[2]<<" ";
				}
			}
			*outf<<endl;
		}		
        //     ----- atomic energy if needed -----
                if ( lengwr )
                {
                        for(i=0;i<header.nmovatm;i++) {
                                if(iddel[i]) {
                                         *outf<<atom[i].eng<<" ";
                                }
                        }
                        *outf<<endl;
                }
		break;
	case 4:
		break;
	case 5:
		break;
	default:
		break;
	}
        return 0;
}

int TRJ::WriteBinarycontent(ostream *outf, int oformat)
{
	int i;
	int size;

	size=sizeof(prp->time)+sizeof(prp->step)+sizeof(prp->T)+sizeof(avetem)
		+sizeof(dstep)+sizeof(firstt)+sizeof(finalt)+sizeof(prp->Ep)
		+sizeof(prp->Ebond)+sizeof(prp->Eangle)+sizeof(prp->Etorsion)
                +sizeof(prp->Einversion)+sizeof(prp->Evdw)
		+sizeof(prp->Eel)+sizeof(prp->Ehb)+sizeof(ec)+sizeof(eu)+sizeof(tint)
		+sizeof(tnb)+sizeof(ea)+sizeof(eba)+sizeof(eta)+sizeof(epa)
		+sizeof(eia)+sizeof(enba)+sizeof(eela)+sizeof(ehba)+sizeof(eca)
		+sizeof(eua)+sizeof(tinta)+sizeof(tnba)+sizeof(prp->Et)
		+sizeof(prp->Ek)+sizeof(totea)+sizeof(tkea)+12*sizeof(dum[0])
                +12*sizeof(duma[0])+sizeof(iconmp)+sizeof(imstep)
		+sizeof(lvelwr)+sizeof(lfrcwr)+sizeof(iconfs)+sizeof(icstep);
        if(header.version>=2012) size+= sizeof(lengwr);

 	outf->write((char *)&size,sizeof(size));

	outf->write((char *)&(prp->time),sizeof(prp->time));
	outf->write((char *)&(prp->step),sizeof(prp->step));
	outf->write((char *)&(prp->T),sizeof(prp->T));
	outf->write((char *)&avetem,sizeof(avetem));
	outf->write((char *)&dstep,sizeof(dstep));
	outf->write((char *)&firstt,sizeof(firstt));
	outf->write((char *)&finalt,sizeof(finalt));
	outf->write((char *)&(prp->Ep),sizeof(prp->Ep));
	outf->write((char *)&(prp->Ebond),sizeof(prp->Ebond));
	outf->write((char *)&(prp->Eangle),sizeof(prp->Eangle));
	outf->write((char *)&(prp->Etorsion),sizeof(prp->Etorsion));
	outf->write((char *)&(prp->Einversion),sizeof(prp->Einversion));
	outf->write((char *)&(prp->Evdw),sizeof(prp->Evdw));
	outf->write((char *)&(prp->Eel),sizeof(prp->Eel));
	outf->write((char *)&(prp->Ehb),sizeof(prp->Ehb));
	outf->write((char *)&ec,sizeof(ec));
        for(i=0;i<12;i++) outf->write((char *)&dum[i],sizeof(dum[i]));
	outf->write((char *)&eu,sizeof(eu));
	outf->write((char *)&tint,sizeof(tint));
	outf->write((char *)&tnb,sizeof(tnb));
	outf->write((char *)&ea,sizeof(ea));
	outf->write((char *)&eba,sizeof(eba));
	outf->write((char *)&eta,sizeof(eta));
	outf->write((char *)&epa,sizeof(epa));
	outf->write((char *)&eia,sizeof(eia));
	outf->write((char *)&enba,sizeof(enba));
	outf->write((char *)&eela,sizeof(eela));
	outf->write((char *)&ehba,sizeof(ehba));
	outf->write((char *)&eca,sizeof(eca));
	outf->write((char *)&eua,sizeof(eua));
	outf->write((char *)&tinta,sizeof(tinta));
        for(i=0;i<12;i++) outf->write((char *)&duma[i],sizeof(duma[i]));
	outf->write((char *)&tnba,sizeof(tnba));
	outf->write((char *)&(prp->Et),sizeof(prp->Et));
	outf->write((char *)&(prp->Ek),sizeof(prp->Ek));
	outf->write((char *)&totea,sizeof(totea));
	outf->write((char *)&tkea,sizeof(tkea));
	outf->write((char *)&iconmp,sizeof(iconmp));
	outf->write((char *)&imstep,sizeof(imstep));
	outf->write((char *)&lvelwr,sizeof(lvelwr));
	outf->write((char *)&lfrcwr,sizeof(lfrcwr));
        if(header.version>=2012) outf->write((char *)&lengwr,sizeof(lengwr));
	outf->write((char *)&iconfs,sizeof(iconfs));
	outf->write((char *)&icstep,sizeof(icstep));
 	outf->write((char *)&size,sizeof(size));

	if (header.version>=155) 
	{
		size=sizeof(prp->P)+sizeof(prp->V)+sizeof(pvtota)+sizeof(pvkina)+sizeof(pvpota)+sizeof(radgyra)+sizeof(pressur)+sizeof(vol)+sizeof(pvtot)+sizeof(pvkin)+sizeof(pvpot)+sizeof(radgyr);
		outf->write((char *)&size,sizeof(size));
		outf->write((char *)&pressur,sizeof(pressur));
		outf->write((char *)&vol,sizeof(vol));
		outf->write((char *)&pvtot,sizeof(pvtot));
		outf->write((char *)&pvkin,sizeof(pvkin));
		outf->write((char *)&pvpot,sizeof(pvpot));
		outf->write((char *)&radgyr,sizeof(radgyr));
		outf->write((char *)&(prp->P),sizeof(prp->P));
		outf->write((char *)&(prp->V),sizeof(prp->V));
		outf->write((char *)&pvtota,sizeof(pvtota));
		outf->write((char *)&pvkina,sizeof(pvkina));
		outf->write((char *)&pvpota,sizeof(pvpota));
		outf->write((char *)&radgyra,sizeof(radgyra));
	 	outf->write((char *)&size,sizeof(size));
	}

	if (header.lcanon) 
	{
		  if (header.version<300) 
		  {
			size=sizeof(signose)+sizeof(zfrict)+sizeof(zprfrict);
		 	outf->write((char *)&size,sizeof(size));
			outf->write((char *)&signose,sizeof(signose));
			outf->write((char *)&zfrict,sizeof(zfrict));
			outf->write((char *)&zprfrict,sizeof(zprfrict));
		 	outf->write((char *)&size,sizeof(size));
		  } else
		  {
			  if(header.lnose) 
			  {
				  size=sizeof(snose)+sizeof(snoseh)+sizeof(ssdot)+sizeof(qcanon);
			      outf->write((char *)&size,sizeof(size));
				  outf->write((char *)&snose,sizeof(snose));
				  outf->write((char *)&snoseh,sizeof(snoseh));
				  outf->write((char *)&ssdot,sizeof(ssdot));
				  outf->write((char *)&qcanon,sizeof(qcanon));
                  outf->write((char *)&size,sizeof(size));
			  }	else 
			  {
				  size=sizeof(signose)+sizeof(zfrict)+sizeof(zprfrict)+sizeof(qcanon);
			      outf->write((char *)&size,sizeof(size));
				  outf->write((char *)&signose,sizeof(signose));
				  outf->write((char *)&zfrict,sizeof(zfrict));
				  outf->write((char *)&zprfrict,sizeof(zprfrict));
				  outf->write((char *)&qcanon,sizeof(qcanon));
                  outf->write((char *)&size,sizeof(size));
			  }
		  }
	}
	
	if(header.version>=220)
	{
		if(header.period) 
		{
			size=sizeof(tcel)+sizeof(tcela)+6*sizeof(cell->H[0][0])+6*sizeof(s2rdot[i])+sizeof(ucel)+sizeof(ucela)+6*sizeof(s2ra[0]);
			outf->write((char *)&size,sizeof(size));
			outf->write((char *)&tcel,sizeof(tcel));
			outf->write((char *)&tcela,sizeof(tcela));
			//for(i=0;i<6;i++) outf->write((char *)&s2r[i],sizeof(s2r[i]));
                        outf->write((char *)&(cell->H[0][0]),sizeof(cell->H[0][0]));
                        outf->write((char *)&(cell->H[1][1]),sizeof(cell->H[1][1]));
                        outf->write((char *)&(cell->H[2][2]),sizeof(cell->H[2][2]));
                        outf->write((char *)&(cell->H[2][1]),sizeof(cell->H[2][1]));
                        outf->write((char *)&(cell->H[2][0]),sizeof(cell->H[2][0]));
                        outf->write((char *)&(cell->H[1][0]),sizeof(cell->H[1][0]));
			for(i=0;i<6;i++) outf->write((char *)&s2rdot[i],sizeof(s2rdot[i]));
			outf->write((char *)&ucel,sizeof(ucel));
			outf->write((char *)&ucela,sizeof(ucela));
			for(i=0;i<6;i++) outf->write((char *)&s2ra[i],sizeof(s2ra[i]));
            		outf->write((char *)&size,sizeof(size));
		} 
/*stlin 8/27/2002 */
	} else if (header.defcel) 
		{
			size=sizeof(tcela)+6*sizeof(cell->H[0][0]);
			outf->write((char *)&size,sizeof(size));
			outf->write((char *)&tcela,sizeof(tcela));
			//for(i=0;i<6;i++) outf->write((char *)&s2r[i],sizeof(s2r[i]));
                        outf->write((char *)&(cell->H[0][0]),sizeof(cell->H[0][0]));
                        outf->write((char *)&(cell->H[1][1]),sizeof(cell->H[1][1]));
                        outf->write((char *)&(cell->H[2][2]),sizeof(cell->H[2][2]));
                        outf->write((char *)&(cell->H[2][1]),sizeof(cell->H[2][1]));
                        outf->write((char *)&(cell->H[2][0]),sizeof(cell->H[2][0]));
                        outf->write((char *)&(cell->H[1][0]),sizeof(cell->H[1][0]));
			outf->write((char *)&size,sizeof(size));
	}

	if (header.period) 
	{
		  if ( (header.version==210) || (header.version>=300) ) 
		  {
			  size=sizeof(natmcel)+6*sizeof(strsa[0])+sizeof(extstrsa)
                              +6*sizeof(prp->strs[0])+sizeof(extstrs);
			  outf->write((char *)&size,sizeof(size));
			  outf->write((char *)&natmcel,sizeof(natmcel));
			  for(i=0;i<6;i++) outf->write((char *)&(prp->strs[i]),sizeof(prp->strs[i]));
			  outf->write((char *)&extstrs,sizeof(extstrs));
			  for(i=0;i<6;i++) outf->write((char *)&strsa[i],sizeof(strsa[i]));
			  outf->write((char *)&extstrsa,sizeof(extstrsa));
			  outf->write((char *)&size,sizeof(size));
		  } else
		  {
			  size=sizeof(natmcel)+6*sizeof(prp->strs[0]);
			  outf->write((char *)&size,sizeof(size));
			  outf->write((char *)&natmcel,sizeof(natmcel));
			  for(i=0;i<6;i++) outf->write((char *)&(prp->strs[i]),sizeof(prp->strs[i]));
			  outf->write((char *)&size,sizeof(size));
		  }
	}

	if ( header.version>=300 ) 
	{
		  if (header.period && header.lnpecan ) 
		  {
			  size=2*sizeof(sigdyn[0])+sizeof(qcanon);
			  outf->write((char *)&size,sizeof(size));
			  outf->write((char *)&sigdyn[0],sizeof(sigdyn[0]));
			  outf->write((char *)&sigdyn[1],sizeof(sigdyn[1]));
			  outf->write((char *)&qcanon,sizeof(qcanon));
			  outf->write((char *)&size,sizeof(size));
		  }
		  if (header.ltmpdamp) 
                          {
								  size=sizeof(gamtmp);
								  outf->write((char *)&size,sizeof(size));
                                  outf->write((char *)&gamtmp,sizeof(gamtmp));
								  outf->write((char *)&size,sizeof(size));
                          }
	}

	if (header.prtthrm)
	{
		size=sizeof(eabtota)+sizeof(eabvala)+sizeof(eabelha)+sizeof(eabnba)
                    +sizeof(eabmisa)+sizeof(dltfaba)+sizeof(expprta);
		outf->write((char *)&size,sizeof(size));
		outf->write((char *)&eabtota,sizeof(eabtota));
		outf->write((char *)&eabvala,sizeof(eabvala));
		outf->write((char *)&eabelha,sizeof(eabelha));
		outf->write((char *)&eabnba,sizeof(eabnba));
		outf->write((char *)&eabmisa,sizeof(eabmisa));
		outf->write((char *)&dltfaba,sizeof(dltfaba));
		outf->write((char *)&expprta,sizeof(expprta));
		outf->write((char *)&size,sizeof(size));
	}
        size=header.nmovatm*sizeof(atom[0].pv[0]);
        outf->write((char *)&size,sizeof(size));
	for(i=0;i<header.nmovatm;i++) outf->write((char *)&(atom[i].pv[0]),sizeof(atom[i].pv[0]));
        outf->write((char *)&size,sizeof(size));
 
        size=header.nmovatm*sizeof(atom[0].pv[1]);
        outf->write((char *)&size,sizeof(size));
	for(i=0;i<header.nmovatm;i++) outf->write((char *)&(atom[i].pv[1]),sizeof(atom[i].pv[1]));
        outf->write((char *)&size,sizeof(size));

        size=header.nmovatm*sizeof(atom[0].pv[2]);
        outf->write((char *)&size,sizeof(size));
	for(i=0;i<header.nmovatm;i++) outf->write((char *)&(atom[i].pv[2]),sizeof((atom[i].pv[2])));
        outf->write((char *)&size,sizeof(size));

//     ----- velocities if needed -----

	if ( lvelwr ) 
	{
                size=header.nmovatm*sizeof(atom[0].vel[0]);
                outf->write((char *)&size,sizeof(size));
		for(i=0;i<header.nmovatm;i++) outf->write((char *)&(atom[i].vel[0]),sizeof(atom[i].vel[0]));
                outf->write((char *)&size,sizeof(size));
 
                size=header.nmovatm*sizeof(atom[0].vel[1]);
                outf->write((char *)&size,sizeof(size));
		for(i=0;i<header.nmovatm;i++) outf->write((char *)&(atom[i].vel[1]),sizeof((atom[i].vel[1])));
                outf->write((char *)&size,sizeof(size));

                size=header.nmovatm*sizeof(atom[0].vel[2]);
                outf->write((char *)&size,sizeof(size));
		for(i=0;i<header.nmovatm;i++) outf->write((char *)&(atom[i].vel[2]),sizeof((atom[i].vel[2])));
		outf->write((char *)&size,sizeof(size));
	}

//     ----- atomic energy if needed -----
        if ( lengwr )
        {
                size=header.nmovatm*sizeof(atom[0].eng);
                outf->write((char *)&size,sizeof(size));
                for(i=0;i<header.nmovatm;i++) outf->write((char *)&(atom[i].eng),sizeof(atom[i].eng));
                outf->write((char *)&size,sizeof(size));
        }

    return 0;
}

int TRJ::init_trj(char *trjfname)
{
    int i;
    switch(itrjformat) {
      case c2trj:
           i=init_c2bintrj(trjfname);
           break;
      case asctrj:
           i=init_c2asctrj(trjfname);
           break;
      case cpmdtrj:
           i=init_cpmdtrj(trjfname);
           break;
      case usrtrj:
           i=init_usrtrj(trjfname);
           break;
      case lmptrj:
           i=init_lmptrj(trjfname);
           break;
      case none:
      default:
           i=0;
           break;
    }
    return i;
}

int TRJ::init_c2bintrj(char *trjfname)
{
    char null[1024];
    cout<<"Analyzing trajectory file "<<trjfname<<endl;
    intrj.open(trjfname,ios::in);
    if(!intrj.is_open()) {
       cout<<" Error: Trj file "<<trjfname<<" cannot be opened."<<endl;
       return filenotopen;
    }

    rd_header();  //read header

    location0=intrj.tellg();
    init_contentDouble();
    CleanTRJcontentDouble();
    ReadBinFrame(); 
    location1=intrj.tellg();
    datasize=location1-location0;

    //find the total frame number
    intrj.seekg(0, ios::end);
    location2=intrj.tellg();
    totframe=(location2-location0+4)/datasize;
    //printf("stlin l0 %d l1 %d size %d l2 %d total frame %d \n",location0,location1,datasize,location2,totframe);
    return 0;
}

int TRJ::init_c2asctrj(char *trjfname)
{
    char null[1024];
    cout<<"Analyzing trajectory file "<<trjfname<<endl;
    intrj.open(trjfname,ios::in);
    if(!intrj.is_open()) {
       cout<<" Error: Trj file "<<trjfname<<" cannot be opened."<<endl;
       return filenotopen;
    }

    rd_header();  //read header

    location0=intrj.tellg();
    init_contentDouble();
    CleanTRJcontentDouble();
    ReadAscFrame();
    
    location1=intrj.tellg();
    datasize=location1-location0;
    ascfid=1; datasize2fid=datasize;

    //find the total frame number
    totframe=1;
    while(!intrj.eof()) {
      ReadAscFrame();
      totframe++;
    }
    totframe--;
    //printf("stlin l0 %d l1 %d size %d l2 %d total frame %d \n",location0,location1,datasize,location2,totframe);
    return 0;
}

int TRJ::init_usrtrj(char *trjfname)
{

    char null[1024];
    cout<<"Analyzing user trajectory file "<<trjfname<<endl;
    intrj.open(trjfname,ios::in);
    if(!intrj.is_open()) {
       cout<<" Error: Trj file "<<trjfname<<" cannot be opened."<<endl;
       return filenotopen;
    }

    rd_header();  //read header

    location0=intrj.tellg();
    init_contentDouble();
    CleanTRJcontentDouble();
    ReadUSRFrame();
    location1=intrj.tellg();
    datasize=location1-location0;

    //find the total frame number
    intrj.seekg(0, ios::end);
    location2=intrj.tellg();
    totframe=(location2-location0+4)/datasize;
    printf("stlin usrtrj l0 %d l1 %d size %d l2 %d total frame %d \n",location0,location1,datasize,location2,totframe);

    return 0;
}

int TRJ::init_cpmdtrj(char *trjfname)
{
    char null[1024];
    cout<<"Analyzing trajectory files : "<<endl;
    strcpy(null,trjfname); strcat(null,"/TRAJECTORY");
    cout<<"TRAJECTORY : "<<null<<endl;
    intrj.open(null,ios::in);
    if(!intrj.is_open()) {
       cout<<" Error: Trj file "<<null<<" cannot be opened."<<endl;
       return filenotopen;
    }
    strcpy(null,trjfname); strcat(null,"/ENERGIES");
    cout<<"ENERGIES   : "<<null<<endl;
    ineng.open(null,ios::in);
    if(!ineng.is_open()) {
       cout<<" Error: ENERGY file "<<null<<" cannot be opened."<<endl;
       return filenotopen;
    }
    strcpy(null,trjfname); strcat(null,"/STRESS");
    cout<<"STRESS     : "<<null<<endl;
    instr.open(null,ios::in);
    if(!instr.is_open()) {
       cout<<" Error: STRESS file "<<null<<" cannot be opened."<<endl;
       return filenotopen;
    }

    rd_header();  //read header

    location0=intrj.tellg();
    location0e=ineng.tellg();
    location0s=instr.tellg();
    init_contentDouble();
    CleanTRJcontentDouble();
    ReadCPMDFrame();
    location1=intrj.tellg();
    datasize=location1-location0;
    ascfid=1; datasize2fid=datasize;
    location1=ineng.tellg();
    datasize2fide=location1-location0e;
    location1=instr.tellg();
    datasize2fids=location1-location0s;
    //printf("trj l0 %d ds %d\n",location0,datasize2fid);
    //printf("eng l0 %d ds %d\n",location0e,datasize2fide);
    //printf("str l0 %d ds %d\n",location0s,datasize2fids);
    //find the total frame number
    totframe=1;
    int tmp=0;
    while(intrj.eof()+instr.eof()+ineng.eof()==0) {
       //cout<<"t "<<intrj.tellg()<<" e "<<ineng.tellg()<<" s "<<instr.tellg()<<endl;
       ReadCPMDFrame();
       //cout<<"t "<<intrj.tellg()<<" e "<<ineng.tellg()<<" s "<<instr.tellg()<<endl;
       if(prp->step > tmp) tmp=prp->step;
       //cout<<totframe<<" "<<intrj.eof()<<" "<<instr.eof()<<" "<<ineng.eof()<<endl;
       //cout<<totframe<<" "<<intrj.good()<<" "<<instr.good()<<" "<<ineng.good()<<endl;
       //cout<<totframe<<" "<<intrj.bad()<<" "<<instr.bad()<<" "<<ineng.bad()<<endl;
       //cout<<totframe<<" "<<intrj.fail()<<" "<<instr.fail()<<" "<<ineng.fail()<<endl;
       totframe++;
    }
    totframe--;
    header.totaccustep=tmp;
    //printf("stlin l0 %d l1 %d size %d l2 %d total frame %d \n",location0,location1,datasize,location2,totframe);
    return 0;
}


int TRJ::init_lmptrj(char *trjfname)
{
    char null[1024];
    cout<<"Analyzing trajectory files : "<<endl;
    strcpy(null,trjfname); strcat(null,"/TRAJECTORY");
    cout<<"TRAJECTORY : "<<null<<endl;
    intrj.open(null,ios::in);
    if(!intrj.is_open()) {
       cout<<" Error: Trj file "<<null<<" cannot be opened."<<endl;
       return filenotopen;
    }
    strcpy(null,trjfname); strcat(null,"/ENERGIES");
    cout<<"ENERGIES   : "<<null<<endl;
    ineng.open(null,ios::in);
    if(!ineng.is_open()) {
       cout<<" Error: ENERGY file "<<null<<" cannot be opened."<<endl;
       return filenotopen;
    }
    strcpy(null,trjfname); strcat(null,"/STRESS");
    cout<<"STRESS     : "<<null<<endl;
    instr.open(null,ios::in);
    if(!instr.is_open()) {
       cout<<" Error: STRESS file "<<null<<" cannot be opened."<<endl;
       return filenotopen;
    }

    rd_header();  //read header

    location0=intrj.tellg();
    location0e=ineng.tellg();
    location0s=instr.tellg();
    init_contentDouble();
    CleanTRJcontentDouble();
    ReadCPMDFrame();
    location1=intrj.tellg();
    datasize=location1-location0;
    ascfid=1; datasize2fid=datasize;
    location1=ineng.tellg();
    datasize2fide=location1-location0e;
    location1=instr.tellg();
    datasize2fids=location1-location0s;
    //printf("trj l0 %d ds %d\n",location0,datasize2fid);
    //printf("eng l0 %d ds %d\n",location0e,datasize2fide);
    //printf("str l0 %d ds %d\n",location0s,datasize2fids);
    //find the total frame number
    totframe=1;
    int tmp=0;
    while(intrj.eof()+instr.eof()+ineng.eof()==0) {
       //cout<<"t "<<intrj.tellg()<<" e "<<ineng.tellg()<<" s "<<instr.tellg()<<endl;
       ReadCPMDFrame();
       //cout<<"t "<<intrj.tellg()<<" e "<<ineng.tellg()<<" s "<<instr.tellg()<<endl;
       if(prp->step > tmp) tmp=prp->step;
       //cout<<totframe<<" "<<intrj.eof()<<" "<<instr.eof()<<" "<<ineng.eof()<<endl;
       //cout<<totframe<<" "<<intrj.good()<<" "<<instr.good()<<" "<<ineng.good()<<endl;
       //cout<<totframe<<" "<<intrj.bad()<<" "<<instr.bad()<<" "<<ineng.bad()<<endl;
       //cout<<totframe<<" "<<intrj.fail()<<" "<<instr.fail()<<" "<<ineng.fail()<<endl;
       totframe++;
    }
    totframe--;
    header.totaccustep=tmp;
    //printf("stlin l0 %d l1 %d size %d l2 %d total frame %d \n",location0,location1,datasize,location2,totframe);
    return 0;
}

int TRJ::rd_header()
{
    switch(itrjformat) {
      case c2trj:
           header.ReadBinHeader(&intrj);
           break;
      case asctrj:
           header.ReadAscHeader(&intrj);
           break;
      case cpmdtrj:
           //find the number of atoms
           header.ReadCPMDHeader(&intrj,&ineng,&instr);
           break;
      case lmptrj:
           //find the number of atoms
           header.ReadLMPHeader(&intrj,&ineng);
           break;
      case usrtrj:
           //find the number of atoms
           header.ReadUSRHeader(&intrj);
           break;
      case none:
      default:
           break;
    }
   
    return 0;
}

int TRJ::wt_header(ostream *outf)
{
    switch(otrjformat) {
      case c2trj:
      case fgtrj:
           header.WriteBinaryheader(outf,0);
           break;
      case asctrj:
           header.WriteASCIIheader(outf,0,0,0);
           break;
      case cpmdtrj:
           break;
      case usrtrj:
           break;
      case lmptrj:
           break;
      case none:
      default:
           break;
    }
    return 0;
}

int TRJ::rd_frame(int iframe)
{
    fcount=(iframe-1);
    intrj.clear();
    switch(itrjformat) {
      case c2trj:
           intrj.seekg(location0+datasize*fcount, ios::beg);
           ReadBinFrame();
           break;
      case usrtrj:
           intrj.seekg(location0+datasize*fcount, ios::beg);
           ReadUSRFrame();
           break;
      case asctrj:
           int i;
           if(iframe>ascfid) {
             intrj.seekg(location0+datasize2fid, ios::beg);
             for(i=ascfid;i<iframe;i++) ReadAscFrame();
             ascfid =iframe;
             datasize2fid  = intrj.tellg();
             datasize2fid -= location0;
           } else {
             intrj.seekg(location0, ios::beg);
             for(i=0;i<iframe;i++) ReadAscFrame();
             ascfid =iframe;
             datasize2fid = intrj.tellg();
             datasize2fid-= location0;
           }
           break;
      case cpmdtrj:
           ineng.clear();
           instr.clear();
           if(iframe>ascfid) {
             intrj.seekg(location0 +datasize2fid , ios::beg);
             ineng.seekg(location0e+datasize2fide, ios::beg);
             instr.seekg(location0s+datasize2fids, ios::beg);
             for(i=ascfid;i<iframe;i++) ReadCPMDFrame();
             ascfid =iframe;
             datasize2fid = intrj.tellg();
             datasize2fid-= location0;
             datasize2fide = ineng.tellg();
             datasize2fide-= location0e;
             datasize2fids = instr.tellg();
             datasize2fids-= location0s;
           } else {
             intrj.seekg(location0 , ios::beg);
             ineng.seekg(location0e, ios::beg);
             instr.seekg(location0s, ios::beg);
             for(i=0;i<iframe;i++) ReadCPMDFrame();
             ascfid =iframe;
             datasize2fid = intrj.tellg();
             datasize2fid-= location0;
             datasize2fide = ineng.tellg();
             datasize2fide-= location0e;
             datasize2fids = instr.tellg();
             datasize2fids-= location0s;
           }
           break;
      case lmptrj:
           break;
      case none:
      default:
           break;
    }

    cell->H2others();
    prp->rho=prp->mass/(prp->V*1E-24*Na); //density in g/cc
    return 0;
}

int TRJ::wt_frame(ostream *outf)
{
    switch(otrjformat) {
      case c2trj:
      case fgtrj:
           WriteBinarycontent(outf,0);
           break;
      case asctrj:
           WriteASCIIcontent(outf,0,0,0);
           break;
      case cpmdtrj:
           break;
      case usrtrj:
           break;
      case lmptrj:
           break;
      case none:
      default:
           break;
    }
    return 0;
}

int TRAJECTORY::init_trj(int ntrj,char (*trjfname)[1024],int iformat,int oformat,CELL *in_cell,ATOM *in_atom,PROPERTY *in_prp)
{
    int filestatus=0;
    ntrjf=ntrj;
    itrjf=iformat;
    otrjf=oformat;
    if(ntrjf==0) return filenotopen;
    if(strj!=NULL) delete [] strj;
    if(iframe!=NULL) delete [] iframe;
    if(fframe!=NULL) delete [] fframe;
    strj=new TRJ [ntrjf];
    iframe=new int [ntrjf];
    fframe=new int [ntrjf];
    int i;
    totframe=0;

    for(i=0;i<ntrjf;i++) {
        strj[i].itrjformat=itrjf;
        strj[i].otrjformat=otrjf;
        strj[i].cell=in_cell;
        strj[i].prp=in_prp;
        strj[i].atom=in_atom;
        filestatus+=strj[i].init_trj(trjfname[i]);
        totframe+=strj[i].totframe;
        fframe[i]=totframe;
    }
    //cout<<filestatus<<" "<<filenotopen<<endl;
    if(filestatus/filenotopen) return filenotopen;

    if(iformat==cpmdtrj) {
       int tmp=strj[0].header.totaccustep;
       int maxxyzfreq=strj[0].header.xyzfreq,maxstrfreq=strj[0].header.strfreq;
       
       for(i=1;i<ntrjf;i++) {
         tmp+=strj[i].header.totaccustep;
         strj[i].header.totaccustep=tmp-strj[i].header.totaccustep;
         if(strj[i].header.xyzfreq>maxxyzfreq) maxxyzfreq=strj[i].header.xyzfreq;
         if(strj[i].header.strfreq>maxstrfreq) maxstrfreq=strj[i].header.strfreq;
       }
       if(maxxyzfreq%maxstrfreq) {
          maxxyzfreq*=maxstrfreq;
       } 
       for(i=0;i<ntrjf;i++) {
          strj[i].header.totaccustep+=maxxyzfreq-1;
          strj[i].header.nxyz=strj[i].header.nvel=maxxyzfreq/strj[i].header.xyzfreq;
          strj[i].header.nstr=maxxyzfreq/strj[i].header.strfreq;
          strj[i].header.neng=maxxyzfreq/strj[i].header.engfreq;
          printf("file %d nxyz %d nvel %d neng %d nstr %d accustep %d\n",i+1,strj[i].header.nxyz,strj[i].header.nvel,strj[i].header.neng,strj[i].header.nstr,strj[i].header.totaccustep);
       }
       strj[0].header.totaccustep=0;
       for(i=0;i<ntrjf;i++) {
          printf("file %d nxyz %d nvel %d neng %d nstr %d accustep %d\n",i+1,strj[i].header.nxyz,strj[i].header.nvel,strj[i].header.neng,strj[i].header.nstr,strj[i].header.totaccustep);
       }
    } else if(iformat==lmptrj) {
    } else if(iformat==usrtrj) {
    } else {
       //correction for lammps trj (remove the 1st frame for 2nd,3rd,..trj)
       for(i=1;i<ntrjf;i++) {
          fframe[i]-= i;
       }
    }

    totframe -= (ntrjf-1);
    if(1) {
      for(i=0;i<ntrjf;i++) {
         printf("stlin trj %d name %s totframe %d last frame id %d\n",i+1,trjfname[i],strj[i].totframe,totframe);
      }
    }
    return itrjf;
}

int TRAJECTORY::rd_frame(int ith)
{
    int i;
    cframe=ith;
    //find current trj file
    for(i=0;i<ntrjf;i++) {
       if( fframe[i]>=ith) break;
    }
    ctrjf=i;
    if(ctrjf==0) { 
       strj[ctrjf].rd_frame(ith); 
    } else { 
       if(itrjf==cpmdtrj) strj[ctrjf].rd_frame(ith-fframe[ctrjf-1]); 
       else if(itrjf==lmptrj) strj[ctrjf].rd_frame(ith-fframe[ctrjf-1]);
       else strj[ctrjf].rd_frame(ith-fframe[ctrjf-1]+1); //adding one for lammps trj
    }
    return 0;
}

int TRAJECTORY::wt_frame(ostream *outf)
{
    strj[ctrjf].wt_frame(outf);
    return 0;
}

int TRAJECTORY::wt_header(ostream *outf)
{
    strj[0].wt_header(outf);
    return 0;
}

