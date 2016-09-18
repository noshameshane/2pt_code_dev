#include "control.h"

int CONTROL::rd_ctl(char *ctlf)
{
    cout<<"Reading control file "<<ctlf<<endl;
    ifstream inctl(ctlf,ios::in);
    if(!inctl.is_open()) return filenotopen;

    char opt[1024],*charptr;
    char delims[] = {" ;,\t"};
    int i;
    inctl>>opt;
    while(!inctl.eof()) {
        if(opt[0]=='#'||opt[0]=='*'||opt[0]=='!'||opt[0]=='/') { inctl.getline(opt,1024); }
        else if (strcmp(opt,"IN_LMPDATA")==0) {inctl>>in_strt; in_strt_flag=strtlmp; }
        else if (strcmp(opt,"IN_BGF")==0) {inctl>>in_strt; in_strt_flag=strtbgf; }
        else if (strcmp(opt,"IN_GRO")==0) {inctl>>in_strt; in_strt_flag=strtgro; }
        else if (strcmp(opt,"IN_G96")==0) {inctl>>in_strt; in_strt_flag=strtg96; }
        else if (strcmp(opt,"IN_CPMDIN")==0) {inctl>>in_strt; in_strt_flag=strtcpmd; }
        else if (strcmp(opt,"IN_C2TRJ")==0||strcmp(opt,"IN_ASCTRJ")==0||strcmp(opt,"IN_CPMDTRJ_PATH")==0||strcmp(opt,"IN_USRTRJ")==0) { 
          if (strcmp(opt,"IN_C2TRJ")==0) in_trj_flag=c2trj; 
          else if (strcmp(opt,"IN_ASCTRJ")==0)  in_trj_flag=asctrj;
          else if (strcmp(opt,"IN_USRTRJ")==0)  in_trj_flag=usrtrj;
          else if (strcmp(opt,"IN_CPMDTRJ_PATH")==0)  in_trj_flag=cpmdtrj;
          inctl.getline(opt,1024);
          char null[1024]; strcpy(null,opt);
          in_trj_n=0;
          charptr = strtok(opt, delims);
          while(charptr!=NULL) {
             in_trj_n++;
             charptr = strtok(NULL, delims);
          }
          charptr = strtok(null, delims);
          in_trj=new char [in_trj_n][1024];
          for(i=0;i<in_trj_n;i++) {
             strcpy(in_trj[i],charptr);
             charptr = strtok(NULL, delims);
             //printf("stlin ntrj %d name %s\n",i,in_trj[i]);
          }
        } 
        else if (strcmp(opt,"IN_GROUPFILE")==0) { inctl>>in_grp; in_grp_flag=1; }
        else if (strcmp(opt,"ANALYSIS_FRAME_INITIAL")==0) { inctl>>ana_iframe; }
        else if (strcmp(opt,"ANALYSIS_FRAME_FINAL")==0) { inctl>>ana_fframe;  }
        else if (strcmp(opt,"ANALYSIS_FRAME_STEP")==0) { inctl>>ana_sframe;  }
        else if (strcmp(opt,"ANALYSIS_VAC_CORLENGTH")==0) { inctl>>ana_vac_corlen; ana_flag=ana_vac_flag=1; }
        else if (strcmp(opt,"ANALYSIS_VAC_MEMORYMB")==0) { inctl>>ana_vac_mem; }
        else if (strcmp(opt,"ANALYSIS_VAC_2PT")==0) { inctl>>ana_vac_2pt; if(ana_vac_2pt>5||ana_vac_2pt<0) ana_vac_2pt=0; }
        else if (strcmp(opt,"ANALYSIS_VAC_FIXED_DF")==0) { inctl>>ana_vac_const; }
        else if (strcmp(opt,"ANALYSIS_VAC_ROTN_SYMMETRY")==0) { inctl>>ana_vac_rotsym; }
        else if (strcmp(opt,"ANALYSIS_VAC_LINEAR_MOL")==0) { inctl>>ana_vac_linear; }
        else if (strcmp(opt,"ANALYSIS_VAC_INC_CNT")==0) { inctl>>ana_vac_cnt; }
        else {
              cout<<" Error in rd_ctl "<<opt<<endl;
              return unregkey;
        }
        inctl>>opt;
    }
    //consistency settings
    if(out_grp_flag>0) { in_grp_flag=1; strcpy(in_grp,out_grp); } //use out_grp as in_grp in the same session
    return 0;
}

int CONTROL::out_ctl(ostream *outf)
{
    char null[1024];
    int i;
    if(in_strt_flag==strtlmp) {
       sprintf(null,"%-40s %s","IN_LMPDATA",in_strt);
       *outf<<null<<endl;
    } else if(in_strt_flag==strtbgf) {
       sprintf(null,"%-40s %s","IN_BGF",in_strt);
       *outf<<null<<endl;
    } else if(in_strt_flag==strtgro) {
       sprintf(null,"%-40s %s","IN_GRO",in_strt);
       *outf<<null<<endl;
    } else if(in_strt_flag==strtg96) {
       sprintf(null,"%-40s %s","IN_G96",in_strt);
       *outf<<null<<endl;
    } else if(in_strt_flag==strtcpmd) {
       sprintf(null,"%-40s %s","IN_CPMDIN",in_strt);
       *outf<<null<<endl;
    }
    if(in_ff_flag==ffdreid) {
       sprintf(null,"%-40s %s","IN_FORCEFIELD",in_ff);
       *outf<<null<<endl;
    }
    if(in_trj_flag==c2trj) {
       sprintf(null,"%-40s %s","IN_C2TRJ",in_trj[0]);
       for(i=1;i<in_trj_n;i++) sprintf(null,"%s %s",null,in_trj[i]);
       *outf<<null<<endl;
    } else if(in_trj_flag==asctrj) {
       sprintf(null,"%-40s %s","IN_ASCTRJ",in_trj[0]);
       for(i=1;i<in_trj_n;i++) sprintf(null,"%s %s",null,in_trj[i]);
       *outf<<null<<endl;
    } else if(in_trj_flag==usrtrj) {
       sprintf(null,"%-40s %s","IN_USRTRJ",in_trj[0]);
       for(i=1;i<in_trj_n;i++) sprintf(null,"%s %s",null,in_trj[i]);
       *outf<<null<<endl;
    } else if(in_trj_flag==cpmdtrj) {
       sprintf(null,"%-40s %s","IN_CPMDTRJ_PATH",in_trj[0]);
       for(i=1;i<in_trj_n;i++) sprintf(null,"%s %s",null,in_trj[i]);
       *outf<<null<<endl;
    }

    if(in_grp_flag) {
       sprintf(null,"%-40s %s","IN_GROUPFILE",in_grp);
       *outf<<null<<endl;
    }

    if(ana_flag) {
      sprintf(null,"%-40s %d","ANALYSIS_FRAME_INITIAL",ana_iframe);
      *outf<<null<<endl;
      sprintf(null,"%-40s %d","ANALYSIS_FRAME_FINAL",ana_fframe);
      *outf<<null<<endl;
      sprintf(null,"%-40s %d","ANALYSIS_FRAME_STEP",ana_sframe);
      *outf<<null<<endl;

      if(ana_vac_flag) {
         sprintf(null,"%-40s %f","ANALYSIS_VAC_CORLENGTH",ana_vac_corlen);
         *outf<<null<<endl;
         sprintf(null,"%-40s %f","ANALYSIS_VAC_MEMORYMB",ana_vac_mem);
         *outf<<null<<endl;
         sprintf(null,"%-40s %d","ANALYSIS_VAC_2PT",ana_vac_2pt);
         *outf<<null<<endl;
         sprintf(null,"%-40s %s","ANALYSIS_VAC_FIXED_DF",ana_vac_const);
         *outf<<null<<endl;
         sprintf(null,"%-40s %s","ANALYSIS_VAC_ROTN_SYMMETRY",ana_vac_rotsym);
         *outf<<null<<endl;
         sprintf(null,"%-40s %s","ANALYSIS_VAC_LINEAR_MOL",ana_vac_linear);
         *outf<<null<<endl;
         sprintf(null,"%-40s %s","ANALYSIS_VAC_INC_CNT",ana_vac_cnt);
         *outf<<null<<endl;
      }
    }
    return 0;
}

int CONTROL::printsameline(char *c1,int num,char*c2)
{
    if (isatty(1)){
       int i,len,f;
       len=1;f=10;
       while(num/f) { f*=10; len++; };
       len+=strlen(c1)+strlen(c2);
       for(i=0;i<len+2;i++) printf("\b");
       printf("%s %d %s",c1,num,c2);
       fflush(stdout);
    }
    return 0;
}

int CONTROL::printsameline(char *c1,int num1,char*c2, int num2)
{
    if (isatty(1)){
       int i,len,f;
       len=1;f=10;
       while(num1/f) { f*=10; len++; };
       f=10;
       while(num2/f) { f*=10; len++; };
       len+=strlen(c1)+strlen(c2);
       for(i=0;i<len+4;i++) printf("\b");
       printf("%s %d %s %d",c1,num1,c2,num2);
       fflush(stdout);
    }
    return 0;
}

