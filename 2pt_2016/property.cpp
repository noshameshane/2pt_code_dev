#include "property.h"

int PROPERTY::zero()
{
    time=T=P=V=rho=Et=Ek=Ep=fvol=Eke=Ete=0;
    Ebond=Eangle=Etorsion=Einversion=Evdw=Eel=Elong=Ehb=0;
    strs[0]=strs[1]=strs[2]=strs[3]=strs[4]=strs[5]=0;
    cmpv[0]=cmpv[1]=cmpv[2]=0;

    return 0;
}
int PROPERTY::out_prp(ostream *outf)
{
    char null[1024];
    *outf<<endl;
    sprintf(null," time(ps)= %14.4f Total_E = %14.4f E_bond  = %14.4f",time,Et,Ebond);
    *outf<<null<<endl;
    sprintf(null," T__(K)  = %14.4f Total_KE= %14.4f E_angle = %14.4f",T,Ek,Eangle);
    *outf<<null<<endl;
    sprintf(null," P__(GPa)= %14.4f Total_PE= %14.4f E_dihed = %14.4f",P,Ep,Etorsion);
    *outf<<null<<endl;
    sprintf(null," V__(A3) = %14.4f E_vdw   = %14.4f E_inv   = %14.4f",V,Evdw,Einversion);
    *outf<<null<<endl;
    sprintf(null," d_(g/cc)= %14.4f E_coul  = %14.4f E_long  = %14.4f",rho,Eel,Elong);
    *outf<<null<<endl;
    return 0;
}
