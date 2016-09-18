#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <iomanip>
#include <ctime>
#include "model.h"
#include "analysis.h"

int main(int argc, char *argv[])
{
  if(argc<2) {
     cout<<" Usage: "<<argv[0]<<" driver.in"<<endl;
     return 0;
  }


  MODEL model;
  ANALYSIS analysis;

  model.ctl.rd_ctl(argv[1]); //read control file
  model.ctl.out_ctl(&cout);
  model.rd_strt(); //read molecular structure
  printf("stlin natom %d nmol %d nbond %d\n",model.natom,model.nmol,model.nbond);

  model.init_trj(); //read trj file
  analysis.doit(&model); //perform analysis

  model.timing.gettimetotal(&cout);


  return 0;
}
