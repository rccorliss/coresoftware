/*
This is a stripped down version of the annular field sim handling code, designed to ONLY generate the rossegger files for use in the other system.  Note that it will only work in the coresoftware directory version, as of apr 3, 2023.

 */



#include "AnnularFieldSim.h"
#include "TTree.h" //this prevents a lazy binding issue and/or is a magic spell.
R__LOAD_LIBRARY(build/.libs/libfieldsim)

void generate_rossegger(int nr, int nphi, int nz, int njobs=1, int thisjob=0){

  TTime now, start;
  start=now=gSystem->Now();
  printf("the time is %lu\n",(unsigned long)now);


  //step 1:  specify the physical parameters, scales, etc of the model, either ALICE or sPHENIXL
  //load the sPHENIX space charge model parameters
  const float tpc_rmin=20.0;
  const float tpc_rmax=78.0;
  float tpc_deltar=tpc_rmax-tpc_rmin;
  const float tpc_z=105.5;
  const float tpc_cmVolt=-400*tpc_z; //V =V_CM-V_RDO -- volts per cm times the length of the drift volume.
  //const float tpc_magField=0.5;//T -- The old value used in carlos's studies.
  //const float tpc_driftVel=4.0*1e6;//cm per s  -- used in carlos's studies
  const float tpc_driftVel=8.0*1e6;//cm per s  -- 2019 nominal value
  //const float tpc_magField=1.4;//T -- 2019 nominal value
  const float tpc_magField=-1.5;//T -- 2019 nominal value in the fieldmap
  const char detgeoname[]="sphenix";

   //step 2: specify the parameters of the field simulation.  Larger numbers of bins will rapidly increase the memory footprint and compute times.
  //there are some ways to mitigate this by setting a small region of interest, or a more parsimonious lookup strategy, specified when AnnularFieldSim() is actually constructed below.
  int nr_roi_min=0;
  int nr_roi=nr;
  int nr_roi_max=nr_roi_min+nr_roi;
  int nphi_roi_min=0;
  int nphi_roi=nphi;
  int nphi_roi_max=nphi_roi_min+nphi_roi;
  int nz_roi_min=0;
  int nz_roi=nz;
  int nz_roi_max=nz_roi_min+nz_roi;

 
 
  //look for the greens functions:
  char lookup_string[200], stripe_string[100];
 sprintf(lookup_string,"ross_phi1_%s_phislice_lookup_r%dxp%dxz%d",detgeoname,nr,nphi,nz);
  char lookupFilename[300];
  //sprintf(lookupFilename,"%s.root",lookup_string);
  if (njobs>1){
    sprintf(stripe_string,"_%dof%d", thisjob,njobs);
  sprintf(lookupFilename,"/sphenix/user/rcorliss/rossegger/%s%s.root",lookup_string,stripe_string); //hardcoded for racf
  TFile *fileptr=TFile::Open(lookupFilename,"READ");

  if (fileptr){
    printf("file %s already exists.  Nothing to do.\n", lookupFilename);
	return;
  } else{

  

    //step 4:  create the fieldsim object.  different choices of the last few arguments will change how it builds the lookup table spatially, and how it loads the spacecharge.  The various start-up steps are exposed here so they can be timed in the macro.
    //To set rossegger vs trivial greens functions you either call or do not call load_rossegger().
    start=now;
    AnnularFieldSim *tpc=
      new  AnnularFieldSim(tpc_rmin,tpc_rmax,tpc_z,
			   nr, nr_roi_min,nr_roi_max,1,2,
			   nphi,nphi_roi_min, nphi_roi_max,1,2,
			   nz, nz_roi_min, nz_roi_max,1,2,
			   tpc_driftVel, AnnularFieldSim::PhiSlice, AnnularFieldSim::FromFile);
    tpc->UpdateEveryN(10);//show reports every 10%.
  
    now=gSystem->Now();
    printf("created sim obj.  the dtime is %lu\n",(unsigned long)(now-start));
    start=now;

    tpc->load_rossegger();
    now=gSystem->Now();
    printf("loaded rossegger greens functions. the dtime is %lu\n",(unsigned long)(now-start));
    start=now;
    tpc->populate_phislice_lookup(njobs,thisjob);
    tpc->save_phislice_lookup(lookupFilename);
    printf("saved greens functions to %s\n",lookupFilename);
 
  }

  
  
  return;
}
