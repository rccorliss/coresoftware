

#include "AnnularFieldSim.h"
#include "TTree.h" //this prevents a lazy binding issue and/or is a magic spell.
#include "TCanvas.h" //this prevents a lazy binding issue and/or is a magic spell.

// cppcheck-suppress unknownMacro
R__LOAD_LIBRARY(libfieldsim.so)

char field_string[200];
char lookup_string[200];

AnnularFieldSim *SetupDefaultSphenixTpc(bool twinMe=false, bool useSpacecharge=true, float xshift=0, float yshift=0, float zshift=0);
AnnularFieldSim *SetupDigitalCurrentSphenixTpc(bool twinMe=false, bool useSpacecharge=true);
void TestSpotDistortion(AnnularFieldSim *t); //
void SurveyFiles(TFileCollection* filelist);

struct DistortionMapParameters {
  //files
  TString outputName;

  //input adc/ibf file and hist name:
  TString inputName;
  TString gainName;
  TString gainHistName[2];

  TString ibfName;
  TString primName;

  // external field maps
  TString bfieldName;
  TString bfieldTreeName;
  TString efieldName;
  TString efieldTreeName;

  //B field positioning
  float xshift;
  float yshift;
  float zshift;

  // flags and parameters
  bool usesChargeDensity;
  bool hasSpacecharge;
  bool isAdc;
  int nSteps;
};

void tidy_generate_distortion_map(){
//I am modifying this entire code, from the generate_distortion_map.C macro, to make it much tidier and simpler, so that the parameters are all in one place, and modifications are hopefully easier to see.

DistortionMapParameters params;
params.inputName="/sphenix/user/shulga/Work/IBF/DistortionMap/Files/Summary_hist_mdc2_UseFieldMaps_AA_event_0_bX10556072.root";
params.gainName="no_gain";
params.gainHistName[0]="hIbfGain_posz"; //0=north=positive z.  
params.gainHistName[1]="hIbfGain_negz";
params.outputName="fill_me_in";
params.ibfName="_h_SC_ibf_0";
params.primName="_h_SC_prim_0";
params.hasSpacecharge=true;
params.isAdc=false;
params.nSteps=450;
params.xshift=0;
params.yshift=0;
params.zshift=0;
params.bfieldName="/sphenix/user/rcorliss/rossegger/sphenix3dmaprhophiz.root";
params.bfieldTreeName="fieldmap";
params.efieldName="/sphenix/user/rcorliss/field/externalEfield.ttree.root";
params.efieldTreeName="fTree"; 
params.usesChargeDensity=false; //true if source hists contain charge density per bin.  False if hists are charge per bin.
params.tpc_chargescale=1.6e-19;//Coulombs per bin unit.
params.spacecharge_cm_per_axis_unit=0.1;//cm per histogram axis unit (mm), matching the MDC2 sample from Evgeny.
params.hasTwin=true; //this flag prompts the code to build both a positive-half and a negative-half for the TPC, reusing as much of the calculations as possible.  It is more efficient to 'twin' one half of the TPC than to recalculate/store the greens functions for both.

std::vector<DistortionMapParameters> paramsets;
params.outputName="/direct/star+u/rcorliss/sphenix/coresoftware_official/brokentpcE/stock_0_0_0";

paramsets.push_back(params);
params.efieldName="/direct/star+u/rcorliss/sphenix/coresoftware_official/brokentpcE/nominal.root.500000.root";
params.efieldTreeName="field_ntuple";
params.outputName="/direct/star+u/rcorliss/sphenix/coresoftware_official/brokentpcE/nominal_0_0_0";
paramsets.push_back(params);
params.efieldName="/direct/star+u/rcorliss/sphenix/coresoftware_official/brokentpcE/ifcmod.root.500000.root";
params.efieldTreeName="field_ntuple";
params.outputName="/direct/star+u/rcorliss/sphenix/coresoftware_official/brokentpcE/ifcmod_0_0_0";
paramsets.push_back(params);
params.efieldName="/direct/star+u/rcorliss/sphenix/coresoftware_official/brokentpcE/stpmod.221.root.180000.root";
params.efieldTreeName="field_ntuple";
params.outputName="/direct/star+u/rcorliss/sphenix/coresoftware_official/brokentpcE/stpmod_0_0_0";
paramsets.push_back(params);
params.efieldName="/direct/star+u/rcorliss/sphenix/coresoftware_official/brokentpcE/stp_and_ifcmod.221.root.180000.root";
params.efieldTreeName="field_ntuple";
params.outputName="/direct/star+u/rcorliss/sphenix/coresoftware_official/brokentpcE/stp_and_ifcmod_0_0_0";
paramsets.push_back(params);


//now build the time-consuming part:
//(making the params that are used explicit.  Those are the ones that won't change if we repeat the below.)
AnnularFieldSim *tpc;
tpc=SetupDefaultSphenixTpc(params.hasTwin,params.hasSpacecharge);//sets the dimensions and rules.

//and pick a consistent location to plot the fieldslices about:
TVector3 pos=0.5*(tpc->GetOuterEdge()+tpc->GetInnerEdge());;
pos.SetPhi(3.14159);
  
//then try to load the greens functions: 
LoadGreensFunctions(tpc);
//(because I can't imagine wanting to _change_ the Greens functions.)


//and now we update Fields and Spacecharge settings for each params in the vector, and generate distortion maps + field maps for each.
for (int i=0; i<paramsets.size(); i++){
  params=paramsets[i];
  if (paramsets.size()>1){
    printf("paramset %d\n", i);
  }
  //try to load the external fields:
  UpdateFields(tpc, params);
  printf("set fields.\n");
  //then load the spacecharge histograms from file: 
  UpdateSpaceCharge(tpc, params);
  printf("set spacecharge.\n");

  //generate the maps and the field slices:
  tpc->GenerateSeparateDistortionMaps(std::string(params.outputName.Data()),params.nSteps,1,1,1,1,false);
  printf("distortions mapped.\n");
  tpc->PlotFieldSlices(params.outputName.Data(),pos, 'E'); //plot the electric field
  tpc->PlotFieldSlices(params.outputName.Data(),pos,'B'); //plot the magnetic field
  printf("fieldslices plotted.\n");
}

printf("done.\n");
return;
}

void UpdateSpaceCharge(AnnularFieldSim *tpc, DistortionMapParameters params){
//load the spacecharge histograms from the file in the current parameters set,
//using the scaling parameters, density rules, and twinning rules.
//if that goes well, also recalculate the fieldmap from the chargemap.

  TFile *infile=TFile::Open(params.inputName.Data(),"READ");
  TFile *gainfile; //if we have one.

  //the total charge is prim + IBF
  //if we are doing ADCs, though, we only read the one.
  TH3* hCharge=(TH3*)(infile->Get(params.ibfName));
  if (!params.isAdc){
    hCharge->Add((TH3*)(infile->Get(params.primName)));
  }   
    //hCharge->Scale(70);//Scaleing the histogram spacecharge by 100 times

  TString chargestring;
	       
  //load the spacecharge into the distortion map generator:
  //the signature is: void load_spacecharge(TH3F *hist, float zoffset, float chargescale, float cmscale, bool isChargeDensity);
  if (!params.isAdc){
    chargestring=Form("%s:(%s+%s)",params.inputName.Data(),params.ibfName.Data(),params.primName.Data());
    tpc->load_spacecharge(hCharge,0,tpc_chargescale,spacecharge_cm_per_axis_unit, usesChargeDensity, chargestring.Data());
    if (params.hasTwin) tpc->twin->load_spacecharge(hCharge,0,tpc_chargescale,spacecharge_cm_per_axis_unit, usesChargeDensity);
  }
  if (params.isAdc){ //load digital current using the scaling:
    gainfile=TFile::Open(params.gainName,"READ");
    TH2* hGain[2];
    hGain[0]=(TH2*)(gainfile->Get(params.gainHistName[0]));
    chargestring=Form("%s:(dc:%s g:%s:%s)",params.inputName.Data(),params.ibfName.Data(),params.gainName.Data(),params.gainHistName[0].Data());
    tpc->load_digital_current(hCharge,hGain[0],tpc_chargescale,spacecharge_cm_per_axis_unit,chargestring.Data());
    if (params.hasTwin) {
      hGain[1]=(TH2*)(gainfile->Get(params.gainHistName[1]));
      tpc->twin->load_digital_current(hCharge,hGain[1],tpc_chargescale,spacecharge_cm_per_axis_unit,chargestring.Data());
    }
  }
  //build the electric fieldmap from the chargemap
  tpc->populate_fieldmap();
  if (params.hasTwin)  tpc->twin->populate_fieldmap();
}

void UpdateFields(AnnularFieldSim *tpc, DistortionMapParameters params){
  //load the external E field:
  tpc->loadEfield(params.efieldName.Data(),params.efieldTreeName.Data());
  if (params.hasTwin) tpc->twin->loadEfield(params.efieldName.Data(),params.efieldTreeName.Data(),-1);//final '-1' tells it to flip z and the field z coordinate. r and phi won't change.

  //load the external B field:
  tpc->load3dBfield(params.bfieldName.Data(),params.bfieldTreeName.Data(),1,-1.4/1.5, params.xshift, params.yshift, params.zshift);
  if (params.hasTwin) tpc->twin->load3dBfield(params.bfieldName.Data(),params.bfieldTreeName.Data(),1,-1.4/1.5, params.xshift, params.yshift, -params.zshift);//z shift is inverted for the twin.
// really need to think about that z shift... but not a problem for the 0,0,0 case.
return;
}

void LoadGreensFunctions(AnnularFieldSim *tpc){
  //load the greens functions:
  const float tpc_rmin=20.0;
  const float tpc_rmax=78.0;
  const float tpc_z=105.5;
  const char detgeoname[]="sphenix";

  sprintf(lookup_string,"ross_phi1_%s_phislice_lookup_r%dxp%dxz%d",detgeoname,
    tpc->GetNRbins(),tpc->GetNPhibins(),tpc->GetNZbins());
  char lookupFilename[200];
  sprintf(lookupFilename,"%s.root",lookup_string);
  TFile *fileptr=TFile::Open(lookupFilename,"READ");

  if (!fileptr){ //generate the lookuptable if it's not where we expect it to be on disk.
    tpc->load_rossegger();
    printf("loaded rossegger greens functions.\n");
    tpc->populate_lookup();
    tpc->save_phislice_lookup(lookupFilename);
  } else{ //load it from a file
    fileptr->Close();
    tpc->load_phislice_lookup(lookupFilename);
  }

  printf("populated lookup.\n");


  if (params.hasTwin==false) return; //no twin to set up.
  //borrow the greens functions:
  tpc->twin->borrow_rossegger(tpc->green,tpc_z);//use the original's green's functions, shift our internal coordinates by tpc_z when querying those functions.
  tpc->twin->borrow_epartial_from(tpc,tpc_z);//use the original's epartial.  Note that those values ought to be symmetric about z, and since our boundary conditions are translated along with our coordinates, they're completely unchanged.  (they're on top of the static solution from the fieldcage)
  printf("twin greens set up.\n");
  return;
}

AnnularFieldSim *SetupDefaultSphenixTpc(bool twinMe, bool useSpacecharge, float xshift, float yshift, float zshift){
// Set up an AnnularFieldSim object for the sPHENIX TPC with default parameters (uses spacecharge, not ADC).



  //step1:  specify the sPHENIX space charge model parameters
  const float tpc_rmin=20.0;
  const float tpc_rmax=78.0;
  float tpc_deltar=tpc_rmax-tpc_rmin;
  const float tpc_z=105.5;
  const float tpc_cmVolt=-400*tpc_z; //V =V_CM-V_RDO -- volts per cm times the length of the drift volume.
  //const float tpc_magField=0.5;//T -- The old value used in carlos's studies.
  //const float tpc_driftVel=4.0*1e6;//cm per s  -- used in carlos's studies
  const float tpc_driftVel=8.0*1e6;//cm per s  -- 2019 nominal value
  const float tpc_magField=1.4;//T -- 2019 nominal value
  const char detgeoname[]="sphenix";


  //step 2: specify the parameters of the field simulation.  Larger numbers of
  // bins will rapidly increase the memory footprint and compute times.  There
  // are some ways to mitigate this by setting a small region of interest, or a
  // more parsimonious lookup strategy, specified when AnnularFieldSim() is
  // actually constructed below.
  int nr=26;//10;//24;//159;//159 nominal
  int nr_roi_min=0;
  int nr_roi=nr;//10;
  int nr_roi_max=nr_roi_min+nr_roi;
  int nphi=40;//38;//360;//360 nominal
  int nphi_roi_min=0;
  int nphi_roi=nphi;//38;
  int nphi_roi_max=nphi_roi_min+nphi_roi;
  int nz=40;//62;//62 nominal
  int nz_roi_min=0;
  int nz_roi=nz;
  int nz_roi_max=nz_roi_min+nz_roi;

  bool realB=true;
  bool realE=true;
  

  
  //step 3:  create the fieldsim object.  different choices of the last few arguments will change how it
  // builds the lookup table spatially, and how it loads the spacecharge.  The various start-up steps
  // are exposed here so they can be timed in the macro.
  
  AnnularFieldSim *tpc;
  if (useSpacecharge){
    tpc=    new  AnnularFieldSim(tpc_rmin,tpc_rmax,tpc_z,
			 nr, nr_roi_min,nr_roi_max,1,2,
			 nphi,nphi_roi_min, nphi_roi_max,1,2,
			 nz, nz_roi_min, nz_roi_max,1,2,
				 tpc_driftVel, AnnularFieldSim::PhiSlice, AnnularFieldSim::FromFile);
  }else{
    tpc=    new  AnnularFieldSim(tpc_rmin,tpc_rmax,tpc_z,
			 nr, nr_roi_min,nr_roi_max,1,2,
			 nphi,nphi_roi_min, nphi_roi_max,1,2,
			 nz, nz_roi_min, nz_roi_max,1,2,
				 tpc_driftVel, AnnularFieldSim::PhiSlice, AnnularFieldSim::NoSpacecharge);
  }
    
  tpc->UpdateEveryN(10);//show reports every 10%.






    
  //make our twin:
  if(twinMe){
    AnnularFieldSim *twin;
      if (useSpacecharge){
	twin=      new  AnnularFieldSim(tpc_rmin,tpc_rmax,-tpc_z,
			   nr, nr_roi_min,nr_roi_max,1,2,
			   nphi,nphi_roi_min, nphi_roi_max,1,2,
			   nz, nz_roi_min, nz_roi_max,1,2,
			   tpc_driftVel, AnnularFieldSim::PhiSlice, AnnularFieldSim::FromFile);
      } else{
	twin=      new  AnnularFieldSim(tpc_rmin,tpc_rmax,-tpc_z,
			   nr, nr_roi_min,nr_roi_max,1,2,
			   nphi,nphi_roi_min, nphi_roi_max,1,2,
			   nz, nz_roi_min, nz_roi_max,1,2,
			   tpc_driftVel, AnnularFieldSim::PhiSlice, AnnularFieldSim::NoSpacecharge);
      }    twin->UpdateEveryN(10);//show reports every 10%.



    tpc->set_twin(twin);
  }

  return tpc;
}
  

AnnularFieldSim *SetupDigitalCurrentSphenixTpc(bool twinMe, bool useSpacecharge){
//build a TPC fieldsim set up for modeling distortions from digital current (uses a gain map)


  //step1:  specify the sPHENIX space charge model parameters
  const float tpc_rmin=20.0;
  const float tpc_rmax=78.0;
  float tpc_deltar=tpc_rmax-tpc_rmin;
  const float tpc_z=105.5;
  const float tpc_cmVolt=-400*tpc_z; //V =V_CM-V_RDO -- volts per cm times the length of the drift volume.
  //const float tpc_magField=0.5;//T -- The old value used in carlos's studies.
  //const float tpc_driftVel=4.0*1e6;//cm per s  -- used in carlos's studies
  const float tpc_driftVel=8.0*1e6;//cm per s  -- 2019 nominal value
  const float tpc_magField=1.4;//T -- 2019 nominal value
  const char detgeoname[]="sphenix";
  
  //step 2: specify the parameters of the field simulation.  Larger numbers of
  // bins will rapidly increase the memory footprint and compute times.  There
  // are some ways to mitigate this by setting a small region of interest, or a
  // more parsimonious lookup strategy, specified when AnnularFieldSim() is
  // actually constructed below.
  
  int nr=8;//10;//24;//159;//159 nominal
  int nr_roi_min=0;
  int nr_roi=nr;//10;
  int nr_roi_max=nr_roi_min+nr_roi;
  int nphi=3*12;//38;//360;//360 nominal
  int nphi_roi_min=0;
  int nphi_roi=nphi;//38;
  int nphi_roi_max=nphi_roi_min+nphi_roi;
  int nz=40;//62;//62 nominal
  int nz_roi_min=0;
  int nz_roi=nz;
  int nz_roi_max=nz_roi_min+nz_roi;

  bool realB=true;
  bool realE=true;
  

  //step 3:  create the fieldsim object.  different choices of the last few arguments will change how it
  // builds the lookup table spatially, and how it loads the spacecharge.  The various start-up steps
  // are exposed here so they can be timed in the macro.
  AnnularFieldSim *tpc;
  if (useSpacecharge){
    tpc=    new  AnnularFieldSim(tpc_rmin,tpc_rmax,tpc_z,
			 nr, nr_roi_min,nr_roi_max,1,2,
			 nphi,nphi_roi_min, nphi_roi_max,1,2,
			 nz, nz_roi_min, nz_roi_max,1,2,
				 tpc_driftVel, AnnularFieldSim::PhiSlice, AnnularFieldSim::FromFile);
  }else{
    tpc=    new  AnnularFieldSim(tpc_rmin,tpc_rmax,tpc_z,
			 nr, nr_roi_min,nr_roi_max,1,2,
			 nphi,nphi_roi_min, nphi_roi_max,1,2,
			 nz, nz_roi_min, nz_roi_max,1,2,
				 tpc_driftVel, AnnularFieldSim::PhiSlice, AnnularFieldSim::NoSpacecharge);
  }
    
  tpc->UpdateEveryN(10);//show reports every 10%.




    
  //make our twin:
  if(twinMe){
    AnnularFieldSim *twin;
      if (useSpacecharge){
	twin=      new  AnnularFieldSim(tpc_rmin,tpc_rmax,-tpc_z,
			   nr, nr_roi_min,nr_roi_max,1,2,
			   nphi,nphi_roi_min, nphi_roi_max,1,2,
			   nz, nz_roi_min, nz_roi_max,1,2,
			   tpc_driftVel, AnnularFieldSim::PhiSlice, AnnularFieldSim::FromFile);
      } else{
	twin=      new  AnnularFieldSim(tpc_rmin,tpc_rmax,-tpc_z,
			   nr, nr_roi_min,nr_roi_max,1,2,
			   nphi,nphi_roi_min, nphi_roi_max,1,2,
			   nz, nz_roi_min, nz_roi_max,1,2,
			   tpc_driftVel, AnnularFieldSim::PhiSlice, AnnularFieldSim::NoSpacecharge);
      }    twin->UpdateEveryN(10);//show reports every 10%.

    tpc->set_twin(twin);
  }

  return tpc;
}
  
void  SurveyFiles(TFileCollection *filelist){
  //list all the files in the filelist, all the histograms in each file,
  //and sum up the total charge in all the Primaries histograms.
   TFile *infile;
 
  TString sourcefilename;
  TString outputfilename;
 //run a check of the files we'll be looking at:
  float integral_sum=0;
  int nhists=0;
  for (int i=0;i<filelist->GetNFiles();i++){
   //for each file, find all histograms in that file.
    sourcefilename=((TFileInfo*)(filelist->GetList()->At(i)))->GetCurrentUrl()->GetFile();//gross
    printf("file %d: %s\n", i, sourcefilename.Data());
    infile=TFile::Open(sourcefilename.Data(),"READ");
    TList *keys=infile->GetListOfKeys();
    //keys->Print();
    int nKeys=infile->GetNkeys();

    for (int j=0;j<nKeys;j++){
      TObject *tobj=infile->Get(keys->At(j)->GetName());
      //if this isn't a 3d histogram, skip it:
      bool isHist=tobj->InheritsFrom("TH3");
      if (!isHist) continue;
      TString objname=tobj->GetName();
      printf(" hist %s ",objname.Data());
      if (objname.Contains("IBF")) {
	printf(" is IBF only.\n");
	continue; //this is an IBF map we don't want.
      }
      float integral=((TH3D*)tobj)->Integral();
      integral_sum+=integral;
      nhists+=1;
	printf(" will be used.  Total Q=%3.3E (ave=%3.3E)\n",integral,integral_sum/nhists);
      //assume this histogram is a charge map.
      //load just the averages:
 
      //break; //rcc temp -- uncomment this to process one hist per file.
      //if (i>maxmaps) return;
    }
      infile->Close();
  }
  return;
}
