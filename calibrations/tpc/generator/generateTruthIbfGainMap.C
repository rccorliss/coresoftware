//Generate a pair of IBF gain maps from an pinut ADC and ibf-truth file.
// if scale factor or scaleuncertainty are different from zero, we rescale bin by bin with random factors drawn from that gaussian.

void generateTruthIbfGainMap(const char* adcFile, const char *adcName, const char* ionFile, const char* ibfName, const char* primName,
			     const char* outputFile, float scalefactor=1.0, float scaleuncertainty=0.0){

  TRandom *rng=new TRandom();
  //printf("scalefactor=%f, scaleuncertainty=%f\n",scalefactor,scaleuncertainty);
  //return;
  
 //load the adc-per-bin data from the specified file.
  TFile* adcInputFile = TFile::Open(adcFile, "READ");
  TH3* hAdc=nullptr;
  adcInputFile->GetObject(adcName, hAdc);
  if (hAdc == nullptr)  {
    printf("ADC hist %s or file %s not found!\n",adcName, adcFile);
    return false;
  }
  
  //load the ions-per-bin data from the specified file
  TFile* ibfInputFile = TFile::Open(ionFile, "READ");
  TH3* hIbf=nullptr;
  TH3* hPrimaries=nullptr;
  ibfInputFile->GetObject(ibfName,hIbf);
  ibfInputFile->GetObject(primName,hPrimaries);
  if (hIbf == nullptr) {
    printf("IBF hist %s or file %s not found!\n",ibfName, ionFile);
    return false;
  }
  if (hPrimaries == nullptr) {
    printf("IBF hist %s IN file %s not found!\n",primName, ionFile);
    return false;
  }

  TFile *output=TFile::Open(outputFile,"RECREATE");
  //generate 2D histograms with the IBF, and IBF + Primaries, per side:

  int nSections=3;
  TH2 *hFlatTotal[nSections],*hFlatIbf[nSections]; //flat charge sums in north and south
  int zbins,neg,cm,pos;
  zbins=hAdc->GetZaxis()->GetNbins();
  neg=1;
  cm=zbins/2;
  pos=zbins;

  TString suffix[]={"","_negz","_posz"};
  int low[]={neg,neg,cm+1};
  int high[]={pos,cm,pos};

  for (int i=0;i<nSections;i++){
    hPrimaries->GetZaxis()->SetRange(low[i],high[i]);
    hFlatTotal[i]=(TH2*)hPrimaries->Project3D(Form("xy%s",suffix[i].Data()));
    hIbf->GetZaxis()->SetRange(low[i],high[i]);
    hFlatIbf[i]=(TH2*)hIbf->Project3D(Form("xy%s",suffix[i].Data()));
    hFlatTotal[i]->Add(hFlatIbf[i]);
  }
  
  
  //sanity check that the bounds are the same

  TAxis* ax[3] = {nullptr, nullptr, nullptr};
  ax[0] = hPrimaries->GetXaxis();
  ax[1] = hPrimaries->GetYaxis();
  ax[2] = hPrimaries->GetZaxis();

  TAxis* ax2[3] = {nullptr, nullptr, nullptr};
  ax2[0] = hAdc->GetXaxis();
  ax2[1] = hAdc->GetYaxis();
  ax2[2] = hAdc->GetYaxis();


  int nbins[3];
  for (int i = 0; i < 3; i++)//note these are dimensions, not sections.
  {
    nbins[i] = ax[i]->GetNbins();  //number of bins, not counting under and overflow.
    if (nbins[i]!=ax2[i]->GetNbins()){
      printf("Primaries and Adc bins are different in axis %d. (%d vs %d).  Aborting unless in z.\n",i,nbins[i],ax2[i]->GetNbins());
      if (i!=2) return;
      printf("this is a z difference, which we can recover.\n");
    }
  }

  
  //project into 2D per rphi, and divide Ions by Adc to get the Ion/Adc gain.
  TH2* hIonGain[nSections], *hIbfGain[nSections],*hFlatAdc[nSections];

  for (int i=0;i<nSections;i++){
    hIonGain[i]=static_cast<TH2*>(hFlatTotal[i]->Clone(Form("hIonGain%s",suffix[i].Data())));//with primaries
    hIbfGain[i]=static_cast<TH2*>(hFlatIbf[i]->Clone(Form("hIbfGain%s",suffix[i].Data())));//with only ibf
    hAdc->GetZaxis()->SetRange(low[i],high[i]);
    hFlatAdc[i]=static_cast<TH2*>(hAdc->Project3D(Form("xy%s",suffix[i].Data())));

    hIonGain[i]->Divide(hFlatAdc[i]);
    hIbfGain[i]->Divide(hFlatAdc[i]);

    //draw gaussians for every entry, only if the input is not zero.
    if (scaleuncertainty>0.0){
      int nCells=hIonGain[i]->GetNcells();
      for (int j=0;j<nCells;j++){
	float val=hIonGain[i]->GetBinContent(j)*rng->Gaus(scalefactor,scaleuncertainty);
	printf("sanity: binorig=%f,newval=%f\n",hIonGain[i]->GetBinContent(j),val);
	hIonGain[i]->SetBinContent(j,val);
	val=hIbfGain[i]->GetBinContent(j)*rng->Gaus(scalefactor,scaleuncertainty);
	hIbfGain[i]->SetBinContent(j,val);
      }
    } else {
      if (scalefactor!=1.0){
	hIonGain[i]->Scale(scalefactor);
	hIbfGain[i]->Scale(scalefactor);
      }
    }

    hIonGain[i]->Write();
    hIbfGain[i]->Write();
  }
  
  //
  return;
}
