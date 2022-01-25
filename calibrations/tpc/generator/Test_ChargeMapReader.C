//This macro loads a charge hist, feeds it into ChargeMapReader, and the compares the resulting fixed-binsize hist to the original.


#include "ChargeMapReader.h"
#include "TH3F.h"
#include "TAxis.h" //this prevents a lazy binding issue and/or is a magic spell.
#include "TCanvas.h" //this prevents a lazy binding issue and/or is a magic spell.

R__LOAD_LIBRARY(build/.libs/libfieldsim)

void CompareHistograms(TH3* hLead, TH3* hFollow, TH1* hOutput){
  TAxis *ax[3]={nullptr,nullptr,nullptr};
  ax[0]=hLead->GetXaxis();
  ax[1]=hLead->GetYaxis();
  ax[2]=hLead->GetZaxis();

  int nbins[3];
  for (int i=0;i<3;i++){
    nbins[i]=ax[i]->GetNbins();//number of bins, not counting under and overflow.
  }

    //   0     1     2   ...   n-1    n    n+1
  // under|first|   ..|  ..  |  .. |last| over
  int i[3];
  int a;
  float low[3],high[3],mid[3];
  float dphi,dr,dz; //bin widths in each dimension.  Got too confusing to make these an array.
  for ( i[0]=1;i[0]<=nbins[0];i[0]++){//phi
    a=0;
    low[a]=ax[a]->GetBinLowEdge(i[a]);
    high[a]=ax[a]->GetBinUpEdge(i[a]);
    mid[a]=(high[a]+low[a])*0.5;
    for ( i[1]=1;i[1]<=nbins[1];i[1]++){//r
      a=1;
      low[a]=ax[a]->GetBinLowEdge(i[a]);
      high[a]=ax[a]->GetBinUpEdge(i[a]);
      mid[a]=(high[a]+low[a])*0.5;
      float rphiterm=dphi*(low[1]+0.5*dr)*dr;
      for ( i[2]=1;i[2]<=nbins[2];i[2]++){//z
	a=1;
	low[a]=ax[a]->GetBinLowEdge(i[a]);
	high[a]=ax[a]->GetBinUpEdge(i[a]);
	mid[a]=(high[a]+low[a])*0.5;
	float leadValue=hLead->GetBinContent(hLead->FindBin(mid[0],mid[1],mid[2]));//take the exact center of the bin.
	float followValue;
	if (ChargeMapReader::CanInterpolateAt(mid[0],mid[1],mid[2],hFollow)){
	  followValue=hFollow->Interpolate(mid[0],mid[1],mid[2]);//interpolate to match the center of the leader bin.
	} else { //if interpolation would fail, fall back to the next best thing:
	  followValue=hFollow->GetBinContent(hLead->FindBin(mid[0],mid[1],mid[2]));//take the exact center of the bin.
	}
	hOutput->Fill((followValue/leadValue-1));
      }
    }
  }
  return;
}
void Test_ChargeMapReader(){

  TFile *originalFile=TFile::Open("/sphenix/user/shulga/Work/IBF/DistortionMap/Files/Summary_hist_mdc2_UseFieldMaps_AA_event_0_bX10556072.root","READ");
  TH3* hOriginalCharge=(TH3*)(originalFile->Get("_h_SC_ibf_0")); //this is only the IBF, so don't be surprised if it looks flatter than expected.

  //get the bounds of the original histogram:
  
  TAxis *ax[3]={nullptr,nullptr,nullptr};
  ax[0]=hOriginalCharge->GetXaxis();
  ax[1]=hOriginalCharge->GetYaxis();
  ax[2]=hOriginalCharge->GetZaxis();

  int nbins[3];
  float low[3],high[3];
  //bin conventions for TH's:
  //   0     1     2   ...   n-1    n    n+1
  // under|first|   ..|  ..  |  .. |last| over
  for (int i=0;i<3;i++){
    nbins[i]=ax[i]->GetNbins();//number of bins, not counting under and overflow.
    low[i]=ax[i]->GetBinLowEdge(1);
    high[i]=ax[i]->GetBinUpEdge(nbins[i]);
  }
  printf("Original bins: %d x %d x %d.  New bins = %d x %d x %d\n", nbins[0],nbins[1],nbins[2], nbins[0],nbins[1],nbins[2]);
  //chargemapreader takes parameters in r,phi,z because that's sane.
  ChargeMapReader *r=new ChargeMapReader(nbins[1],low[1],high[1],nbins[0],low[0],high[0],nbins[2],low[2],high[2]);
  printf("Macro requesting r to readSourceCharge\n");
  r->ReadSourceCharge(hOriginalCharge);
  printf("Returned to macro\n");

  //get the density map from the reader, and populate the resampled charge as well:
  TH3* hOriginalDensity=r->GetDensityHistogram();
  //but it fills output histograms in phi,r,z, because that's convention.
  TH3F* hResampledCharge=new TH3F("hResampledCharge","Resampled Charge using ChargeMapReader",nbins[0],low[0],high[0],nbins[1],low[1],high[1],nbins[2],low[2],high[2]);
  printf("Macro requesting r to fillChargeHistogram\n");
  r->FillChargeHistogram(hResampledCharge);
  printf("Returned to macro\n");

  
  //build a new reader so we can re-use the density map generation.  Binning doesn't matter here.
  printf("Macro building rSneaky\n");
  ChargeMapReader *rSneaky=new ChargeMapReader(nbins[1],low[1],high[1],nbins[0],low[0],high[0],nbins[2],low[2],high[2]);
  rSneaky->ReadSourceCharge(hResampledCharge);
  TH3* hResampledDensity=rSneaky->GetDensityHistogram();
  printf("Returned to macro\n");

  //build a third reader to make sure same-in makes same-out
  bool checkConsistency=false;
  TH3F* hCheckCharge;
  TH3* hCheckDensity;
  if (checkConsistency){
  printf("Macro building rCheck\n");
  hCheckCharge=new TH3F("hCheckCharge","Check Charge using ChargeMapReader",nbins[0],low[0],high[0],nbins[1],low[1],high[1],nbins[2],low[2],high[2]);
  rSneaky->FillChargeHistogram(hCheckCharge);
  ChargeMapReader *rCheck=new ChargeMapReader(nbins[1],low[1],high[1],nbins[0],low[0],high[0],nbins[2],low[2],high[2]);
  rCheck->ReadSourceCharge(hCheckCharge);
  hCheckDensity=rCheck->GetDensityHistogram();
  printf("Returned to macro\n");
  }

  

  TH1F* hFracChargeDiff=new TH1F("hFracChargeDiff","(Resampled-Original)/(Original) charge, should differ",100,-2,2);
  TH1F* hFracDensityDiff=new TH1F("hFracDensityDiff","(Resampled-Original)/(Original) density, should be the same",100,-2,2);
  TH1F* hFracChargeDiffCheck=new TH1F("hFracChargeDiffCheck","(Check-Original)/(Original) charge, should differ",100,-2,2);
  TH1F* hFracDensityDiffCheck=new TH1F("hFracDensityDiffCheck","(Check-Original)/(Original) density, should be the same",100,-2,2);
 
  // CompareHistograms(hOriginalDensity,hResampledDensity,hFracDensityDiff);
  CompareHistograms(hResampledDensity,hOriginalDensity,hFracDensityDiff);
  //  CompareHistograms(hOriginalCharge,hResampledCharge,hFracChargeDiff);
  CompareHistograms(hResampledCharge,hOriginalCharge,hFracChargeDiff);
  if (checkConsistency){
    CompareHistograms(hOriginalDensity,hCheckDensity,hFracDensityDiffCheck);
    CompareHistograms(hOriginalCharge,hResampledCharge,hFracChargeDiffCheck);
  }

 TAxis *resampledax[3]={nullptr,nullptr,nullptr};
  resampledax[0]=hResampledCharge->GetXaxis();
  resampledax[1]=hResampledCharge->GetYaxis();
  resampledax[2]=hResampledCharge->GetZaxis();

  float pos[]={2.2,500,500};//phi,r,z
  int sliceBin[2][3];
  for (int i=0;i<3;i++){
    sliceBin[0][i]=ax[i]->FindBin(pos[i]);
    sliceBin[1][i]=resampledax[i]->FindBin(pos[i]);
  }

  
  TCanvas *c=new TCanvas("c","c",1600,800);
  hOriginalCharge->SetName("hOriginalCharge");
  hOriginalCharge->SetTitle("hOriginalCharge");
  hOriginalDensity->SetName("hOriginalDensity");
  hOriginalDensity->SetTitle("hOriginalDensity");
   hResampledCharge->SetLineColor(kRed);
  hResampledDensity->SetLineColor(kRed);
  hResampledDensity->SetName("hResampledDensity");
  hResampledDensity->SetTitle("hResampledDensity");
  if (checkConsistency) {
    hCheckCharge->SetLineColor(kBlue);
    hCheckDensity->SetLineColor(kBlue);
    hCheckDensity->SetName("hCheckDensity");
    hCheckDensity->SetTitle("hCheckDensity");
  }

  
  c->Divide(4,2);
  int iPad=0;
  c->cd(++iPad);
  hOriginalCharge->ProjectionX("_phi",sliceBin[0][1],sliceBin[0][1],sliceBin[0][2],sliceBin[0][2])->Draw("hist");
  hResampledCharge->ProjectionX("_phir",sliceBin[1][1],sliceBin[1][1],sliceBin[1][2],sliceBin[1][2])->Draw("same,hist");
  if (checkConsistency) hCheckCharge->ProjectionX()->Draw("same,hist");
  c->cd(++iPad);
  hOriginalCharge->ProjectionY("_r",sliceBin[0][0],sliceBin[0][0],sliceBin[0][2],sliceBin[0][2])->Draw("hist");
  hResampledCharge->ProjectionY("_rr",sliceBin[1][0],sliceBin[1][0],sliceBin[1][2],sliceBin[1][2])->Draw("same,hist");
  if (checkConsistency) hCheckCharge->ProjectionY()->Draw("same,hist");
  c->cd(++iPad);
  hOriginalCharge->ProjectionZ("_z",sliceBin[0][0],sliceBin[0][0],sliceBin[0][1],sliceBin[0][1])->Draw("hist");
  hResampledCharge->ProjectionZ("_zr",sliceBin[1][0],sliceBin[1][0],sliceBin[1][1],sliceBin[1][1])->Draw("same,hist");
  if (checkConsistency) hCheckCharge->ProjectionZ()->Draw("same,hist");
  c->cd(++iPad);
  hFracChargeDiff->Draw();

  c->cd(++iPad);
  hOriginalDensity->ProjectionX("_phi",sliceBin[0][1],sliceBin[0][1],sliceBin[0][2],sliceBin[0][2])->Draw("hist");
  hResampledDensity->ProjectionX("_phir",sliceBin[1][1],sliceBin[1][1],sliceBin[1][2],sliceBin[1][2])->Draw("same,hist");
  if (checkConsistency) hCheckDensity->ProjectionX()->Draw("same,hist");
  c->cd(++iPad);
  hOriginalDensity->ProjectionY("_r",sliceBin[0][0],sliceBin[0][0],sliceBin[0][2],sliceBin[0][2])->Draw("hist");
  hResampledDensity->ProjectionY("_rr",sliceBin[1][0],sliceBin[1][0],sliceBin[1][2],sliceBin[1][2])->Draw("same,hist");
  if (checkConsistency) hCheckDensity->ProjectionY()->Draw("same,hist");
  c->cd(++iPad);
  hOriginalDensity->ProjectionZ("_z",sliceBin[0][0],sliceBin[0][0],sliceBin[0][1],sliceBin[0][1])->Draw("hist");
  hResampledDensity->ProjectionZ("_zr",sliceBin[1][0],sliceBin[1][0],sliceBin[1][1],sliceBin[1][1])->Draw("same,hist");
  if (checkConsistency) hCheckDensity->ProjectionZ()->Draw("same,hist");
  c->cd(++iPad);
  hFracDensityDiff->Draw();
 

  c->SaveAs("Test_ChargeMapReader.output.pdf");
  
  printf("All done.  Errors past here are root's problem.\n");
  return;
}
