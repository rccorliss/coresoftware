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
	hOutput->Fill((leadValue-followValue)/leadValue);
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

  return;
  //build a new reader so we can re-use the density map generation.  Binning doesn't matter here.
  ChargeMapReader *rSneaky=new ChargeMapReader(2,low[1],high[1],2,low[0],high[0],2,low[2],high[2]);
  TH3* hResampledDensity=rSneaky->GetDensityHistogram();

  TH1F* hFracChargeDiff=new TH1F("hFracChargeDiff","(Resampled-Original)/(Original) charge, should differ",100,-1,2);
  TH1F* hFracDensityDiff=new TH1F("hFracDensityDiff","(Resampled-Original)/(Original) density, should be the same",100,-1,2);
  
  CompareHistograms(hOriginalDensity,hResampledDensity,hFracDensityDiff);
  CompareHistograms(hOriginalCharge,hResampledCharge,hFracChargeDiff);

  TCanvas *c=new TCanvas("c","c",1000,1000);
  hResampledCharge->SetLineColor(kRed);
  hResampledDensity->SetLineColor(kRed);
  c->Divide(4,2);
  c->cd(1);
  hOriginalCharge->ProjectionX()->Draw("hist");
  hResampledCharge->ProjectionX()->Draw("same,hist");
  c->cd(2);
  hOriginalCharge->ProjectionY()->Draw("hist");
  hResampledCharge->ProjectionY()->Draw("same,hist");
  c->cd(3);
  hOriginalCharge->ProjectionZ()->Draw("hist");
  hResampledCharge->ProjectionZ()->Draw("same,hist");
  c->cd(4);
  hFracChargeDiff->Draw();
  
  c->cd(5);
  hOriginalDensity->ProjectionX()->Draw("hist");
  hResampledDensity->ProjectionX()->Draw("same,hist");
  c->cd(6);
  hOriginalDensity->ProjectionY()->Draw("hist");
  hResampledDensity->ProjectionY()->Draw("same,hist");
  c->cd(7);
  hOriginalDensity->ProjectionZ()->Draw("hist");
  hResampledDensity->ProjectionZ()->Draw("same,hist");
  c->cd(8);
  hFracDensityDiff->Draw();

  c->SaveAs("Test_ChargeMapReader.output.pdf");
  return;
}
