

vector<std::pair<std::string,std::string>> fileName;
vector<std::string> histName;

const int resample_factor=11; //how many samples, evenly spaced but centered in the bin boundaries so that they're never on the edge unless n=inf, to resample.  may have aliasing issues if we're not careful.

void Resample(std::vector<TH3*> hin, std::vector<TH3*> hout);
void CheckClosure(std::vector<TH3*> hdistort, std::vector<TH3*> hcorrect, bool rFirst=false);//rFirst means apply rphi shift at the shifted r coord, instead of the original r coord
void ClosureTest(const char* originalfilename, const char* invertfilename, const char* closurefilename, bool rFirst=false); //opens appropriate files, then calls CheckClosure.
void invertHistograms(const char* originalfilename, const char* invertfilename, bool checkClosure=false);



void invertHistograms(){
  printf("proper usage:  root -b -q invertHistograms(const char* originalfilename, const char* invertfilename, bool checkClosure=false)\n");
  return;
}

void invertHistograms(const char* originalfilename, const char* invertfilename, bool checkClosure){
  //do the inversion and save the check histogram to the invert file:
  histName.push_back("hIntDistortionR_negz");
  histName.push_back("hIntDistortionP_negz");
  histName.push_back("hIntDistortionZ_negz");
  histName.push_back("hIntDistortionR_posz");
  histName.push_back("hIntDistortionP_posz");
  histName.push_back("hIntDistortionZ_posz");

  printf("Inverting:  %s --> %s\n",originalfilename,invertfilename);
  
  TFile *infile;
  TFile *outfile;
  std::vector<TH3*> hin;
  std::vector<TH3*> hout;
  
  infile=TFile::Open(originalfilename,"READ");
  outfile=TFile::Open(invertfilename,"RECREATE");
  outfile->cd();


  hin.clear();
  hout.clear();

  for (int i=0;i<3;i++){
    hin.push_back((TH3*)infile->Get(histName[i].data()));
    hout.push_back((TH3*)hin[i]->Clone(histName[i].data()));
    hout[i]->Reset();
  }
  Resample(hin,hout);
  for (int i=0;i<3;i++){
    hout[i]->Write();
  }
  if (checkClosure) CheckClosure(hin,hout);

  hin.clear();
  hout.clear();
    
  for (int i=0;i<3;i++){
    hin.push_back((TH3*)infile->Get(histName[i+3].data()));
    hout.push_back((TH3*)hin[i]->Clone(histName[i+3].data()));
    hout[i]->Reset();

  }
  Resample(hin,hout);
  for (int i=0;i<3;i++){
    hout[i]->Write();
  }
  

  
  return;
}


void ClosureTest(const char* originalfilename, const char* invertfilename, const char* closurefilename, bool rFirst){
  histName.push_back("hIntDistortionR_negz");
  histName.push_back("hIntDistortionP_negz");
  histName.push_back("hIntDistortionZ_negz");
  histName.push_back("hIntDistortionR_posz");
  histName.push_back("hIntDistortionP_posz");
  histName.push_back("hIntDistortionZ_posz");

 
  TFile *infile;
  TFile *invertfile;
  TFile *closurefile;

  std::vector<TH3*> hin;
  std::vector<TH3*> hout;

  printf("Checking Closure of Set:  Distort:%s Correct:%s ----> Residual:%s\n",originalfilename,invertfilename,closurefilename);
  infile=TFile::Open(originalfilename,"READ");
  invertfile=TFile::Open(invertfilename,"READ");
  closurefile=TFile::Open(closurefilename,"RECREATE");
  closurefile->cd();
    
  hin.clear();
  hout.clear();

  for (int i=0;i<3;i++){
    hin.push_back((TH3*)infile->Get(histName[i].data()));
    hout.push_back((TH3*)invertfile->Get(histName[i].data()));
  }
  for (int i=0;i<3;i++){
    hout[i]->Write();
  }
  CheckClosure(hin,hout);

  hin.clear();
  hout.clear();
    
  for (int i=0;i<3;i++){
    hin.push_back((TH3*)infile->Get(histName[i+3].data()));
    hout.push_back((TH3*)invertfile->Get(histName[i+3].data()));
  }
  for (int i=0;i<3;i++){
    hout[i]->Write();
  }

  CheckClosure(hin,hout, rFirst);
  //CheckClosure(hin,hin);

    
    infile->Close();
    invertfile->Close();
    closurefile->Close();
    return;
}



void Resample(std::vector<TH3*> hin, std::vector<TH3*> hout){
  TH3* hhits=(TH3*)hin[0]->Clone("hhits"); //number of elements in each output bin, for normalization purposes.
  TH1* hnhits=new TH1F("hnhits","number of hits in each output bin",resample_factor*resample_factor+1,-0.5,resample_factor*resample_factor+0.5);
  TH1* hdist=new TH1F("hdist","distortion sanity check",1000,-1,1);
  hhits->Reset();
  
  TAxis *ax[3]={nullptr,nullptr,nullptr};
  ax[0]=hhits->GetXaxis();
  ax[1]=hhits->GetYaxis();
  ax[2]=hhits->GetZaxis();

  int nbins[3];
  for (int i=0;i<3;i++){
    nbins[i]=ax[i]->GetNbins();//number of bins, not counting under and overflow.
  }


  unsigned long int ntotalsteps=nbins[0]*nbins[1]*nbins[2]*resample_factor*resample_factor*resample_factor;
  printf("Resampling.  %lu steps\n",ntotalsteps);
  int nWaypoints=20;
  printf("|");
  for (int i=0;i<nWaypoints-2;i++) printf("-");
  printf("|\n");
  unsigned long int ntotalstepWaypoint=ntotalsteps/nWaypoints;
  unsigned long int itotal=0;

  //remember that these histograms have an extra 'buffer' set of bins at their edges so that interpolation works correctly/
  //we do not wish to apply this procedure beyond the buffer bins, so need to ignore both the over/underflows and the bins adjacent
  //we will rebuild those edge cells separately.


  //   0     1     2   ...   n-1    n    n+1
  // under|guard|first|..|..|last|guard|over

  float low[3],high[3], step[3], sample_pos[3];
  float distortion[3], distorted_pos[3];
  int a=0;
  for (int i=2;i<nbins[0];i++){
    a=0;
    low[a]=ax[a]->GetBinLowEdge(i);
    high[a]=ax[a]->GetBinUpEdge(i);
    step[a]=(high[a]-low[a])/(1.*resample_factor);
    for (int isub=0;isub<resample_factor;isub++){
      a=0;
      sample_pos[a]=low[a]+step[a]*(isub+0.5);
      for (int j=2;j<nbins[1];j++){
	a=1;
	low[a]=ax[a]->GetBinLowEdge(j);
	high[a]=ax[a]->GetBinUpEdge(j);
	step[a]=(high[a]-low[a])/(1.*resample_factor);
	for (int jsub=0;jsub<resample_factor;jsub++){
	  a=1;
	  sample_pos[a]=low[a]+step[a]*(jsub+0.5);
	  for (int k=2;k<nbins[2];k++){
	    a=2;
	    low[a]=ax[a]->GetBinLowEdge(k);
	    high[a]=ax[a]->GetBinUpEdge(k);
	    step[a]=(high[a]-low[a])/(1.*resample_factor);
	    for (int ksub=0;ksub<resample_factor;ksub++){
	      itotal++;
	      unsigned long int progress=100*itotal/ntotalsteps;
		if (itotal%ntotalstepWaypoint==0){
		  printf("."); fflush(stdout);
		  //printf("Resampling... (itotal=%lu, waypoint=%lu) %lu%%\n",itotal,ntotalstepWaypoint,progress);
		}
	      a=2;
	      sample_pos[a]=low[a]+step[a]*(ksub+0.5);

	      //get the distorted position
	      for (int m=0;m<3;m++){
		distortion[m]=hin[m]->Interpolate(sample_pos[0],sample_pos[1],sample_pos[2]);

		if (m==0){
		  //because the map stores the phi-hat distortion, and we want the resulting coordinate:
		  distorted_pos[m]=sample_pos[m]+distortion[m]/sample_pos[1];
		  if (distorted_pos[m]>6.28) distorted_pos[m]-=6.28;
		  if (distorted_pos[m]<0) distorted_pos[m]+=6.28;
		} else {
		  distorted_pos[m]=sample_pos[m]+distortion[m];
		}
	
		hdist->Fill(distortion[m]);
	      }

	      //histogram the distortion in the distorted position.
	      hhits->Fill(distorted_pos[0],distorted_pos[1],distorted_pos[2],1.);
	      for (int m=0;m<3;m++){
		hout[m]->Fill(distorted_pos[0],distorted_pos[1],distorted_pos[2],distortion[m]);
	      }
	    }
	  }
	}
      }
    }
  }
  //we finished resampling all the useful bins.
  printf("\nCleaning up...\n"); //close off our progress dots.

  
  //   0     1     2   ...   n-1    n    n+1
  // under|guard|first|..|..|last|guard|over
   //normalize based on how many hits we put into each.  yell if we have less than 1/2 the sampling.
  float half_sampling=0.5*resample_factor*resample_factor*resample_factor;
  for (int i=2;i<nbins[0];i++){
    a=0;
    distorted_pos[a]=ax[a]->GetBinCenter(i);
    for (int j=2;j<nbins[1];j++){
      a=1;
      distorted_pos[a]=ax[a]->GetBinCenter(j);
      for (int k=2;k<nbins[2];k++){
	a=2;
	distorted_pos[a]=ax[a]->GetBinCenter(k);
	//histogram the distortion in the distorted position.
	int global_bin=hin[0]->FindBin(distorted_pos[0],distorted_pos[1],distorted_pos[2]);
	float global_hits=hhits->GetBinContent(global_bin);
	hnhits->Fill(global_hits);
	if (global_hits<half_sampling){
	  printf("(%2.2f,%2.2f,%2.2f)(glob=%d) has %1.2f entries\n",distorted_pos[0],distorted_pos[1],distorted_pos[2],global_bin,global_hits);
	  printf("this is less than half the expected entries, suggesting there are regions with poorly defined corrections.  Contact an expert.  Aborting.\n");
	  assert(false);
	} else {
	  //average the contents in the bin
	  for (int m=0;m<3;m++){
	    hout[m]->SetBinContent(global_bin,hout[m]->GetBinContent(global_bin)/global_hits);
	  }
	  hnhits->SetBinContent(global_bin,0);
	}     
      }
    }
  }


  
  //   0     1     2   ...   n-1    n    n+1
  // under|guard|first|..|..|last|guard|over
  //rebuild the edge bins:

  //first the sides:
  float target_pos[3], source_pos[3];
  for (int side_axis=0;side_axis<3;side_axis++){
    //the axis we're on the side of:
    for (int i=1;i<nbins[side_axis]+1;i+=nbins[side_axis]){
      a=side_axis;
      target_pos[a]=ax[a]->GetBinCenter(i);
      int source_i=i+1;
      if (i>1) source_i=i-1;
      source_pos[a]=ax[a]->GetBinCenter(source_i);

      //and the two axes we're in the middle of:
      for (int j=2;j<nbins[(1+side_axis)%3];j++){
	a=(1+side_axis)%3;
	target_pos[a]=ax[a]->GetBinCenter(j);
	source_pos[a]=ax[a]->GetBinCenter(j);
	for (int k=2;k<nbins[(2+side_axis)%3];k++){
	  a=(2+side_axis)%3;
	  source_pos[a]=ax[a]->GetBinCenter(k);
	  target_pos[a]=ax[a]->GetBinCenter(k);

	  //take the next bin 'inward' from our target bin and copy it.
	  int target_bin=hin[0]->FindBin(target_pos[0],target_pos[1],target_pos[2]);
	  int source_bin=hin[0]->FindBin(source_pos[0],source_pos[1],source_pos[2]);
	  for (int m=0;m<3;m++){
	    hout[m]->SetBinContent(target_bin,hout[m]->GetBinContent(source_bin));
	  }
	}
      }
    }
  }

 //now the edges:
  for (int side_axis=0;side_axis<3;side_axis++){
    //the axes we're on the side of:
    for (int i=1;i<nbins[side_axis]+1;i+=nbins[side_axis]){
      a=side_axis;
      target_pos[a]=ax[a]->GetBinCenter(i);
      int source_i=i+1;
      if (i>1) source_i=i-1;
      source_pos[a]=ax[a]->GetBinCenter(source_i);

      for (int j=1;j<nbins[(1+side_axis)%3]+1;j+=nbins[(1+side_axis)%3]){
	a=(1+side_axis)%3;
	target_pos[a]=ax[a]->GetBinCenter(j);
	int source_j=j+1;
	if (j>1) source_j=j-1;
	source_pos[a]=ax[a]->GetBinCenter(source_j);
	//and the axis we're in the middle of:

	for (int k=2;k<nbins[(2+side_axis)%3];k++){
	  a=(2+side_axis)%3;
	  source_pos[a]=ax[a]->GetBinCenter(k);
	  target_pos[a]=ax[a]->GetBinCenter(k);

	  //take the next bin 'inward' from our target bin and copy it.
	  int target_bin=hin[0]->FindBin(target_pos[0],target_pos[1],target_pos[2]);
	  int source_bin=hin[0]->FindBin(source_pos[0],source_pos[1],source_pos[2]);
	  for (int m=0;m<3;m++){
	    hout[m]->SetBinContent(target_bin,hout[m]->GetBinContent(source_bin));
	  }
	}
      }
    }
  }


 //now the corners:
  for (int side_axis=0;side_axis<3;side_axis++){
    //the axes we're on the side of:
    for (int i=1;i<nbins[side_axis]+1;i+=nbins[side_axis]){
      a=side_axis;
      target_pos[a]=ax[a]->GetBinCenter(i);
      int source_i=i+1;
      if (i>1) source_i=i-1;
      source_pos[a]=ax[a]->GetBinCenter(source_i);

      for (int j=1;j<nbins[(1+side_axis)%3]+1;j+=nbins[(1+side_axis)%3]){
	a=(1+side_axis)%3;
	target_pos[a]=ax[a]->GetBinCenter(j);
	int source_j=j+1;
	if (j>1) source_j=j-1;
	source_pos[a]=ax[a]->GetBinCenter(source_j);
	//and the axis we're in the middle of:

	for (int k=1;k<nbins[(2+side_axis)%3]+1;k+=nbins[(2+side_axis)%3]){
	  a=(2+side_axis)%3;
	  target_pos[a]=ax[a]->GetBinCenter(k);
	  int source_k=k+1;
	  if (k>1) source_k=k-1;
	  source_pos[a]=ax[a]->GetBinCenter(source_k);

	  //take the next bin 'inward' from our target bin and copy it.
	  int target_bin=hin[0]->FindBin(target_pos[0],target_pos[1],target_pos[2]);
	  int source_bin=hin[0]->FindBin(source_pos[0],source_pos[1],source_pos[2]);
	  for (int m=0;m<3;m++){
	    hout[m]->SetBinContent(target_bin,hout[m]->GetBinContent(source_bin));
	  }
	}
      }
    }
  }

  
  
  hdist->Write();
  hnhits->Write();
  return;
}



void CheckClosure(std::vector<TH3*> hdistort, std::vector<TH3*> hcorrect, bool rFirst){
  printf("Caution: Checking closure directly here may make a different (but internally consistent) assumption that the phi component is denominated in cm, and hence has a (minor) ordering issue.  Doing this through the actual Fun4All process would be better.\n");
  TH3* hclosure[3];
  TH1F* hresidual[3];
  TH2F* hresidual2D[3][3];
  static int ncalls=0;
  ncalls++;
  const char axischar[3]={'p','r','z'};
  for(int i=0;i<3;i++){
    hclosure[i]=(TH3*)hdistort[0]->Clone(Form("hclosure_%c_%d",axischar[i],ncalls));
    hclosure[i]->SetTitle(Form("%c-component of closure residual (f(g(x0))-x0 );phi;r;z;residual [cm]",axischar[i]));
    hclosure[i]->Reset();
    hresidual[i]=new TH1F(Form("hresidual%c_%d",axischar[i],ncalls),Form("%c-component of closure residual ( f(g(x0))-x0 ;residual [cm]",axischar[i]),200,-0.1,0.1);
  }
  TH3* hSuccesses=(TH3*)hdistort[0]->Clone(Form("hsuccesses_%d",ncalls));
  hSuccesses->Reset();
  TH3* hFailed=(TH3*)hdistort[0]->Clone(Form("hfailed_%d",ncalls));
 hFailed->Reset();


  TAxis *ax[3]={nullptr,nullptr,nullptr};
  ax[0]=hdistort[0]->GetXaxis();
  ax[1]=hdistort[0]->GetYaxis();
  ax[2]=hdistort[0]->GetZaxis();

  int nbins[3];
  for (int i=0;i<3;i++){
    nbins[i]=ax[i]->GetNbins();
  }

  //create our 2D plots which are the residual in direction A plotted as a function of position on axis B, and hence needs the binning of axis B
    double  lowEdges[1000];
  for (int i=0;i<3;i++){
    for (int j=0;j<3;j++){
      ax[j]->GetLowEdge(lowEdges);
      lowEdges[nbins[j]]=ax[j]->GetXmax();
      if (false){ //print debug information about the binning of the source file.
	printf("%c-axis reports %d bins, edges=",axischar[j],nbins[j]);
	for (int k=0;k<nbins[j]+1;k++){
	  printf("  %1.2f",lowEdges[k]);
	}
	printf("\n");
      }
      hresidual2D[i][j]=new TH2F(Form("hresidual%d_%d_%d",i,j,ncalls),Form("%c-component of closure residual;%c-position;%c-residual [cm]",axischar[i],axischar[j],axischar[i]),nbins[j],lowEdges,100,-0.1,0.1);
    }
  }



  //calculate the total number of bins, so we can give an accurate progress bar.
  int nWaypoints=40; //number of progress pips to display.
  unsigned long int ntotalsteps=nbins[0]*nbins[1]*nbins[2]*resample_factor*resample_factor*resample_factor;
  printf("Checking Closure.  %lu steps\n",ntotalsteps);
  printf("|");
  for (int i=0;i<nWaypoints-2;i++) printf("-");
  printf("|\n");
  unsigned long int ntotalstepWaypoint=ntotalsteps/nWaypoints;
  unsigned long int itotal=0;


  
  float low[3],high[3], step[3], sample_pos[3];
  float distortion[3], distorted_pos[3];
  float correction[3], corrected_pos[3];
  float residual[3];
  int a=0;


  
  //   0     1     2   ...   n-1    n    n+1
  // under|guard|first|..|..|last|guard|over
  for (int i=2;i<nbins[0];i++){
    a=0;
    low[a]=ax[a]->GetBinLowEdge(i);
    high[a]=ax[a]->GetBinUpEdge(i);
    step[a]=(high[a]-low[a])/(1.*resample_factor);
    for (int isub=0;isub<resample_factor;isub++){
      a=0;
      sample_pos[a]=low[a]+step[a]*(isub+0.5);
      for (int j=2;j<nbins[1];j++){
	a=1;
	low[a]=ax[a]->GetBinLowEdge(j);
	high[a]=ax[a]->GetBinUpEdge(j);
	step[a]=(high[a]-low[a])/(1.*resample_factor);
	for (int jsub=0;jsub<resample_factor;jsub++){
	  a=1;
	  sample_pos[a]=low[a]+step[a]*(jsub+0.5);
	  for (int k=2;k<nbins[2];k++){
	    a=2;
	    low[a]=ax[a]->GetBinLowEdge(k);
	    high[a]=ax[a]->GetBinUpEdge(k);
	    step[a]=(high[a]-low[a])/(1.*resample_factor);
	    for (int ksub=0;ksub<resample_factor;ksub++){

	      itotal++;
	      unsigned long int progress=100*itotal/ntotalsteps;
		if (itotal%ntotalstepWaypoint==0){
		  printf("."); fflush(stdout);
		}

	      
	      a=2;
	      sample_pos[a]=low[a]+step[a]*(ksub+0.5);

	      //check bounds and yell if we're weird:
	      for (int m=0;m<3;m++){
		int bin=ax[m]->FindBin(sample_pos[m]);
		if (bin<2 || bin>nbins[m]-1){
		  printf(" bin%d out of bounds:  input coord=%1.2f =%d/%d (bin edges (%1.2f to %1.2f))\n",
			 m,sample_pos[m],bin,nbins[m],ax[m]->GetBinLowEdge(bin),ax[m]->GetBinLowEdge(bin+1));
		  assert(false);
		}
	      }
	      //get the distorted position
	      for (int m=0;m<3;m++){//
		//bin center:  distortion[m]=hdistort[m]->GetBinContent(hdistort[m]->FindBin(sample_pos[0],sample_pos[1],sample_pos[2]));
		//interpolation version:
		distortion[m]=hdistort[m]->Interpolate(sample_pos[0],sample_pos[1],sample_pos[2]);
		if (m!=0) distorted_pos[m]=sample_pos[m]+distortion[m];
	      }
	      //because the map stores the phi-hat distortion, and we want the resulting coordinate, pos[0](phi) must be done separately:
	      if (rFirst){//distort r first, means use the distorted position to calculate the phi distorted position.
		distorted_pos[0]=sample_pos[0]+distortion[0]/distorted_pos[1];
	      } else {//distort r second, means use the sample position to calculate the phi distorted position.
		distorted_pos[0]=sample_pos[0]+distortion[0]/sample_pos[1];
	      }
	      //handle wrap-around in phi
	      if (distorted_pos[0]>6.28) distorted_pos[0]-=6.28;
	      if (distorted_pos[0]<0) distorted_pos[0]+=6.28;
		
	      //get the corrected position
	      //check bounds and yell if we're weird:
	      bool within_bounds=true;
	      for (int m=0;m<3;m++){
		int bin=ax[m]->FindBin(distorted_pos[m]);
		if (bin<2 || bin>nbins[m]-1){
		  //printf(" bin%d out of bounds:  distorted coord=%1.2f =%d/%d (bin edges (%1.2f to %1.2f))\n",
		  //	 m,sample_pos[m],bin,nbins[m],ax[m]->GetBinLowEdge(bin),ax[m]->GetBinLowEdge(bin+1));
		  within_bounds=false;
		  hFailed->Fill(sample_pos[0],sample_pos[1],sample_pos[2],1);

		}
	      }
	      if (within_bounds){
		//don't try to correct if the distortion should have crashed us.
		for (int m=0;m<3;m++){
		  //bin center: correction[m]=hcorrect[m]->GetBinContent(hcorrect[m]->FindBin(distorted_pos[0],distorted_pos[1],distorted_pos[2]));
		  correction[m]=hcorrect[m]->Interpolate(distorted_pos[0],distorted_pos[1],distorted_pos[2]);
		  if(m!=0)   corrected_pos[m]=distorted_pos[m]-correction[m];//MINUS! by convention
		}

	      //because the map stores the phi-hat correction, and we want the resulting coordinate, pos[0](phi) must be done separately:
	      if (rFirst){//correct r first, means use the corrected position to calculate the phi corrected position.
		corrected_pos[0]=distorted_pos[0]-correction[0]/corrected_pos[1];
	      } else {//correct r second, means use the distorted position to calculate the phi corrected position.
		corrected_pos[0]=distorted_pos[0]-correction[0]/distorted_pos[1];
	      }
	      //handle wrap-around in phi
	      if (corrected_pos[0]>6.28) corrected_pos[0]-=6.28;
	      if (corrected_pos[0]<0) corrected_pos[0]+=6.28;


		//fill the position residual
		for (int m=0;m<3;m++){
		
		  residual[m]=distortion[m]-correction[m];
		  hclosure[m]->Fill(sample_pos[0],sample_pos[1],sample_pos[2],residual[m]);
		  hresidual[m]->Fill(residual[m]);
		}
		hSuccesses->Fill(sample_pos[0],sample_pos[1],sample_pos[2],1);
		for (int f=0;f<3;f++){
		  for (int g=0;g<3;g++){
		    hresidual2D[f][g]->Fill(sample_pos[g],residual[f]);
		  }
		}
	      }//within bounds only

	      
	    }
	  }

	      
	}
      }
    }
  }

  //now we have to average the residual plots as well!

  //   0     1     2   ...   n-1    n    n+1
  // under|guard|first|..|..|last|guard|over
   //normalize based on how many hits we put into each.  yell if we have less than 1/2 the sampling.
  float half_sampling=0.5*resample_factor*resample_factor*resample_factor;
  for (int i=2;i<nbins[0];i++){
    a=0;
    distorted_pos[a]=ax[a]->GetBinCenter(i);
    for (int j=2;j<nbins[1];j++){
      a=1;
      distorted_pos[a]=ax[a]->GetBinCenter(j);
      for (int k=2;k<nbins[2];k++){
	a=2;
	distorted_pos[a]=ax[a]->GetBinCenter(k);
	//histogram the distortion in the distorted position.
	int global_bin=hSuccesses->FindBin(distorted_pos[0],distorted_pos[1],distorted_pos[2]);
	float global_hits=hSuccesses->GetBinContent(global_bin);
	if (global_hits<half_sampling){
	  printf("residual (%2.2f,%2.2f,%2.2f)(glob=%d) has %1.2f entries\n",distorted_pos[0],distorted_pos[1],distorted_pos[2],global_bin,global_hits);
	} else {
	  //average the contents in the bin
	  for (int m=0;m<3;m++){
	    hclosure[m]->SetBinContent(global_bin,hclosure[m]->GetBinContent(global_bin)/global_hits);
	  }
	}     
      }
    }
  }

  
  //TFile f=TFile::Open("ClosurePlots.hist.root");
  //note that we assume a file has been opened for our use!
  for (int i=0;i<3;i++){
    hclosure[i]->Write();
    hresidual[i]->Write();
  }
  //f->Close();
  return;
}

//int main(){
//  invertHistograms();
//  return 0;
//}
