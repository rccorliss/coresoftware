//Loads a distortion map and generates a somewhat-realistic correction map following the basic procedure expected from each correction type.

enum inversionType{kPerfect=0,kStatic,kAverage,kFluctuation};
vector<std::string> histName;


//the guts:  Reads the distortion triplets at the truePos, generates the readoutPos, and assumes it came from the measPos.  Derives correction from this.
void ResampleOnGrid(std::vector<TH3*> hdistort, std::vector<TH3*> hcorrect, std::vector<TVector3>truePos, std::vector<TVector3>measPos);


//perfect correction:  generate the sampling points in a specified, fixed, fine grid
void mockup_perfect_correction(char* distortionfilename, char *correctionfilename,float rstep ,float phistep, float zstep);

//static correction:  generate the sampling points by throwing rays from the laser heads with uncertain direction expressed as a fraction of the step size, sample in z steps along these rays.
void mockup_static_correction(char* distortionfilename, char *correctionfilename,float phistep,float thetastep, float zstep, float phierr,float thetaerr, float zerr);


void mockup_correction_map(const char* invtype, char* distortionfilename, char* correctionfilename){
  std::string typestring[]={"perfect","static","average","fluctuation"};
  int type=-1;
  for (int i=0;i<4;i++){
    if (typestring[i].compare(invtype)==0){//yup, here '0' means 'they match'
      type=i;
      break;
    }
  }

  //declare the histogram names we'll be inverting
  histName.push_back("hIntDistortionR_negz");
  histName.push_back("hIntDistortionP_negz");
  histName.push_back("hIntDistortionZ_negz");
  // histName.push_back("hIntDistortionR_posz");
  // histName.push_back("hIntDistortionP_posz");
  // histName.push_back("hIntDistortionZ_posz");

  switch(type){
  case kPerfect:
    mockup_perfect_correction(distortionfilename,correctionfilename,1,0.2,1);
    return;
  case kStatic:
    mockup_static_correction(distortionfilename,correctionfilename,2*6.28/4096.0,2*2*6.28/4096.0,2,1,1,0.5);
    return;
  case kAverage:
    printf("not implemented yet\n");
    return;
  case kFluctuation:
    printf("not implemented yet\n");
    return;
  default:
    printf("invalid correction type.  arg must be:");
    for (int i=0;i<4;i++){
      printf("%s ",typestring[i].c_str());
    }
    printf("\n");
    return;
  }
  

  return;
}

void mockup_perfect_correction(char* distortionfilename, char *correctionfilename,float rstep, float phistep, float zstep){
  //take steps computed as a rotation in phi and theta, with fractional uncertainty as specified.
  TFile *infile;
  TFile *outfile;

  std::vector<TH3*> hin;
  std::vector<TH3*> hout;

  infile=TFile::Open(distortionfilename,"READ");
  outfile=TFile::Open(correctionfilename,"RECREATE");
  outfile->cd();
    
  hin.clear();
  hout.clear();


  std::vector<TVector3> truePos;
  std::vector<TVector3> measPos;
  for (float r=20;r<78;r+=rstep){
    for (float phi=0;phi<6.28;phi+=phistep){
      for (float z=-105.5;z>0;z+=zstep){
	TVector3 pos(1,1,z);
	pos.SetPerp(r);
	pos.SetPhi(phi);
	truePos.push_back(pos);
	measPos.push_back(pos);
	//and mirror them to do the same on the positive side
	// truePos.push_back(trueLaser*-1.);
	// measPos.push_back(measLaser*-1.);
      }
    }
  }

	
	
  

  //clone each histogram to copy its dimensions, then clear the data out:
  
  for (int i=0;i<histName.size();i++){
    hin.push_back((TH3*)infile->Get(histName[i].data()));
    hout.push_back((TH3*)hin[i]->Clone(histName[i].data()));
    hout[i]->Reset();
  }
  ResampleOnGrid(hin,hout,truePos,measPos);

  //gotta put some code in here to save the result.

  for (int i=0;i<hout.size();i++){
    hout[i]->Write();
  }
  outfile->Close();
  return;

}


void mockup_static_correction(char* distortionfilename, char *correctionfilename,float phistep,float thetastep, float zstep, float phierr,float thetaerr, float zerr){
  //take steps computed as a rotation in phi and theta, with fractional uncertainty as specified.
  TFile *infile;
  TFile *outfile;

  std::vector<TH3*> hin;
  std::vector<TH3*> hout;

  infile=TFile::Open(distortionfilename,"READ");
  outfile=TFile::Open(correctionfilename,"RECREATE");
  outfile->cd();
    
  hin.clear();
  hout.clear();

  TRandom rand;
  printf("generating static correction test points with phistep=%f,thetastep=%f,zstep=%f\n",phistep,thetastep,zstep);
  
  //create the grid of pointswe wish to sample at, and the position we think they are at.
  int nLasers=4;
  TVector3 laserHead[nLasers];
  for (int i=0;i<nLasers;i++){
    laserHead[i].SetXYZ(60,0,-105.5);
    laserHead[i].RotateZ(3.14/2.0*i+15*3.14/360.);
  }
  std::vector<TVector3> truePos;
  std::vector<TVector3> measPos;
  for (float theta=0;theta<3.00/2.0;theta+=thetastep){
    //printf("theta=%f\n",theta);
    for (float phi=0;phi<6.28;phi+=phistep){
      TVector3 laserBeam(0,0,1);
      laserBeam.RotateX(theta);
      laserBeam.RotateZ(phi);
      TVector3 measBeam(0,0,1);
      measBeam.RotateX(theta+thetastep*rand.Gaus());
      measBeam.RotateZ(phi+phistep*rand.Gaus());
      laserBeam=zstep/laserBeam.Z()*laserBeam;//scale laser beam so that its z component is one z step.
      //printf("laserBeam z size=%f\n",laserBeam.Z());
      for (int i=0;i<nLasers;i++){
	TVector3 trueLaser=laserHead[i];
	TVector3 measLaser=laserHead[i];
	for (float z=-105.5;z<0;z+=zstep){
	  trueLaser=trueLaser+laserBeam;
	  measLaser=measLaser+measBeam;
	  if (trueLaser.Perp()>20 && trueLaser.Perp()<80 && trueLaser.Z()<0 && trueLaser.Z()>=-105.6){
	    truePos.push_back(trueLaser);
	    measPos.push_back(measLaser);
	    //and mirror them to do the same on the positive side
	    // truePos.push_back(trueLaser*-1.);
	    // measPos.push_back(measLaser*-1.);
	  } else {
	    //we hit the IFC, OFC, or CM
	    printf("laser reached r=%1.3f,z=%1.3f. Stopping propagation\n",trueLaser.Perp(),trueLaser.Z());
	    break;
	  }
	}
      }
    }
  }
	
	
  

  //clone each histogram to copy its dimensions, then clear the data out:
  
  for (int i=0;i<histName.size();i++){
    hin.push_back((TH3*)infile->Get(histName[i].data()));
    hout.push_back((TH3*)hin[i]->Clone(histName[i].data()));

   hout[i]->Reset();
  }
  printf("Issuing ResampleOnGrid with %d true, %d meas\n",truePos.size(),measPos.size());
  ResampleOnGrid(hin,hout,truePos,measPos);

  //gotta put some code in here to save the result.

  for (int i=0;i<hout.size();i++){
    hout[i]->Write();
  }
  outfile->Close();
  return;

}


void ResampleOnGrid(std::vector<TH3*> hdistort, std::vector<TH3*> hcorrect, std::vector<TVector3>truePos, std::vector<TVector3>measPos){
  //take in a triplet of distortions (phi,r,z) and compute the corrections on the grid points given.
  
  TH3* hhits=(TH3*)hdistort[0]->Clone("hhits"); //number of elements in each output bin, for normalization purposes.
  hhits->Reset();
  
  TAxis *ax[3]={nullptr,nullptr,nullptr};
  ax[0]=hhits->GetXaxis();
  ax[1]=hhits->GetYaxis();
  ax[2]=hhits->GetZaxis();

  int nbins[3];
  float axismin[3], axismax[3];
  for (int i=0;i<3;i++){
    nbins[i]=ax[i]->GetNbins();//number of bins, not counting under and overflow.
    axismin[i]=ax[i]->GetBinLowEdge(2);//since 0 is underflow and 1 is the guard bin
    axismax[i]=ax[i]->GetBinLowEdge(nbins[i]);//since n is the guard bin, its low edge is the top of the last real bin
  }

  //generate coarser and coarser maps
  int ntempbins[3];
  float tempbinwidth[3];
  std::vector<TH3F*> hCoarse[3], hCoarseHits;
  for (int i=0;i<3;i++){
    ntempbins[i]=nbins[i];
  }
  //make as many coarser versions as we can:
  for (int i=0, okay=true;okay;i++){
    //find the next coarser size
    for (int j=0;j<3;j++){
      ntempbins[j]=ntempbins[j]/2;
      tempbinwidth[j]=(axismax[j]-axismin[j])/(1.0*ntempbins[j]);
      if (ntempbins[j]<2) okay=false; //stop when we don't have any meaningful bins in one dimension.
    }
    if(okay){
      printf("Assembling coarseness=%d, bins=%d,%d,%d\n",i,ntempbins[0],ntempbins[1],ntempbins[2]);
    //define the histogram with that size:
    for (int j=0;j<3;j++){//one for each distortion direction:
      hCoarse[j].push_back(new TH3F(Form("hCoarseTemp%d_%d",i,j),Form("hCoarseTemp%d_%d",i,j),
			       ntempbins[0]+2,axismin[0]-tempbinwidth[0],axismax[0]+tempbinwidth[0],
			       ntempbins[1]+2,axismin[1]-tempbinwidth[1],axismax[1]+tempbinwidth[1],
			       ntempbins[2]+2,axismin[2]-tempbinwidth[2],axismax[2]+tempbinwidth[2]));
    }
    hCoarseHits.push_back(new TH3F(Form("hCoarseTempHits%d",i),Form("hCoarseTempHits%d",i),
				   ntempbins[0]+2,axismin[0]-tempbinwidth[0],axismax[0]+tempbinwidth[0],
				   ntempbins[1]+2,axismin[1]-tempbinwidth[1],axismax[1]+tempbinwidth[1],
				   ntempbins[2]+2,axismin[2]-tempbinwidth[2],axismax[2]+tempbinwidth[2]));
    }
  }
  
  printf("Resampling.  %d steps\n",truePos.size());

  float distortion[3];
  TVector3 distortedPos;
  for (int i=0;i<truePos.size();i++){
    float phitrue=truePos[i].Phi();
    if (phitrue<0) phitrue+=6.28;
    float phimeas=measPos[i].Phi();
    if (phimeas<0) phimeas+=6.28;

    //get the real distortion at the true position:
    for (int m=0;m<3;m++){
      distortion[m]=hdistort[m]->Interpolate(phitrue,truePos[i].Perp(),truePos[i].Z());
    }
    //but figure out what the distorted position is:
    distortedPos=truePos[i];
    distortedPos.SetPhi(truePos[i].Phi()+distortion[0]/truePos[i].Perp());
    distortedPos.SetZ(truePos[i].Z()+distortion[2]);
    distortedPos.SetPerp(truePos[i].Perp()+distortion[1]);
    float distortedPhi=distortedPos.Phi();
    if (distortedPhi<0) distortedPhi+=6.28;
    
    //and compute the inferred distortion with the assumption that it came from the measured position instead:
    float measDistortion[3];
    measDistortion[0]=distortedPos.DeltaPhi(measPos[i])*measPos[i].Perp();
    measDistortion[1]=distortedPos.Perp()-measPos[i].Perp();
    measDistortion[2]=distortedPos.Z()-measPos[i].Z();

    //and fill that inferred distortion at the distorted position.
    hhits->Fill(distortedPhi,distortedPos.Perp(),distortedPos.Z());
    for (int c=0;c<hCoarseHits.size();c++){
      hCoarseHits.at(c)->Fill(distortedPhi,distortedPos.Perp(),distortedPos.Z());
    }
    for (int m=0;m<3;m++){
      hcorrect[m]->Fill(distortedPhi,distortedPos.Perp(),distortedPos.Z(),measDistortion[m]);
      for (int c=0;c<hCoarse[m].size();c++){
	hCoarse[m].at(c)->Fill(distortedPhi,distortedPos.Perp(),distortedPos.Z(),measDistortion[m]);
      }
    }
  }


  //now normalize our distortion correction and clean up the histogram:
  //remember that these histograms have an extra 'buffer' set of bins at their edges so that interpolation works correctly/
  //we do not wish to apply this procedure beyond the buffer bins, so need to ignore both the over/underflows and the bins adjacent
  //we will rebuild those edge cells separately.

  
  //   0     1     2   ...   n-1    n    n+1
  // under|guard|first|..|..|last|guard|over
  int a;
  float distorted_pos[3];
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
	int global_bin=hcorrect[0]->FindBin(distorted_pos[0],distorted_pos[1],distorted_pos[2]);
	float global_hits=hhits->GetBinContent(global_bin);
	//hnhits->Fill(global_hits);
	if (global_hits<1){
	  printf("(%2.2f,%2.2f,%2.2f)(glob=%d) has %1.1f entries. ",distorted_pos[0],distorted_pos[1],distorted_pos[2],global_bin,global_hits);
	  int c=-1;
	  int coarse_bin=global_bin;
	  //this is not the best way to do this, but is a decent start.
	  while (global_hits==0){//continue looking into coarser and coarser bins until we find one with nonzero contents
	    c++;
	    coarse_bin=hCoarseHits.at(c)->FindBin(distorted_pos[0],distorted_pos[1],distorted_pos[2]);
	    global_hits=hCoarseHits.at(c)->GetBinContent(coarse_bin);
	  }
	  printf(" ==> fills at coarseness=%d, coarse bin has %1.1f hits\n",c,global_hits);
	  //fill hcorrect with the first nonzero contents, which is an average of the first region that has at least one hit in it.  This is not a great guess, but is better than a zero.
	  //todo:  we could fit the point cloud instead, and sample our fit functions (splines?) at the missing points.
	  for (int m=0;m<3;m++){
	    hcorrect[m]->SetBinContent(global_bin,hCoarse[m].at(c)->GetBinContent(coarse_bin)/global_hits);
	  }					      
	  //exit;
	} else {
	  //average the contents in the bin
	  for (int m=0;m<3;m++){
	    hcorrect[m]->SetBinContent(global_bin,hcorrect[m]->GetBinContent(global_bin)/global_hits);
	  }
	  //hnhits->SetBinContent(global_bin,0);
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
	  int target_bin=hcorrect[0]->FindBin(target_pos[0],target_pos[1],target_pos[2]);
	  int source_bin=hcorrect[0]->FindBin(source_pos[0],source_pos[1],source_pos[2]);
	  for (int m=0;m<3;m++){
	    hcorrect[m]->SetBinContent(target_bin,hcorrect[m]->GetBinContent(source_bin));
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
	  int target_bin=hcorrect[0]->FindBin(target_pos[0],target_pos[1],target_pos[2]);
	  int source_bin=hcorrect[0]->FindBin(source_pos[0],source_pos[1],source_pos[2]);
	  for (int m=0;m<3;m++){
	    hcorrect[m]->SetBinContent(target_bin,hcorrect[m]->GetBinContent(source_bin));
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
	  int target_bin=hcorrect[0]->FindBin(target_pos[0],target_pos[1],target_pos[2]);
	  int source_bin=hcorrect[0]->FindBin(source_pos[0],source_pos[1],source_pos[2]);
	  for (int m=0;m<3;m++){
	    hcorrect[m]->SetBinContent(target_bin,hcorrect[m]->GetBinContent(source_bin));
	  }
	}
      }
    }
  }
  return;
}


    
