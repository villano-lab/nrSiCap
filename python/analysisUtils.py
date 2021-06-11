

##Set of utilities useful for anaylsing rq data in root
##N. Mast, 2019
##
#include <vector> #Probably don't need to replace this, likely built into python.
#If I do it'll probably be numpy or something easy.

import uproot

#######################
##Load RQ's into a TChain
#######################
def loadSeriesRQ(seriesNumber, test="", path='data', zip="1"):
  filepath = path+test+"/byseries/"+seriesNumber+"/umn_"+seriesNumber
  Dump1fname = filepath+"_F0001.root"
  ##Check if file exists
  try:
        uproot.open(Dump1fname)
  except: 
    print("No rq files found for "+seriesNumber)
    print(filename+" does not exist.\n")
  zipName = "rqDir/zip"+zip
  
    # THIS IS WHERE I STOPPED - HOW TO TRANSLATE TCHAIN??
    
    
  TChain *z = new TChain(zipName);
  TChain *e = new TChain("rqDir/eventTree");
  z->AddFriend(e);
  filename.Append("*.root");
  z->Add(filename);
  e->Add(filename);
  return z

#######################
##Load RQ's into a TChain
#######################
TChain* loadSeriesRQFull(TString filepath ,TString zip="1"){
  TString zipName = "rqDir/zip";
  zipName.Append(zip);
  TChain *z = new TChain(zipName);
  TChain *e = new TChain("rqDir/eventTree");
  z->AddFriend(e);
  filepath.Append("*.root");
  z->Add(filepath);
  e->Add(filepath);
  return z;
}

#######################
##Save TCut to file
#######################
void saveCuts(TCut cuts, TString fileName){
  ofstream file(fileName);
  if (file.is_open()){
    cout << "Saving cuts: " << cuts.GetTitle() <<'\n';
    file << cuts.GetTitle() << '\n';
    file.close();
  }else{
    cout<<"Unable to open file.\n";
  }
}

#######################
##Load TCut from File
#######################
TCut loadCuts(TString fileName){
  string line;
  TString cutString;
  ifstream file(fileName);
  if (file.is_open()){
    getline(file,line);
    cout << line <<'\n';
    cutString=line;
    file.close();
    return cutString.Data();
  }else{
    cout<<"Unable to open file.\n";
    return "";
  }
}
#######################
##Save eventlist to file
#######################
Int_t eventListToFile(TChain *z, TCut cuts,TString fileName, Int_t maxEvents){
  cout << "Writing eventList to file.\n";

  z->Draw(">>elist",cuts);
  TEventList *elist = (TEventList*)gDirectory->Get("elist");

  Long64_t entryIndex=0;
  Double_t SeriesNumberCurr,SeriesNumberNew, EventNumber;
  Double_t prevEventNumber=0;
  char rawFileNameCurr[100];
  char rawFileNameNew[100];
  ofstream file;
  file.open(fileName);
  Int_t eventCount=0;
  Int_t i=0;

  entryIndex=elist->GetEntry(i);
  while (entryIndex!=-1&&eventCount<maxEvents){
    z->GetEntry(entryIndex);
    SeriesNumberNew=z->GetLeaf("SeriesNumber")->GetValue();
    strcpy(rawFileNameNew,z->GetFile()->GetName());

    ##Single File for all events, separated by series,dumps
    if (strcmp(rawFileNameNew,rawFileNameCurr)){
      file << "**********\n";
      file << rawFileNameNew << '\n';
      ##file << Form("%.0f\n",SeriesNumberNew);
      strcpy(rawFileNameCurr,rawFileNameNew);
      }
    EventNumber=z->GetLeaf("EventNumber")->GetValue();
    ##Only take 1 of every 10 adjascent events to avoid getting retriggers
    if(1||fabs(EventNumber-prevEventNumber)>10){
      prevEventNumber=EventNumber;
      file << Form("%.0f\n",EventNumber);
      eventCount++;
    }

    i++;
    entryIndex=elist->GetEntry(i);
    }
  file.close();
  cout << "Wrote " << eventCount << " events to: " << fileName <<'\n';
  return eventCount;
}

#######################
##Save RQs to json file
#######################
Int_t RQsToJSON(TChain *z, TCut cuts,TString fileName, vector<TString> RQList, Int_t maxEvents){
  ##This is a cumbersome way to do this, but I can't figure out a better one since you can't
  ##  access values of calculated aliases as you would a normal branch

  cout << "Collecting RQs from Tree.\n";
  vector<vector<double>> RQs;
  ##Generate arrays of desired RQs ordered by Tree Entry$ number
  for(std::vector<TString>::iterator it = RQList.begin(); it != RQList.end(); ++it) {
    int Nentries=z->Draw(Form("%s%s",(*it).Data(),":Entry$"),cuts,"goff");
    if(Nentries==-1){
      ##Got an error drawing
      cout<<"Error, bad RQ."<<endl;
      return -1;
    }
    double* RQarr=z->GetV1();##RQ

    vector<double> RQtemp;
    for(int i=0;i<Nentries;i++){
      RQtemp.push_back(RQarr[i]);
    }
    RQs.push_back(RQtemp);
  }

  ##Loop through arrays, building and outputting a json-style entry for each Tree entry
  cout << "Writing RQs to JSON.\n";

  ofstream fout(fileName);
  Int_t eventCount=0;

  ##Begin overall library
  fout << "{"<<endl;
  ##Begin main pulse data list
  fout << "\"data\":[" << endl;

  for(int iEvent=0;iEvent<RQs[0].size();iEvent++){
    ##Write the current event     {"RQ1":071705231956,"RQ2":10001},
    if(iEvent>0){fout << "," << endl;}##comma after the previous event's block
    fout<<"\t{";##Start this event's block
    for(int iRQ=0; iRQ < RQList.size(); iRQ++){
      if(iRQ>0){fout<<",";}##Comma after previous RQ

      if(RQList[iRQ]=="SeriesNumber"){fout << "\""<<RQList[iRQ].Data()<<"\":" << Form("%.0f",RQs[iRQ][iEvent]);}
      else{fout << "\""<<RQList[iRQ].Data()<<"\":" << RQs[iRQ][iEvent];}
    }
    fout<<"}";##End this event's block
    eventCount++;
    if(eventCount>=maxEvents){break;}
  }
  ##Close data list
  fout << endl << "]" << endl;
  ##Close overall library
  fout << "}" << endl;

  cout<<"Closing file"<<endl;
  ##close the raw data file
  fout.close();

  cout << "Wrote " << eventCount << " events to: " << fileName <<'\n';
  return eventCount;
}



#################################################
##Begin PulseTools-like functions
#################################################

##Remove a trace's baseline
void BSsubtract(TH1D *h, Double_t bs){
  for (Int_t i=1; i<=h->GetNbinsX(); i++){
    h->AddBinContent(i,-bs);
  }
}

##Remove a trace's baseline
##Define pre and post numbers of bins to use in calculateing baseline
##Removes linear background based on beginning and ending of trace
void BSsubtractSlope(TH1D *h, int iStart, int nPre, int nPost){
  int nBins = h->GetNbinsX();
  
  if(nPost>0){##Do the sloped BS subtraction
	  double pre = getBSpre(h, iStart, nPre);
	  double post = getBSpost(h, nPost);
	  double m = (post-pre)/(nBins-nPre/2-nPost/2);
	  double b = ((nBins-nPost/2)*pre-(nPre/2)*post) / (nBins-nPre/2-nPost/2);

	  double val;

	  for (Int_t i=1; i<=nBins; i++){
	    val=h->GetBinContent(i);
	    val=val-(m*i+b);
	    h->SetBinContent(i,val);
	  }
  }else{##Just do the flat, pre pulse BS subtraction
	BSsubtract(h, getBSpre(h,iStart, nPre));
  }
}

##Rescale a trace's x (time) axis
##h is the existing trace, hNew will contain the rescaled version
##Bin centers will start with t0 and be separated by dt.
TH1D* rescaleX(TH1D *h, double t0, double dt){
  int nBins = h->GetNbinsX();
  TH1D *hNew = new TH1D(h->GetName(),h->GetTitle(),nBins,t0-dt/2,t0+dt*(nBins-0.5));
  for(int i=1; i<=nBins; i++){
    hNew->SetBinContent(i,h->GetBinContent(i));
  }
  return hNew;
}

##Rescale a trace's x axis
##h is the existing trace, hNew will contain the rescaled version
##scale x bins by factor f
##Assumes constant bin widths
TH1D* rescaleX(TH1D *h, double f){
  int nBins = h->GetNbinsX();
  double low = f*h->GetBinLowEdge(1);
  double hi = f*h->GetBinLowEdge(nBins+1);
  TH1D *hNew = new TH1D(h->GetName(),h->GetTitle(),nBins,low,hi);
  for(int i=1; i<=nBins; i++){
    hNew->SetBinContent(i,h->GetBinContent(i));
  }
  return hNew;
}


##calculate a trace's beginning baseline
/*double getBSpre(TH1D *h, int Nbins){
  double sum=0;
  for (Int_t i=1; i<=Nbins; i++){
    sum=sum+h->GetBinContent(i);
  }
  return sum/Nbins;
}
*/

double getBSpre(TH1D *h, int binLow, int binHi){
  double sum=0;
  for (Int_t i=binLow; i<=binHi; i++){
    sum=sum+h->GetBinContent(i);
  }
  return sum/(binHi-binLow+1);
}

##calculate a trace's ending baseline
double getBSpost(TH1D *h, int Nbins){
  double sum=0;
  for (Int_t i=h->GetNbinsX()-Nbins+1; i<=h->GetNbinsX(); i++){
    sum=sum+h->GetBinContent(i);
  }
  return sum/Nbins;
}

##Crude low-pass filter
##tau, tStart, tEnd have same units as x-axis of h
##Use average of back and forth filtering
void lpf(TH1D *h, double tau, int tStart, int tEnd){

  double dt=h->GetXaxis()->GetBinCenter(2)-h->GetXaxis()->GetBinCenter(1);
  double alpha = dt/tau;
  TH1D *htemp = (TH1D*)h->Clone();

  int iStart=h->FindBin(tStart);
  int iEnd=h->FindBin(tEnd);

  if(iStart<=0){iStart=1;}
  if(iEnd<=0 || iEnd>h->GetNbinsX()){iEnd=h->GetNbinsX();}

  ##Forward pass
  double yk1=h->GetBinContent(iStart);
  for(int i=iStart; i<=iEnd; i++){
    yk1 = yk1 + alpha*(h->GetBinContent(i)-yk1);
    htemp->SetBinContent(i,yk1);
  }
  ##Backwards pass
  double yk2=h->GetBinContent(h->GetNbinsX());
  for(int i=iEnd; i>=iStart; i--){
    yk2 = yk2 + alpha*(h->GetBinContent(i)-yk2);
    yk1=htemp->GetBinContent(i);
    h->SetBinContent(i,(yk1+yk2)/2);
  }

  return;
}

##using darkpipe symmetric convention for normalization
void RealToComplexFFT(const vector<double>& pulsevector, vector<TComplex>& outComp)
{
    int n = pulsevector.size();
    if(n == 0)
    {
      cerr <<"RealToComplexFFT - ERROR! empty pulse passed into this function" << endl;
      exit(1);
    }

    double *pvector = new double[pulsevector.size()];
    double re,im;
    int i;

   ##copying std vector back into array to pass to FFT (to mirror what is done in PulseTools)
   for (i=0;i<(int)pulsevector.size();i++){
      pvector[i] = pulsevector[i];
   }

   TVirtualFFT *fftr2c = TVirtualFFT::FFT(1,&n,"R2C ES");
   fftr2c->SetPoints(pvector);
   fftr2c->Transform();
   fftr2c->GetPointComplex(0,re,im);
   TComplex tempCopy(re/sqrt((double)n), im/sqrt((double)n));
   outComp.push_back(tempCopy);

   for(i=1; i<n ;i++)
   {
      fftr2c->GetPointComplex(i,re,im);
      TComplex tempCopy(re/sqrt((double)n), im/sqrt((double)n));
      outComp.push_back(tempCopy);
   }

   ##done! so delete new'd objects
   delete[] pvector;
   delete fftr2c;

   return;

}


void ComplexToRealIFFT(const vector<TComplex>& inComp, vector<double>& outRe)
{
    int n = inComp.size();
    if(n == 0)
    {
      cerr <<"PulseTools::ComplexToRealIFFT - ERROR! empty pulse passed into this function" << endl;
      exit(1);
    }

    double *re_pvector = new double[inComp.size()];
    double *im_pvector = new double[inComp.size()];

    double re,im;
    int i;

    ##copying std vector back into arrays to pass to IFFT
    for (i=0; i < n; i++){
       re_pvector[i] = inComp[i].Re();
       im_pvector[i] = inComp[i].Im();
    }

    TVirtualFFT *ifftc2r = TVirtualFFT::FFT(1,&n,"C2R ES"); ##this should be the opposite of FFT
    ifftc2r->SetPointsComplex(re_pvector, im_pvector);
    ifftc2r->Transform();
    ifftc2r->GetPointComplex(0,re,im);
    outRe.push_back(re/sqrt((double)n)); ##DC component

    for(i=1; i<n ;i++)
    {
       ifftc2r->GetPointComplex(i,re, im);
       outRe.push_back(re/sqrt((double)n));
    }

   delete[] re_pvector;
   delete[] im_pvector;
   delete ifftc2r;

   return;
}


#################################################
##End PulseTools
#################################################
void DrawTraces(TH1D **hArr, TString title, Int_t nTraces=9, Double_t xmin=0, Double_t xmax=3.27640000000000000e+06){
  if (nTraces>9){
    cout<<"Not ready for more than 9 traces"<<endl;
    return;
  }
  Int_t goodTraces[]={0,1,2,3,4,5,6,7,10};

  EColor color[]={kCyan,kCyan+2,kBlue,kBlue+2,kMagenta,kMagenta+2,kRed,kRed+2,kBlack};
  ##Set the colors and find the largest amplitude
  Double_t max=0;
  for(Int_t i=0;i<nTraces;i++){
    Int_t index=goodTraces[i];
    if (max<hArr[index]->GetMaximum()){max=hArr[index]->GetMaximum();}
    hArr[index]->SetLineColor(color[i]);
  }
  ##Draw it
  hArr[0]->GetYaxis()->SetRangeUser(-100,100*int(ceil(max/100.0)));
  hArr[0]->GetXaxis()->SetRangeUser(xmin,xmax); 
  hArr[0]->GetXaxis()->SetTitle("time (ns)");
  hArr[0]->SetTitle(title);
  hArr[0]->Draw();
  for(Int_t i=1;i<nTraces;i++){
    Int_t index=goodTraces[i];
    hArr[index]->Draw("same");
  }
  
  TLegend *leg = new TLegend(0.7, 0.5, 0.9, 0.9);
  for(Int_t i=0;i<9;i++){
    Int_t index=goodTraces[i];
    leg->AddEntry(rawPulse[index],rawPulse[index]->GetName());
  }
  leg->Draw();
  
  return;
}



##Shift a trace by some number of bins, appending 0 signal where needed
##Only use this after removing the baseline from a trace
void shiftTrace(TH1D *h, Int_t shift){
	if (shift>0){##Moving the pulse later in time.
	  ##Shift contents towards end of trace
	  for(Int_t i=h->GetNbinsX();i>shift;i--){
	    h->SetBinContent(i,h->GetBinContent(i-shift));
	  }
	  ##Zero the beginning of the trace
	  for(Int_t i =shift;i>0;i--){
	    h->SetBinContent(i,0);
	  }
	}else{##Move pulse earlier in time
	  ##Shift contents toward beginning of trace
	  Double_t value;
	  for(Int_t i=1;i<h->GetNbinsX()-abs(shift);i++){
	    value=h->GetBinContent(i+abs(shift));
	    h->SetBinContent(i,value);
	  }
	  ##Zero the end of the trace
	  for(Int_t i=h->GetNbinsX()-abs(shift);i<=h->GetNbinsX();i++){
	    h->SetBinContent(i,0);
	  }
	}
	return;
}

##Find the 50%RT bin
##Use only on bs subtracted traces
int getRTbin(TH1D* h, Double_t risePercent=0.5, Int_t firstBin=-1, Int_t lastBin=-1){
	if (firstBin==-1) firstBin=1;
	if (lastBin==-1) lastBin=h->GetNbinsX();

	double height = getPulseHeight(h,firstBin,lastBin);
	for (Int_t i=firstBin; i<lastBin; i++){
		if (height>0){##positive pulse
			if (h->GetBinContent(i)>risePercent*height)
				return i;
		}
		else{##negative pulse
			if (h->GetBinContent(i)<risePercent*height)
				return i;
		}
	}
	##Something went wrong if we get here;
	return -1;
}

##Fit rising edge to find 0%RT
int fitRT0Bin(TH1D* h, double lowRiseFrac=0.2, double hiRiseFrac=0.5, int  firstBin=-1, int lastBin=-1){
	TF1 *flin=new TF1("flin","[0]*(x-[1])",0,1);
	double lowLim=h->GetBinCenter(getRTbin(h,lowRiseFrac,firstBin,lastBin));
	double hiLim=h->GetBinCenter(getRTbin(h,hiRiseFrac,firstBin,lastBin));

	double t00=(hiRiseFrac*lowLim-lowRiseFrac*hiLim)/(hiRiseFrac-lowRiseFrac);
	double m0=(hiRiseFrac-lowRiseFrac)*getPulseHeight(h,firstBin,lastBin)/(hiLim-lowLim);
	flin->SetRange(lowLim,hiLim);
	flin->SetParameters(m0,t00);
	h->Fit(flin,"RQC");
	return h->FindBin(flin->GetParameter(1));
}

##Find the bin that trips a 2-window trigger
##dir sets the direction to scan (1= forward, -1=reverse)
##threshCond sets whether we want to exceed (1) or fall below (-1) the threshold
##Trigger bin is the bin in the leading window which is closest to the window division. 
##   I.e. the leftmost bin of the upper window in a forward search or the rightmost bin of the lower window in the reverse.
int findTriggerBin(TH1D* h, int widthLow=5, int widthHi=5, double thresh=0.1, int limLow=-1, int limHi=-1, int dir=1, int threshCond=1, int verbose=0){
	int Nbins=h->GetNbinsX();	

	if(dir>=0){
		##Forward-going search
		if(verbose>0){cout<<"Forward"<<endl;}
		if(limLow==-1 || limLow<=widthLow){limLow=widthLow+1;}
		if(limHi==-1 || limHi>Nbins-widthHi+1){limHi=Nbins-widthHi+1;}

		##Initialize sums
		double sumLow=0, sumHi=0;
		double threshProd=thresh*widthLow*widthHi;
		for(int i=limLow-widthLow; i<limLow; i++){
			sumLow+=h->GetBinContent(i);
		}

		for(int i=limLow; i<limLow+widthHi; i++){
			sumHi+=h->GetBinContent(i);
		}

		if(verbose>0){cout << sumLow <<" "<<sumHi<<endl;}
	
		for(int i=limLow;i<=limHi;i++){
			if(verbose>0){
				cout << i << " "<<sumLow <<" "<<sumHi;
				cout <<" "<<sumHi*widthLow-sumLow*widthHi<<" ?><? "<<threshProd<<endl;
			}
			##Condition check
			if(threshCond>0 && (sumHi*widthLow-sumLow*widthHi)>threshProd){
				return i;
			}else if(threshCond<0 && (sumHi*widthLow-sumLow*widthHi)<threshProd){
				return i;
			}
			##Increment sums
			double oldLow=h->GetBinContent(i-widthLow);
			double oldTrig=h->GetBinContent(i);
			double newHi=h->GetBinContent(i+widthHi);
			sumLow=sumLow-oldLow+oldTrig;
			sumHi=sumHi-oldTrig+newHi;
		
		}
	}else{
		##Backwards-going search
		if(verbose>0){cout<<"Reverse"<<endl;}
		if(limLow==-1 || limLow<widthLow){limLow=widthLow;}
		if(limHi==-1 || limHi>Nbins-widthHi){limHi=Nbins-widthHi;}

		##Initialize sums
		double sumLow=0, sumHi=0;
		double threshProd=thresh*widthLow*widthHi;
		for(int i=limHi-widthLow+1; i<=limHi; i++){
			sumLow+=h->GetBinContent(i);
		}

		for(int i=limHi+1; i<=limHi+widthHi; i++){
			sumHi+=h->GetBinContent(i);
		}

		if(verbose>0){cout << sumLow <<" "<<sumHi<<endl;}
	
		for(int i=limHi;i>=limLow;i--){
			if(verbose>0){
				cout << i << " "<<sumLow <<" "<<sumHi;
				cout <<" "<<sumLow*widthHi-sumHi*widthLow<<" ?><? "<<threshProd<<endl;
			}
			##Condition check
			if(threshCond>0 && (sumLow*widthHi-sumHi*widthLow)>threshProd){
				return i;
			}else if(threshCond<0 && (sumLow*widthHi-sumHi*widthLow)<threshProd){
				return i;
			}
			##Increment sums
			double oldHi=h->GetBinContent(i+widthHi);
			double oldTrig=h->GetBinContent(i);
			double newLow=h->GetBinContent(i-widthLow);
			sumLow=sumLow-oldTrig+newLow;
			sumHi=sumHi+oldTrig-oldHi;
		
		}
	}
	##Didn't find it
	return -1;
}

##Find the max value when using a 2-window trigger
##This is the value that get's compared to threshold
##dir sets the direction to scan (1= forward, -1=reverse)
##Trigger bin is the bin in the leading window which is closest to the window division. 
##   I.e. the leftmost bin of the upper window in a forward search or the rightmost bin of the lower window in the reverse.
double getTrigMax(TH1D* h, int widthLow=5, int widthHi=5, int limLow=-1, int limHi=-1, int verbose=0){
	int Nbins=h->GetNbinsX();	
	double trig=0, trigMax=0;

	if(limLow==-1 || limLow<=widthLow){limLow=widthLow+1;}
	if(limHi==-1 || limHi>Nbins-widthHi+1){limHi=Nbins-widthHi+1;}

	##Initialize sums
	double sumLow=0, sumHi=0;
	for(int i=limLow-widthLow; i<limLow; i++){
		sumLow+=h->GetBinContent(i);
	}

	for(int i=limLow; i<limLow+widthHi; i++){
		sumHi+=h->GetBinContent(i);
	}

	if(verbose>0){cout << sumLow <<" "<<sumHi<<endl;}

	for(int i=limLow;i<=limHi;i++){
		trig=sumHi/widthHi-sumLow/widthLow;
		if(verbose>0){cout << i << ", "<<sumLow <<", "<<sumHi<<", "<<trig<<endl;}
		##Condition check
		trigMax=trig>trigMax?trig:trigMax;
		##Increment sums
		double oldLow=h->GetBinContent(i-widthLow);
		double oldTrig=h->GetBinContent(i);
		double newHi=h->GetBinContent(i+widthHi);
		sumLow=sumLow-oldLow+oldTrig;
		sumHi=sumHi-oldTrig+newHi;
	}

	return trigMax;
}

##Get the trigger level map for a 2-window trigger
##This is the value that get's compared to threshold
##dir sets the direction to scan (1= forward, -1=reverse)
##Trigger bin is the bin in the leading window which is closest to the window division. 
##   I.e. the leftmost bin of the upper window in a forward search or the rightmost bin of the lower window in the reverse.
TH1D* getTrigLvl(TH1D* h, int widthLow=5, int widthHi=5, int limLow=-1, int limHi=-1, int verbose=0){
	int Nbins=h->GetNbinsX();	
	TH1D *hTrigLvl=(TH1D*)h->Clone();
	hTrigLvl->Reset();

	if(limLow==-1 || limLow<=widthLow){limLow=widthLow+1;}
	if(limHi==-1 || limHi>Nbins-widthHi+1){limHi=Nbins-widthHi+1;}

	##Initialize sums
	double trig=0, sumLow=0, sumHi=0;
	for(int i=limLow-widthLow; i<limLow; i++){
		sumLow+=h->GetBinContent(i);
	}

	for(int i=limLow; i<limLow+widthHi; i++){
		sumHi+=h->GetBinContent(i);
	}

	if(verbose>0){cout << sumLow <<" "<<sumHi<<endl;}

	for(int i=limLow;i<=limHi;i++){
		trig=sumHi/widthHi-sumLow/widthLow;
		if(verbose>0){cout << i << ", "<<sumLow <<", "<<sumHi<<", "<<trig<<endl;}
		##Add to histo
		hTrigLvl->SetBinContent(i,trig);
		##Increment sums
		double oldLow=h->GetBinContent(i-widthLow);
		double oldTrig=h->GetBinContent(i);
		double newHi=h->GetBinContent(i+widthHi);
		sumLow=sumLow-oldLow+oldTrig;
		sumHi=sumHi-oldTrig+newHi;
	}

	return hTrigLvl;
}

##Crude check for prepulse digitization glitch
bool hasGlitch(TH1D *h){
  return h->GetMaximumBin()<9||h->GetMinimumBin()<9;  
}

##Calculate the RMS of the trace
double getRMS(TH1D *h, Int_t firstBin=-1, Int_t lastBin=-1){
	if (firstBin==-1) firstBin=1;
	if (lastBin==-1) lastBin=h->GetNbinsX();
	
	Double_t rms=0.0;
	for (Int_t i=firstBin; i<lastBin; i++){
	  rms+=h->GetBinContent(i)*h->GetBinContent(i);
	}
	return sqrt(rms/(lastBin-firstBin+1));
}

##Find the first entry in a trace larger than a threshold
##gtlt sets greater than (gtlt=1) or less than (gtlt=-1)
Int_t getThreshEntry(TH1D *h, Double_t thresh, Int_t gtlt=1,Int_t firstBin=-1, Int_t lastBin=-1){
	if (firstBin==-1) firstBin=1;
        if (lastBin==-1) lastBin=h->GetNbinsX();

	if (gtlt==1){
		for (Int_t i=firstBin; i<lastBin; i++){
			if(h->GetBinContent(i) > thresh) return i;
		}
	}else if(gtlt==-1){
		for (Int_t i=firstBin; i<lastBin; i++){
			if(h->GetBinContent(i) < thresh) return i;
		}
	
	}

	##If we get here, the threshold was not met or the inputs were bad
	return -1;
}
##Find the pulse height from a baseline subtracted trace
##Handles positive or negative traces
Double_t getPulseHeight(TH1D *h, int iStart=-1, int iEnd=-1){

	if(iStart<=0 || iStart>h->GetNbinsX()){iStart=-1;}
	if(iEnd<=0 || iEnd>h->GetNbinsX()){iEnd=-1;}
	if(iStart>iEnd){iStart=-1;iEnd=-1;}

	if(iStart==-1 && iEnd==-1){
		return (fabs(h->GetMaximum()) > fabs(h->GetMinimum()) ? h->GetMaximum() : h->GetMinimum() );
	}else{
		double min=h->GetBinContent(iStart);
		double max=min;
		for(int i=iStart;i<=iEnd;i++){
			double val=h->GetBinContent(i);
			min=val<min?val:min;
			max=val>max?val:max;
		}
		return fabs(max)>fabs(min)?max:min;		
	}
}

##Zero the first 8 bins of a histogram. This is used to eliminate the effects of a DCRC digitization glitch
void zeroFirst8(TH1D *h){
	for(int i=1;i<=8;i++)h->SetBinContent(i,0);
	return;
}

##Draw a horizontal histogram
void DrawHistHor(TH1 *h){

	##Draw histogram h horizontaly with bars
	Double_t ymin = h->GetMinimum();
	Double_t ymax = 1.05*h->GetMaximum();
	TAxis *axis   = h->GetXaxis();
	Double_t xmin = axis->GetXmin();
	Double_t xmax = axis->GetXmax();
	Int_t nbins   = axis->GetNbins();
	
	##Draw each bin as a box
	##This gives lines between each bin
	/*TBox box;
	Int_t Lcolor = h->GetLineColor();
	if (Lcolor == 0) Lcolor = 1;
	box.SetLineColor(Lcolor);
	box.SetLineWidth(h->GetLineWidth());
	box.SetFillStyle(0);
	Double_t dy;
	Double_t x1,y1,x2,y2;
	for (Int_t i=1;i<=nbins;i++) {
	  dy = axis->GetBinWidth(i);
	  x1 = gPad->GetUxmin();
	  y1 = axis->GetBinCenter(i)-0.5*dy;
	  x2 = h->GetBinContent(i);
	  y2 = axis->GetBinCenter(i)+0.5*dy;
	  box.DrawBox(x1,y1,x2,y2);
	}
	*/

	##Try to recreate the look of a hist
	
	##TLine line;
	TBox line;##Using boxes as lines make the vertices nicer
	Int_t Lcolor = h->GetLineColor();
	if (Lcolor == 0) Lcolor = 1;
	line.SetLineColor(Lcolor);
	line.SetLineWidth(h->GetLineWidth()/2);
	line.SetFillStyle(0);
	Double_t dy;
	Double_t xprev,xnow,y1now,y2now,y2prev;
	for (Int_t i=1;i<=nbins;i++) {
	  dy = axis->GetBinWidth(i);
	  xnow = h->GetBinContent(i);
	  y1now = axis->GetBinCenter(i)-0.5*dy;
	  y2now = axis->GetBinCenter(i)+0.5*dy;
	  ##Draw the edge of the bin
	  if(i!=1)line.DrawBox(xprev,y2prev,xnow,y1now);
	  ##Draw the top of the bin
	  line.DrawBox(xnow,y1now,xnow,y2now);

	  xprev=xnow;
	  y2prev=y2now;
	}

}

##Format the axis titles of a histogram
void formatHistTitles(TH1 *h, double offset=1.3){
	h->GetXaxis()->CenterTitle();
	h->GetYaxis()->CenterTitle();
	h->GetXaxis()->SetTitleOffset(offset);
	h->GetYaxis()->SetTitleOffset(offset);
	return;
}

##Write the contents of a matrix to a csv file
void writeMatrixToCSV(TMatrixD *m, TString fileName){
  ofstream file(fileName);
  if (file.is_open()){
    cout << "Saving matrix\n";
    int nR=m->GetNrows();
    int nC=m->GetNcols();

    for(int r=0;r<nR;r++){
      if(r>0){file<<'\n';}
      for(int c=0;c<nC;c++){
        if(c>0){file<<'\t';}
        file << m(r,c);
      }
    }

    file.close();
  }else{
    cout<<"Unable to open file.\n";
  }
}
##Read the contents of a matrix from a csv file
##Only works if m's shape is already correct
void readMatrixFromCSV(TMatrixD *m, TString fileName){
    ##using namespace std;
    ifstream in(fileName);

    int nR=m->GetNrows();
    int nC=m->GetNcols();

    for(int r=0;r<nR;r++){
      for(int c=0;c<nC;c++){
        in>>m(r,c);
      }
    }
  
  return;
}



##Get the cumulative distribution between binlow and binhi
##Direction of integration is set by forward

TH1D* GetCumulative(TH1D *h, int binlow=-1, int binhi=-1, Bool_t forward=true){
    TH1D *hintegrated = (TH1D*)h->Clone();
    hintegrated->SetName(h->GetName()+TString("_int"));
    hintegrated->SetTitle(h->GetTitle()+TString(" integrated"));

    if (binlow==-1) binlow=1;
    if (binhi==-1) binhi=h->GetNbinsX();
    hintegrated->Reset();

    if (forward) { ## Forward computation
       Double_t sum = 0.;
        for (Int_t binx = binlow; binx <= binhi; ++binx) {
           sum += h->GetBinContent(binx);
           hintegrated->SetBinContent(binx, sum);
        }

    } else { ## Backward computation
       Double_t sum = 0.;
        for (Int_t binx = binhi; binx >= binlow; --binx) {
           sum += h->GetBinContent(binx);
           hintegrated->SetBinContent(binx, sum);
        }
    }
    return hintegrated;
}

##Build trigger logic filter histogram
##N0 is the first bin of the Hi side
TH1D* buildTrigLogicHist(int Nsamples=4096, int widthLow=5, int widthHi=5, int N0=2048){
	TH1D *htrig = new TH1D("htrig","htrig",Nsamples,0,Nsamples);
	for(int i=N0-widthLow;i<=N0;i++){htrig->SetBinContent(i,-1.0/widthLow);}
	for(int i=N0;i<N0+widthHi;i++){htrig->SetBinContent(i,1.0/widthHi);}
	return htrig;
}

##Convert a series name to a double
double SeriesName2Double(TString name){
	##name is of the form: FFYYMMDD_HHMM
	##Where FF is facility code and the rest is date/time
	double number=0;

	name.ReplaceAll("_","");
	int L=name.Length();
	for(int i=0;i<L;i++){
		number=number+((int)name[L-1-i]-48)*(10.**i);
	}

	return number;
}

##Convert a double representation of a SeriesNumber to a string
TString Double2SeriesName(double number){
	##name should be of the form: FFYYMMDD_HHMM
	int FF=int(number/1e10);
	int YYMMDD = int((number-FF*1e10)/1e4);
	int HHMM = number-FF*1e10-YYMMDD*1e4;

	return TString(Form("%02i%06i_%04i",FF,YYMMDD,HHMM));
}

##Calculate sigmaPTOFamps from a filterfile
double calcSigmaPTOFamps(TFile *f){
	TH1D *hPSD;
	TH1D *hTempFFTRe;
	TH1D *hTempFFTIm;
	double re,im,sPn,sum=0;

	double sampleRate = 1000000/0.8;
	double fScale = 0.5*sampleRate; ##Nyquist frequency = sampling freq/2
	
	f->GetObject("zip1/PTNoisePSD",hPSD);
	f->GetObject("zip1/PTTemplateFFTRe",hTempFFTRe);
	f->GetObject("zip1/PTTemplateFFTIm",hTempFFTIm);

	int N=(int)(2*(hPSD->GetNbinsX()-1));
	double binWidth = fScale/(N/2);

	hTempFFTRe->Scale(1/(0.8e-6)**0.5);
	hTempFFTIm->Scale(1/(0.8e-6)**0.5);

	for(int i=2;i<=N/2+1;i++){
		re=hTempFFTRe->GetBinContent(i);
		im=hTempFFTIm->GetBinContent(i);
		sPn=N*binWidth*(hPSD->GetBinContent(i))**2;
		sum=sum+2*(re**2+im**2)/(sPn/2);
	}
	return 1.0/(sum**(0.5));
}

