##Set of utilities useful for analysing rq data in root
##N. Mast, 2019
##
#include <vector> #Probably don't need to replace this, likely built into python.
#If I do it'll probably be numpy or something easy.

import uproot

#Workaround for ROOT path missing
import sys
sys.path.append('/home/gerudo7/anaconda3/envs/root_env/lib')

import ROOT
import math


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
    print(filepath+" does not exist.\n")
  zipName = "rqDir/zip"+zip

  z = ROOT.TChain(zipName)
  e = ROOT.TChain("rqDir/eventTree")
  z.AddFriend(e)
  filepath += "/*.root"
  z.Add(filepath)
  e.Add(filepath)
  return z

#######################
##Load RQ's into a TChain
#######################
def loadSeriesRQFull(filepath, zip="1"):
  zipName = "rqDir/zip"+zip
  z = ROOT.TChain(zipName)
  e = ROOT.TChain("rqDir/eventTree")
  z.AddFriend(e)
  filepath.Append("*.root")
  z.Add(filepath)
  e.Add(filepath)
  return z

#######################
##Save TCut to file
#######################
def saveCuts(cuts, fileName):
  try:
    file = open(fileName,'wt')
    print("Saving cuts: " + cuts.GetTitle() +'\n')
    file.write(cuts.GetTitle() + '\n')
    file.close()
  except:
    print("Unable to open or write file.\n")

#######################
##Load TCut from File
#######################
def loadCuts(fileName):
  try:
    file = open(fileName)
    cutString = file.read #get data from file and save it to `line`
    print(cutString + '\n')
    file.close()
    return cutString.Data()
  except:
    print("Unable to open or read file.\n")
    return ""

#######################
##Save eventlist to file
#######################
def eventListToFile(z, cuts, fileName, maxEvents):
  print("Writing eventList to file.\n")

  z.Draw(">>elist",cuts)
  elist = ROOT.gDirectory.Get("elist")

  entryIndex=0
  #SeriesNumberCurr,SeriesNumberNew, EventNumber
  prevEventNumber=0
  rawFileNameCurr = "{:<100}".format("")
  file = open(fileName,'wt')
  eventCount=0
  i=0

  entryIndex=elist.GetEntry(i)
  while (entryIndex !=-1) & (eventCount<maxEvents):
    z.GetEntry(entryIndex)
    SeriesNumberNew=z.GetLeaf("SeriesNumber").GetValue()
    rawFileNameNew = z.GetFile().GetName()

    ##Single File for all events, separated by series,dumps
    if (rawFileNameNew != rawFileNameCurr):
      file.write("**********\n")
      file.write(rawFileNameNew + '\n')
      rawFileNameCurr = rawFileNameNew.value() #trying to copy the value over without telling it that they are the same variable.
      
    EventNumber=z.GetLeaf("EventNumber").GetValue()
    ##Only take 1 of every 10 adjascent events to adef getting retriggers
    if 1 or abs(EventNumber-prevEventNumber)>10:
      prevEventNumber=EventNumber
      file.write(ROOT.Form("%.0f\n",EventNumber))
      eventCount += 1

    i += 1
    entryIndex=elist.GetEntry(i)
  file.close()
  print("Wrote " + eventCount + " events to: " + fileName + '\n')
  return eventCount

#######################
##Save RQs to json file
#######################
def RQsToJSON(z, cuts, fileName, RQList, maxEvents):
  ##This is a cumbersome way to do this, but I can't figure out a better one since you can't
  ##  access values of calculated aliases as you would a normal branch

  print("Collecting RQs from Tree.\n")
  RQs = None #initialize multidimensional vector? is that what this does?

  ##Generate arrays of desired RQs ordered by Tree Entry$ number
  it = RQList.begin()
  while(it != RQList.end()):
    Nentries=z.Draw(ROOT.Form("%s%s",it.Data(),":Entry$"),cuts,"goff")
    if(Nentries==-1):
      ##Got an error drawing
      print("Error, bad RQ. \n")
      return -1
    
    RQarr=z.GetV1()##RQ

    RQtemp = []
    i = 0
    while(i<Nentries):
      RQtemp.push_back(RQarr[i])
      i += 1
    RQs.push_back(RQtemp)

    it += 1

  ##Loop through arrays, building and outputting a json-style entry for each Tree entry
  print("Writing RQs to JSON.\n")

  file = open(fileName,'wt')
  eventCount=0

  ##Begin overall library
  file.write("{ \n")
  ##Begin main pulse data list
  file.write("\"data\":[ \n")

  iEvent = 0
  while(iEvent<RQs[0].size()):
    ##Write the current event     {"RQ1":071705231956,"RQ2":10001},
    if(iEvent>0):
      file.write(", \n")##comma after the previous event's block
    file.write("\t{")##Start this event's block
    iRQ = 0
    while(iRQ < RQList.size()):
      if(iRQ>0):
        file.write(",")##Comma after previous RQ

      if(RQList[iRQ]=="SeriesNumber"):
        file.write("\"" + RQList[iRQ].Data() + "\":" + ROOT.Form("%.0f",RQs[iRQ][iEvent]))
      else:
        file.write("\"" + RQList[iRQ].Data() + "\":" + RQs[iRQ][iEvent])

      iRQ += 1

    file.write("}\n")##End this event's block
    eventCount += 1
    if(eventCount>=maxEvents):
      break
    iEvent += 1
  
  ##Close data list
  file.write("\n ] \n")
  ##Close overall library
  file.write("\n } \n")

  print("Closing file \n")
  ##close the raw data file
  file.close()

  print("Wrote " + eventCount + " events to: " + fileName + '\n')
  return eventCount

#################################################
##Begin PulseTools-like functions
#################################################

##Remove a trace's baseline
def BSsubtract(h, bs):
  i = 1
  while (i<=h.GetNbinsX()):
    h.AddBinContent(i,-bs)
    i += 1 

##Remove a trace's baseline
##Define pre and post numbers of bins to use in calculateing baseline
##Removes linear background based on beginning and ending of trace
def BSsubtractSlope(h, iStart, nPre, nPost):
  nBins = h.GetNbinsX()
  
  if(nPost>0):##Do the sloped BS subtraction
    pre = getBSpre(h, iStart, nPre)
    post = getBSpost(h, nPost)
    m = (post-pre)/(nBins-nPre/2-nPost/2)
    b = ((nBins-nPost/2)*pre-(nPre/2)*post) / (nBins-nPre/2-nPost/2)
    i = 1
    while(i<=nBins):
      val=h.GetBinContent(i)
      val=val-(m*i+b)
      h.SetBinContent(i,val)
      i += 1
    
  else:##Just do the flat, pre pulse BS subtraction
    BSsubtract(h, getBSpre(h,iStart, nPre))

##Rescale a trace's x (time) axis
##h is the existing trace, hNew will contain the rescaled version
##Bin centers will start with t0 and be separated by dt.
def rescaleX(h, t0, dt):
  nBins = h.GetNbinsX()
  hNew = ROOT.TH1D(h.GetName(),h.GetTitle(),nBins,t0-dt/2,t0+dt*(nBins-0.5))
  i = 1
  while(i<=nBins):
    hNew.SetBinContent(i,h.GetBinContent(i))
    i += 1
  return hNew

##Rescale a trace's x axis
##h is the existing trace, hNew will contain the rescaled version
##scale x bins by factor f
##Assumes constant bin widths
def rescaleX(h, f):
  nBins = h.GetNbinsX()
  low = f*h.GetBinLowEdge(1)
  hi = f*h.GetBinLowEdge(nBins+1)
  hNew = ROOT.TH1D(h.GetName(),h.GetTitle(),nBins,low,hi)
  i = 1
  while(i<=nBins):
    hNew.SetBinContent(i,h.GetBinContent(i))
    i += 1
  return hNew


##calculate a trace's beginning baseline
def getBSpre(h, binLow, binHi):
  sum=0
  i = binLow
  while (i<=binHi):
    sum=sum+h.GetBinContent(i)
    i += 1
  
  return sum/(binHi-binLow+1)

##calculate a trace's ending baseline
def getBSpost(h, Nbins):
  sum=0
  i = h.GetNBinsX()-Nbins+1
  while (i<=h.GetNbinsX()):
    sum=sum+h.GetBinContent(i)
    i += 1
  return sum/Nbins

##Crude low-pass filter
##tau, tStart, tEnd have same units as x-axis of h
##Use average of back and forth filtering
def lpf(h, tau, tStart, tEnd):

  dt=h.GetXaxis().GetBinCenter(2)-h.GetXaxis().GetBinCenter(1)
  alpha = dt/tau
  htemp = h.Clone()

  iStart=h.FindBin(tStart)
  iEnd=h.FindBin(tEnd)

  if(iStart<=0):
    iStart=1
  if(iEnd<=0) or (iEnd>h.GetNbinsX()):
    iEnd=h.GetNbinsX()

  ##Forward pass
  yk1=h.GetBinContent(iStart)
  i = iStart
  while(i<=iEnd):
    yk1 = yk1 + alpha*(h.GetBinContent(i)-yk1)
    htemp.SetBinContent(i,yk1)
    i += 1
  
  ##Backwards pass
  yk2=h.GetBinContent(h.GetNbinsX())
  i = iEnd
  while(i>=iStart):
    yk2 = yk2 + alpha*(h.GetBinContent(i)-yk2)
    yk1=htemp.GetBinContent(i)
    h.SetBinContent(i,(yk1+yk2)/2)
    i -= 1
  return

##using darkpipe symmetric convention for normalization
def RealToComplexFFT( pulsevector, outComp):
  n = pulsevector.size()
  if(n == 0):
    print("RealToComplexFFT - ERROR! empty pulse passed into this function \n")
    exit(1)

    pvector = []#double[pulsevector.size()]
    re,im = None#re,im

  ##copying std vector back into array to pass to FFT (to mirror what is done in PulseTools)
  i = 0
  while i<pulsevector.size():
    pvector[i] = pulsevector[i]
    i += 1

  fftr2c = ROOT.TVirtualFFT.FFT(1,n,"R2C ES")
  fftr2c.SetPoints(pvector)
  fftr2c.Transform()
  fftr2c.GetPointComplex(0,re,im)
  tempCopy = ROOT.TComplex.operator(re/(n**(1/2)), im/(n**(1/2)))
  outComp.push_back(tempCopy)

  i = 1
  while(i<n):
    fftr2c.GetPointComplex(i,re,im)
    tempCopy = ROOT.TComplex.operator(re/(n**(1/2)), im/(n**(1/2)))
    outComp.push_back(tempCopy)
    i += 1

  ##done! so delete new'd objects
  del pvector
  del fftr2c

  return

def ComplexToRealIFFT(inComp, outRe):
    n = inComp.size()
    if(n == 0):
      print("PulseTools::ComplexToRealIFFT - ERROR! empty pulse passed into this function \n")
      exit(1)

    re_pvector = [] #double[inComp.size()]
    im_pvector = [] #double[inComp.size()]

    re,im = None

    ##copying std vector back into arrays to pass to IFFT
    i = 0
    while (i < n):
       re_pvector[i] = inComp[i].Re()
       im_pvector[i] = inComp[i].Im()
       i += 1

    ifftc2r = ROOT.TVirtualFFT.FFT(1, n,"C2R ES") ##this should be the opposite of FFT
    ifftc2r.SetPointsComplex(re_pvector, im_pvector)
    ifftc2r.Transform()
    ifftc2r.GetPointComplex(0,re,im)
    outRe.push_back(re/(n**(1/2))) ##DC component

    i = 1
    while(i<n):
       ifftc2r.GetPointComplex(i,re, im)
       outRe.push_back(re/(n**(1/2)))
       i += 1

    del re_pvector
    del im_pvector
    del ifftc2r

    return


#################################################
##End PulseTools
#################################################
def DrawTraces(hArr, title, nTraces=9, xmin=0, xmax=3.27640000000000000e+06):
  if (nTraces>9):
    print("Not ready for more than 9 traces \n")
    return

  goodTraces=[0,1,2,3,4,5,6,7,10]

  color=[ROOT.EColor.kCyan,ROOT.EColor.kCyan+2,ROOT.EColor.kBlue,ROOT.EColor.kBlue+2,ROOT.EColor.kMagenta,ROOT.EColor.kMagenta+2,ROOT.EColor.kRed,ROOT.EColor.kRed+2,ROOT.EColor.kBlack]
  ##Set the colors and find the largest amplitude
  max=0
  i = 0
  while(i<nTraces):
    index=goodTraces[i]
    if (max<hArr[index].GetMaximum()):
      max=hArr[index].GetMaximum()
    hArr[index].SetLineColor(color[i])
    i += 1
  
  ##Draw it
  hArr[0].GetYaxis().SetRangeUser(-100,100*int(math.ceil(max/100.0)))
  hArr[0].GetXaxis().SetRangeUser(xmin,xmax)
  hArr[0].GetXaxis().SetTitle("time (ns)")
  hArr[0].SetTitle(title)
  hArr[0].Draw()
  i = 1
  while(i<nTraces):
    index=goodTraces[i]
    hArr[index].Draw("same")
    i += 1
  
  leg = ROOT.TLegend(0.7, 0.5, 0.9, 0.9)
  i = 0
  while(i<9):
    index=goodTraces[i]
    leg.AddEntry(rawPulse[index],rawPulse[index].GetName())
    i += 1
  
  leg.Draw()
  
  return

##Shift a trace by some number of bins, appending 0 signal where needed
##Only use this after removing the baseline from a trace
def shiftTrace(h, shift):
  if (shift>0):##Moving the pulse later in time.
    ##Shift contents towards end of trace
    i = h.GetNbinsX()
    while(i>shift):
      h.SetBinContent(i,h.GetBinContent(i-shift))
      i -= 1
    
    ##Zero the beginning of the trace
    i = shift
    while(i>0):
      h.SetBinContent(i,0)
      i -= 1
  else:##Move pulse earlier in time
    ##Shift contents toward beginning of trace
    #value
    i = 1
    while(i<h.GetNbinsX()-abs(shift)):
      value=h.GetBinContent(i+abs(shift))
      h.SetBinContent(i,value)
      i += 1
    ##Zero the end of the trace
    i=h.GetNbinsX()-abs(shift)
    while(i<=h.GetNbinsX()):
      h.SetBinContent(i,0)
      i += 1
  return

##Find the 50%RT bin
##Use only on bs subtracted traces
def getRTbin(h, risePercent=0.5, firstBin=-1, lastBin=-1):
  if (firstBin==-1): firstBin=1
  if (lastBin==-1): lastBin=h.GetNbinsX()
  
  height = getPulseHeight(h,firstBin,lastBin)
  i=firstBin
  while (i<lastBin):
    if (height>0):##positive pulse
      if (h.GetBinContent(i)>risePercent*height):
        return i
    else:##negative pulse
      if (h.GetBinContent(i)<risePercent*height):
        return i
    i += 1
  ##Something went wrong if we get here
  return -1

##Fit rising edge to find 0%RT
def fitRT0Bin(h, lowRiseFrac=0.2, hiRiseFrac=0.5, firstBin=-1, lastBin=-1):
  flin = ROOT.TF1("flin","[0]*(x-[1])",0,1)
  lowLim=h.GetBinCenter(getRTbin(h,lowRiseFrac,firstBin,lastBin))
  hiLim=h.GetBinCenter(getRTbin(h,hiRiseFrac,firstBin,lastBin))

  t00=(hiRiseFrac*lowLim-lowRiseFrac*hiLim)/(hiRiseFrac-lowRiseFrac)
  m0=(hiRiseFrac-lowRiseFrac)*getPulseHeight(h,firstBin,lastBin)/(hiLim-lowLim)
  flin.SetRange(lowLim,hiLim)
  flin.SetParameters(m0,t00)
  h.Fit(flin,"RQC")
  return h.FindBin(flin.GetParameter(1))

##Find the bin that trips a 2-window trigger
##dir sets the direction to scan (1= forward, -1=reverse)
##threshCond sets whether we want to exceed (1) or fall below (-1) the threshold
##Trigger bin is the bin in the leading window which is closest to the window division. 
##   I.e. the leftmost bin of the upper window in a forward search or the rightmost bin of the lower window in the reverse.
def findTriggerBin(h, widthLow=5, widthHi=5, thresh=0.1, limLow=-1, limHi=-1, dir=1, threshCond=1, verbose=0):
  Nbins=h.GetNbinsX()  

  if(dir>=0): ##Forward-going search
    if(verbose>0):
      print("Forward \n")
    if(limLow==-1) or (limLow<=widthLow):
      limLow=widthLow+1
    if(limHi==-1) or (limHi>Nbins-widthHi+1):
      limHi=Nbins-widthHi+1

##Initialize sums
    sumLow=0
    sumHi=0
    threshProd=thresh*widthLow*widthHi
    i = limLow-widthLow
    while(i<limLow):
      sumLow+=h.GetBinContent(i)
      i += 1

    i = limLow
    while(i<limLow+widthHi):
      sumHi+=h.GetBinContent(i)
      i += 1

    if(verbose>0):
      print( + sumLow + " " + sumHi + '\n')

    i = limLow
    while(i<=limHi):
      if(verbose>0):
        print(i + " " + sumLow + " " + sumHi)
        print(" " + (sumHi*widthLow-sumLow*widthHi) + " ?><? " + threshProd + '\n')
      ##Condition check
      if(threshCond>0) and ((sumHi*widthLow-sumLow*widthHi)>threshProd):
        return i
      elif(threshCond<0) and ((sumHi*widthLow-sumLow*widthHi)<threshProd):
        return i

      ##Increment sums
      oldLow=h.GetBinContent(i-widthLow)
      oldTrig=h.GetBinContent(i)
      newHi=h.GetBinContent(i+widthHi)
      sumLow=sumLow-oldLow+oldTrig
      sumHi=sumHi-oldTrig+newHi
    
  else:
    ##Backwards-going search
    if(verbose>0):
      print("Reverse \n")
    if(limLow==-1) or (limLow<widthLow):
      limLow=widthLow
    if(limHi==-1) or (limHi>Nbins-widthHi):
      limHi=Nbins-widthHi

    ##Initialize sums
    sumLow=0
    sumHi=0
    threshProd=thresh*widthLow*widthHi
    i = limHi-widthLow+1
    while(i<=limHi):
      sumLow+=h.GetBinContent(i)
      i += 1

    i = limHi+1
    while(i<=limHi+widthHi):
      sumHi+=h.GetBinContent(i)
      i += 1

    if(verbose>0):
      print(sumLow + " " + sumHi + '\n')
  
    i = limHi
    while(i>=limLow):
      if(verbose>0):
        print(i + " " + sumLow + " " + sumHi)
        print(" " + sumLow*widthHi-sumHi*widthLow + " ?><? " + threshProd + '\n')

      ##Condition check
      if(threshCond>0) and ((sumLow*widthHi-sumHi*widthLow)>threshProd):
        return i
      elif(threshCond<0) and ((sumLow*widthHi-sumHi*widthLow)<threshProd):
        return i
      ##Increment sums
      oldHi=h.GetBinContent(i+widthHi)
      oldTrig=h.GetBinContent(i)
      newLow=h.GetBinContent(i-widthLow)
      sumLow=sumLow-oldTrig+newLow
      sumHi=sumHi+oldTrig-oldHi
      i -= 1
  ##Didn't find it
  return -1

##Find the max value when using a 2-window trigger
##This is the value that gets compared to threshold
##dir sets the direction to scan (1= forward, -1=reverse)
##Trigger bin is the bin in the leading window which is closest to the window division. 
##   I.e. the leftmost bin of the upper window in a forward search or the rightmost bin of the lower window in the reverse.
def getTrigMax(h, widthLow=5, widthHi=5, limLow=-1, limHi=-1, verbose=0):
  Nbins=h.GetNbinsX()
  trig=0
  trigMax=0
  
  if(limLow==-1) or (limLow<=widthLow):
    limLow=widthLow+1
  if(limHi==-1) or (limHi>Nbins-widthHi+1):
    limHi=Nbins-widthHi+1

  ##Initialize sums
  sumLow=0
  sumHi=0
  i=limLow-widthLow
  while(i<limLow):
    sumLow += h.GetBinContent(i)
    i += 1

  i = limLow
  while(i<limLow+widthHi):
    sumHi+=h.GetBinContent(i)
    i += 1

  while(verbose>0):
    print(sumLow + " " + sumHi + '\n')

  i = limLow
  while(i<=limHi):
    trig=sumHi/widthHi-sumLow/widthLow
    if(verbose>0): 
      print(i + ", " + sumLow + ", " + sumHi + ", " + trig + '\n') ##Condition check
    if trig > trigMax: trigMax = trig
    #trigMax=trig>trigMax?trig:trigMax
    ##Increment sums
    oldLow=h.GetBinContent(i-widthLow)
    oldTrig=h.GetBinContent(i)
    newHi=h.GetBinContent(i+widthHi)
    sumLow=sumLow-oldLow+oldTrig
    sumHi=sumHi-oldTrig+newHi
    i += 1

  return trigMax

##Get the trigger level map for a 2-window trigger
##This is the value that gets compared to threshold
##dir sets the direction to scan (1= forward, -1=reverse)
##Trigger bin is the bin in the leading window which is closest to the window division. 
##   I.e. the leftmost bin of the upper window in a forward search or the rightmost bin of the lower window in the reverse.
def getTrigLvl(h, widthLow=5, widthHi=5, limLow=-1, limHi=-1, verbose=0):
  Nbins=h.GetNbinsX()
  hTrigLvl=h.Clone()
  hTrigLvl.Reset()

  if(limLow==-1) or (limLow<=widthLow):
    limLow=widthLow+1
  if(limHi==-1) or (limHi>Nbins-widthHi+1):
    limHi=Nbins-widthHi+1

  ##Initialize sums
  trig=0
  sumLow=0
  sumHi=0
  i=limLow-widthLow
  while(i<limLow):
    sumLow+=h.GetBinContent(i)
    i += 1

  i=limLow
  while(i<limLow+widthHi):
    sumHi+=h.GetBinContent(i)
    i += 1

  if(verbose>0):
    print(sumLow + " " + sumHi + '\n')

  i=limLow
  while(i<=limHi):
    trig=sumHi/widthHi-sumLow/widthLow
    if(verbose>0):
      print(i + ", " + sumLow + ", " + sumHi + ", " + trig + '\n')
    ##Add to hist
    hTrigLvl.SetBinContent(i,trig)
    ##Increment sums
    oldLow=h.GetBinContent(i-widthLow)
    oldTrig=h.GetBinContent(i)
    newHi=h.GetBinContent(i+widthHi)
    sumLow=sumLow-oldLow+oldTrig
    sumHi=sumHi-oldTrig+newHi

  return hTrigLvl

##Crude check for prepulse digitization glitch
def hasGlitch(h):
  return h.GetMaximumBin()<9 or (h.GetMinimumBin()<9)

##Calculate the RMS of the trace
def getRMS(h, firstBin=-1, lastBin=-1):
  if (firstBin==-1): firstBin=1
  if (lastBin==-1): lastBin=h.GetNbinsX()
  
  rms=0.0
  i = firstBin
  while (i<lastBin):
    rms += h.GetBinContent(i)*h.GetBinContent(i)
    i += 1
  return (rms/(lastBin-firstBin+1))**(1/2)

##Find the first entry in a trace larger than a threshold
##gtlt sets greater than (gtlt=1) or less than (gtlt=-1)
def getThreshEntry(h, thresh, gtlt=1,firstBin=-1, lastBin=-1):
  if (firstBin==-1): firstBin=1
  if (lastBin==-1): lastBin=h.GetNbinsX()

  if (gtlt==1):
    i = firstBin
    while (i<lastBin):
      if(h.GetBinContent(i) > thresh): return i
      i += 1
  elif(gtlt==-1):
    i = firstBin
    while (i<lastBin):
      if(h.GetBinContent(i) < thresh): return i
      i += 1

  ##If we get here, the threshold was not met or the inputs were bad
  return -1

##Find the pulse height from a baseline subtracted trace
##Handles positive or negative traces
def getPulseHeight(h, iStart=-1, iEnd=-1):

  if(iStart<=0) or (iStart>h.GetNbinsX()):
    iStart=-1
  if(iEnd<=0) or (iEnd>h.GetNbinsX()):
    iEnd=-1
  if(iStart>iEnd):
    iStart=-1
    iEnd=-1

  if(iStart==-1) and (iEnd==-1):
    if abs(h.GetMaximum()) > abs(h.GetMinimum()):
      return h.GetMaximum()
    else:
      return h.GetMinimum
  else:
    min=h.GetBinContent(iStart)
    max=min
    i = iStart
    while(i<=iEnd):
      val=h.GetBinContent(i)
      if val<min:
        min = val
      if val>max:
        max = val
      i += 1
    if abs(max)>abs(min):
      return max
    else:
      return min

##Zero the first 8 bins of a histogram. This is used to eliminate the effects of a DCRC digitization glitch
def zeroFirst8(h):
  i = 1
  while(i<=8):
    h.SetBinContent(i,0)
    i += 1
  return

##Draw a horizontal histogram
def DrawHistHor(h):

  ##Draw histogram h horizontally with bars
  ymin = h.GetMinimum()
  ymax = 1.05*h.GetMaximum()
  axis   = h.GetXaxis()
  xmin = axis.GetXmin()
  xmax = axis.GetXmax()
  nbins   = axis.GetNbins()
  
  ##Draw each bin as a box
  ##This gives lines between each bin

  ##Try to recreate the look of a hist
  
  ##TLine line
  line = ROOT.TBox##Using boxes as lines make the vertices nicer
  Lcolor = h.GetLineColor()
  if (Lcolor == 0): Lcolor = 1
  line.SetLineColor(Lcolor)
  line.SetLineWidth(h.GetLineWidth()/2)
  line.SetFillStyle(0)
  #dy = None
  #xprev,xnow,y1now,y2now,y2prev = None
  i = 1
  while (i<=nbins):
    dy = axis.GetBinWidth(i)
    xnow = h.GetBinContent(i)
    y1now = axis.GetBinCenter(i)-0.5*dy
    y2now = axis.GetBinCenter(i)+0.5*dy
    ##Draw the edge of the bin
    if(i!=1): line.DrawBox(xprev,y2prev,xnow,y1now)
    ##Draw the top of the bin
    line.DrawBox(xnow,y1now,xnow,y2now)

    xprev=xnow
    y2prev=y2now

    i += 1

##Format the axis titles of a histogram
def formatHistTitles(h, offset=1.3):
  h.GetXaxis().CenterTitle()
  h.GetYaxis().CenterTitle()
  h.GetXaxis().SetTitleOffset(offset)
  h.GetYaxis().SetTitleOffset(offset)
  return

##Write the contents of a matrix to a csv file
def writeMatrixToCSV(m, fileName):
  file = open(fileName,'wt')
  if (file.is_open()):
    print("Saving matrix\n")
    nR=m.GetNrows()
    nC=m.GetNcols()

    r = 0
    while(r<nR):
      if(r>0): file.write('\n')
      c = 0
      while(c<nC):
        if(c>0): file.write('\t')
        file.write(m(r,c))
        c += 1
      r += 1

    file.close()
  else:
    print("Unable to open file.\n")

##Read the contents of a matrix from a csv file
##Only works if m's shape is already correct
def readMatrixFromCSV(m, fileName):
  ##using namespace std
  infile = open(fileName,'rt')

  nR=m.GetNrows()
  nC=m.GetNcols()

  r = 0
  while(r<nR):
    c = 0
    while(c<nC):
      infile.read(m(r,c))
      c += 1
    r += 1

  return



##Get the cumulative distribution between binlow and binhi
##Direction of integration is set by forward

def GetCumulative(h, binlow=-1, binhi=-1, forward=True):
    hintegrated = h.Clone()
    hintegrated.SetName(h.GetName()+ROOT.TString("_int"))
    hintegrated.SetTitle(h.GetTitle()+ROOT.TString(" integrated"))

    if (binlow==-1): binlow=1
    if (binhi==-1): binhi=h.GetNbinsX()
    hintegrated.Reset()

    if (forward): ## Forward computation
      sum = 0.

      binx = binlow 
      while (binx <= binhi):
        sum += h.GetBinContent(binx)
        hintegrated.SetBinContent(binx, sum)
        binx += 1

    else: ## Backward computation
      sum = 0.
      binx = binhi
      while (binx >= binlow):
        sum += h.GetBinContent(binx)
        hintegrated.SetBinContent(binx, sum)
        binx -= 1

    return hintegrated

##Build trigger logic filter histogram
##N0 is the first bin of the Hi side
def buildTrigLogicHist(Nsamples=4096, widthLow=5, widthHi=5, N0=2048):
  htrig = ROOT.TH1D("htrig","htrig",Nsamples,0,Nsamples)
  i = N0-widthLow
  while(i<=N0):
    htrig.SetBinContent(i,-1.0/widthLow)
    i += 1
  i = N0
  while(i<N0+widthHi):
    htrig.SetBinContent(i,1.0/widthHi)
    i += 1
  return htrig

##Convert a series name to a double
def SeriesName2Double(name):
  ##name is of the form: FFYYMMDD_HHMM
  ##Where FF is facility code and the rest is date/time
  number=0

  name.ReplaceAll("_","")
  L=name.Length()
  i = 0
  while(i<L):
    number=number+(name[L-1-i]-48)*(10.**i)
    i += 1

  return number

##Convert a representation of a SeriesNumber to a string
def Double2SeriesName(number):
  ##name should be of the form: FFYYMMDD_HHMM
  FF=int(number/1e10)
  YYMMDD = int((number-FF*1e10)/1e4)
  HHMM = number-FF*1e10-YYMMDD*1e4

  return ROOT.TString(ROOT.Form("%02i%06i_%04i",FF,YYMMDD,HHMM))

##Calculate sigmaPTOFamps from a filterfile
def calcSigmaPTOFamps(f):
  hPSD = None
  hTempFFTRe = None
  hTempFFTIm = None
  re,im,sPn,sum=0

  sampleRate = 1000000/0.8
  fScale = 0.5*sampleRate ##Nyquist frequency = sampling freq/2
  
  f.GetObject("zip1/PTNoisePSD",hPSD)
  f.GetObject("zip1/PTTemplateFFTRe",hTempFFTRe)
  f.GetObject("zip1/PTTemplateFFTIm",hTempFFTIm)

  N=(int)(2*(hPSD.GetNbinsX()-1))
  binWidth = fScale/(N/2)

  hTempFFTRe.Scale(1/(0.8e-6)**0.5)
  hTempFFTIm.Scale(1/(0.8e-6)**0.5)

  i = 2
  while(i<=N/2+1):
    re=hTempFFTRe.GetBinContent(i)
    im=hTempFFTIm.GetBinContent(i)
    sPn=N*binWidth*(hPSD.GetBinContent(i))**2
    sum=sum+2*(re**2+im**2)/(sPn/2)
    i += 1
  return 1.0/(sum**(0.5))