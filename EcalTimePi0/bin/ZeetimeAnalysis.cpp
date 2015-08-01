#include <iostream>
#include <string>

#include <map>
#include <vector>
#include <algorithm> 
#include <functional>
#include <set>
#include <boost/tokenizer.hpp>

#include <iostream>
#include <math.h> 
#include <assert.h>

#include "FWCore/PluginManager/interface/PluginManager.h"
#include "FWCore/PluginManager/interface/standard.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Math/interface/LorentzVector.h"

#include "CalibCalorimetry/EcalTiming/interface/EcalTimeTreeContent.h"
#include "ECALTime/EcalTimePi0/interface/timeVsAmpliCorrector.h"
#include "ECALTime/EcalTimePi0/interface/EcalObjectTimeCalibCalorimetryTree.h"

#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TF1.h"


typedef std::set<std::pair<int,int> > SetOfIntPairs;

// authors: S. Cooper and G. Franzoni (UMN)
// Maitained by Author: Tambe E. Norbert (UMN)
#define BarrelLimit  1.479
#define EndcapLimit  3.0

#define ADCtoGeVEB   0.039
#define ADCtoGeVEE   0.063

#define numAeffBins     35
#define numAoSigmaBins  25

#define NSlices 54   // # of slices in log(A)
#define LogStep 0.05 // size of Slices log(A)

#define lightSpeed 299792458



// -------- Globals ----------------------------------------
EcalTimeTreeContent treeVars_; 
TFile* saving_;
std::vector<std::string> listOfFiles_;
bool speak_=false;
char buffer_ [75];
std::string bufferTitle_; 
// default settings
std::string outputRootName_ = "outputHistos.root";
int   numEvents_      = -1;
unsigned int  minRun_ = 0;
unsigned int  maxRun_ = 9999999;
unsigned int  minLS_ = 0;
unsigned int  maxLS_ = 9999999;
float eTGammaMinEB_   = 0.2;
float s4s9GammaMinEB_ = 0.85;
float eTPi0MinEB_     = 0.65;
float eTGammaMinEE_   = 0.250;
float s4s9GammaMinEE_ = 0.85;
float eTPi0MinEE_     = 0.800;
float swissCrossMaxEB_ = 0.95; // 1-E4/E1
float swissCrossMaxEE_ = 0.95; // 1-E4/E1
// based on range and bins, the bin width is 50 ps
float rangeTDistro_ = 10; 
int   binsTDistro_  = 120; 
std::vector<std::vector<double> > trigIncludeVector;
std::vector<std::vector<double> > trigExcludeVector;
std::vector<std::vector<double> > ttrigIncludeVector;
std::vector<std::vector<double> > ttrigExcludeVector;



int  minEntriesForFit_ = 7;
int  flagOneVertex_ = 0;
bool limitFit_(true); 
//std::string fitOption_(""); // default: use chi2 method
std::string fitOption_("L"); // use likelihood method

// Ao dependent timing corrections
// By the definition of these corrections, the timing should be zero for the hits in Module 1 
// or Low eta EE within the valid A/sigma ranges.  Earlier data will have positive time due 
// to the graduate timing shifts in the positive direction.
TF1* timeCorrectionEB_ = new TF1("timeCorrectionEB_","pol4(0)",0,1.2);
TF1* timeCorrectionEE_ = new TF1("timeCorrectionEE_","pol4(0)",0,1.2);


// -------- Histograms -------------------------------------
TH1 * nVertices_;
TH1F* mass_;
TH1F* dZvertices_;
TH1F* Zvertices_;
TH1F* gammaEt1_, *gammaEt2_;
TH1F* numSclusters;
// ---------------------------------------------------------------------------------------
// - Function to decide to include/exclude event based on the vectors passed for triggers 
bool includeEvent(int* triggers,
    int numTriggers,
    std::vector<std::vector<double> > includeVector,
    std::vector<std::vector<double> > excludeVector)
{
 bool keepEvent = false;
  if(includeVector.size()==0) keepEvent = true;
  for (int ti = 0; ti < numTriggers; ++ti) {
    for(uint i=0; i!=includeVector.size();++i){
      if(includeVector[i].size()==1 && triggers[ti]==includeVector[i][0]) keepEvent=true;
      else if(includeVector[i].size()==2 && (triggers[ti]>=includeVector[i][0] && triggers[ti]<=includeVector[i][1])) keepEvent=true;
    }
  }
  if(!keepEvent)
    return false;

  keepEvent = true;
  for (int ti = 0; ti < numTriggers; ++ti) {
    for(uint i=0; i!=excludeVector.size();++i){
      if(excludeVector[i].size()==1 && triggers[ti]==excludeVector[i][0]) keepEvent=false;
      else if(excludeVector[i].size()==2 && (triggers[ti]>=excludeVector[i][0] && triggers[ti]<=excludeVector[i][1])) keepEvent=false;
    }
  }

  return keepEvent;
}

// ---------------------------------------------------------------------------------------
// ------- Function to decide to include/exclude event based on the vectors passed -------
bool includeEvent(double eventParameter,
    std::vector<std::vector<double> > includeVector,
    std::vector<std::vector<double> > excludeVector)
{
  bool keepEvent = false;
  if(includeVector.size()==0) keepEvent = true;
  for(uint i=0; i!=includeVector.size();++i){
    if(includeVector[i].size()==1 && eventParameter==includeVector[i][0])
      keepEvent=true;
    else if(includeVector[i].size()==2 && (eventParameter>=includeVector[i][0] && eventParameter<=includeVector[i][1]))
      keepEvent=true;
  }
  if(!keepEvent) // if it's not in our include list, skip it
    return false;

  keepEvent = true;
  for(uint i=0; i!=excludeVector.size();++i){
    if(excludeVector[i].size()==1 && eventParameter==excludeVector[i][0])
      keepEvent=false;
    else if(excludeVector[i].size()==2 && (eventParameter>=excludeVector[i][0] && eventParameter<=excludeVector[i][1]))
      keepEvent=false;
  }

  return keepEvent; // if someone includes and excludes, exseedion will overrule

}



// ---------------------------------------------------------------------------------------
// ------------------ Function to split arg input lists by comma -------------------------
std::vector<std::string> split(std::string msg, std::string separator)
{
  boost::char_separator<char> sep(separator.c_str());
  boost::tokenizer<boost::char_separator<char> > tok(msg, sep );
  std::vector<std::string> token;
  for ( boost::tokenizer<boost::char_separator<char> >::const_iterator i = tok.begin(); i != tok.end(); ++i ) {
    token.push_back(std::string(*i));
  }
  return token;
}


// ---------------------------------------------------------------------------------------
// ------------------ Function to generate include/exclude vectors -----------------------
void genIncludeExcludeVectors(std::string optionString,
			      std::vector<std::vector<double> >& includeVector,
			      std::vector<std::vector<double> >& excludeVector)
{
  std::vector<std::string> rangeStringVector;
  std::vector<double> rangeIntVector;

  if(optionString != "-1"){
    std::vector<std::string> stringVector = split(optionString,",") ;

    for (uint i=0 ; i<stringVector.size() ; i++) {
      bool exclude = false;

      if(stringVector[i].at(0)=='x'){
        exclude = true;
        stringVector[i].erase(0,1);
      }
      rangeStringVector = split(stringVector[i],"-") ;

      rangeIntVector.clear();
      for(uint j=0; j<rangeStringVector.size();j++) {
        rangeIntVector.push_back(atof(rangeStringVector[j].c_str()));
      }
      if(exclude) excludeVector.push_back(rangeIntVector);
      else includeVector.push_back(rangeIntVector);

    }
  }
}


// ---------------------------------------------------------------------------------------
// ------------------ Function to parse the command-line arguments------------------------
void parseArguments(int argc, char** argv)
{
  std::string stringGenericOption    = "--";
  std::string stringHelp             = "--help";
  std::string stringInputFileName    = "--i";
  std::string stringOutFileName      = "--o";
  std::string stringETGammaMinEB     = "--eTGammaMinEB";
  std::string strings4s9GammaMinEB   = "--s4s9GammaMinEB";
  std::string stringeTPi0MinEB       = "--eTPi0MinEB";
  std::string stringETGammaMinEE     = "--eTGammaMinEE";
  std::string strings4s9GammaMinEE   = "--s4s9GammaMinEE";
  std::string stringeTPi0MinEE       = "--eTPi0MinEE";
  std::string stringminAmpliOverSigma= "--minAOverSigma";
  std::string stringNumEvents        = "--n";
  std::string stringMinRun           = "--minRun";
  std::string stringMaxRun           = "--maxRun";
  std::string stringMinLS            = "--minLS";
  std::string stringMaxLS            = "--maxLS";
  std::string vertex                 = "--vertex";
  std::string stringTriggers         = "--trig";
  std::string stringTechTriggers     = "--techTrig";

  // if no arguments are passed, suggest help
  if (argc < 2){
    std::cerr << "\n\tERROR: specify arguments, or use --help\n" << std::endl ;
    exit (1) ;  
  }

  // loop over input options
  for (int v=1; v<argc; v++ )
  {
    //std::cout << "argv number " << v << " is: " << argv[v] << std::endl;

    if (argv[v] == stringHelp) { // help message
      std::cout << " --help : display help" << std::endl ;
      std::cout << " --o : set name of output root file name (e.g. histograms.root)" << std::endl ;
      std::cout << " --n : number of events" << std::endl ;
      std::cout << " --i <list of strings> list of input files" << std::endl ;     
      std::cout << " --eTGammaMinEB: min eT for EB gammas" << std::endl;
      std::cout << " --s4s9GammaMinEB: min EB shower shape" << std::endl;
      std::cout << " --eTPi0MinEB min eT for EB pi0 candidate" << std::endl;
      std::cout << " --eTGammaMinEE: min eT for EE gammas" << std::endl;
      std::cout << " --s4s9GammaMinEE: min EE shower shape" << std::endl;
      std::cout << " --eTPi0MinEE: min eT for EE pi0 candidate" << std::endl;
      std::cout << " --minAOverSigma: min ampli considered for time" << std::endl;
      std::cout << " --minRun: lowest run number considered" << std::endl;
      std::cout << " --maxRun: highest run number considered" << std::endl;
      std::cout << " --minLS: lowest lumi section number considered" << std::endl;
      std::cout << " --maxLS: highest lumi section number considered" << std::endl;
      std::cout << " --vertex: require vertex@IP (1), veto it (2) or either (0, or unset)" << std::endl;
      std::cout << " --trig: L1 triggers to include (exclude with x)" << std::endl;
      std::cout << " --techTrig: L1 technical triggers to include (exclude with x)" << std::endl;
      exit(1);      }


    else if (argv[v] == stringNumEvents) { // set number of events
      std::cout << "events number" << std::endl;
      numEvents_=atoi(argv[v+1]);
      v++;
    }
    else if (argv[v] == stringMaxLS) { // set last LS of interval to be considered 
      std::cout << "max LS number" << std::endl;
      maxLS_=atoi(argv[v+1]);
      v++;
    }
    else if (argv[v] == stringMinLS) { // set first LS of interval to be considered 
      std::cout << "min LS number" << std::endl;
      minLS_=atoi(argv[v+1]);
      v++;
    }
    else if (argv[v] == stringMaxRun) { // set last run of interval to be considered 
      std::cout << "max LS number" << std::endl;
      maxRun_=atoi(argv[v+1]);
      v++;
    }
    else if (argv[v] == stringMinRun) { // set first run of interval to be considered 
      std::cout << "min run number" << std::endl;
      minRun_=atoi(argv[v+1]);
      v++;
    }
    else if (argv[v] == stringETGammaMinEB) { // choose et cut for EB single cluster
      eTGammaMinEB_ = atof(argv[v+1]);
      v++;
    }
    else if (argv[v] == strings4s9GammaMinEB) { // choose cut for EB shower shape
      s4s9GammaMinEB_ = atof(argv[v+1]);
      v++;
    }
    else if (argv[v] == stringeTPi0MinEB) { // choose et cut for EB pi0 candidate
      eTPi0MinEB_ = atof(argv[v+1]);
      v++;
    }
    else if (argv[v] == stringETGammaMinEE) { // choose et cut for EE single cluster
      eTGammaMinEE_ = atof(argv[v+1]);
      v++;
    }
    else if (argv[v] == strings4s9GammaMinEE) { // choose cut for EE shower shape
      s4s9GammaMinEE_ = atof(argv[v+1]);
      v++;
    }
    else if (argv[v] == stringeTPi0MinEE) { // choose et cut for EE pi0 candidate
      eTPi0MinEE_ = atof(argv[v+1]);
      v++;
    }
    else if (argv[v] == stringOutFileName) { // set output file
      outputRootName_ = argv[v+1];
      v++;
    }
    else if (argv[v] == stringminAmpliOverSigma) { // set min amplitude considered for time measurement
     /*
      * minAmpliOverSigma_  = atof(argv[v+1]);
     */
      v++ ;
    }
    else if (argv[v] == vertex) { // collect requirement for one vertex only or not
      flagOneVertex_  = atof(argv[v+1]);
       if (flagOneVertex_!=0 && flagOneVertex_!=1 && flagOneVertex_!=2){
         std::cout << "Not a valid value for flagOneVertex_ (0,1,2). Returning." << std::endl;
	 exit (1);}
       v++;
    } 
    else if (argv[v] == stringTriggers) { // set L1 triggers to include/exclude
      genIncludeExcludeVectors(std::string(argv[v+1]),trigIncludeVector,trigExcludeVector);
      v++;
    }
    else if (argv[v] == stringTechTriggers) { // set L1 technical triggers to include/exclude
      genIncludeExcludeVectors(std::string(argv[v+1]),ttrigIncludeVector,ttrigExcludeVector);
      v++;
    }
    // handle here the case of multiple arguments for input files
    else if (argv[v] == stringInputFileName){// && v<(argc-1) ) 

      for (int u=v+1; u<argc; u++) {

        if ( 0==std::string(argv[u]).find( stringGenericOption ) ){
          if ( 0==listOfFiles_.size())  {std::cout << "no input files listed" << std::cout;}
          //else  {std::cout << "no more files listed, found: " << argv[u] << std::cout;}
          break;
        }

        else {  listOfFiles_.push_back(argv[u]);
          v++;
        }

      }// loop on arguments following --i

      continue;

    }//end 'if input files'

    else
    {std::cout << "input format unrecognized" << std::endl; exit(1);}

    }// loop over arguments input to the program
}



// ---------------------------------------------------------------------------------------
// -------------------- struct holding sed of  histograms  -------------------------------
struct HistSet{
  
  //book histogram set w/ common suffix inside the provided TFileDirectory
  //void book(edm::Service<TFileService>& td,const std::string&);
  void book(TFileDirectory subDir,const std::string&);
  
  // fill all histos of the set with the two electron candidates
  void fill(int sc1, int sc2, int cl1, int cl2);
  
  TH1 * nVertices_;
  TH1F* mass_;
  TH1F* dZvertices_;
  TH1F* Zvertices_;
  TH1F* chi2_;
  TH1F* seedTime_;
  TH1F* secondTime_;
  TH1F* clusterTime_;
  TH1F* seedTimeDiffHist_;
  TH1F* TOFcorrections_;
  TH2F* TOFcorrectionsVSdeltaEta_;
  TH2F* clusTimeDiffHistTOFVSdeltaEtaRightVertex_, *clusTimeDiffHistTOFVSdeltaEtaWrongVertex_;
  TH1F* tColl_;
  TH2F* tCollVSdeltaEtaRightVertex_, * tCollVStimeDiffHistTOF_;
  TH1F* seedTimeDiffHistTOF_;
  TH1F* secondTimeDiffHist_;
  TH1F* secondTimeDiffHistTOF_;
  TH1F* clusTimeDiffHist_;
  TH1F* clusTimeDiffHistTOF_, *clusTimeDiffHistTOFwrongVertex_;
  TH1F* numCryBC1, *numCryBC2;
  TH2F* timeVsEtaLead_, *timeVsEtaSub_, *timeVsEta_, *outliersVsEtaPhi_,*timeVsPhi_,*timeVsSeedE_; 
  TH1F* seedAmpli_;
  TH1F* secondAmpli_;
  TH1F* diffSeedOther_, *diffSeedOtherOverErr_;
  TH1F* diffSeedSecond_, *diffSeedSecondOverErr_;
  TH2F* seedVSSecond_;  
  TH2F* e2e9VsSx_;
  TH1F* e2e9_;
  TH1F* Sx_;
  TH1F* gammaEt1_, *gammaEt2_;
  TH1F* NZmassseedTime_;
  TH1F* NZmass_;
  TH1F* ZmassBkgseedTime_, *ZmassBkg_;
  TH1F* ZmassChi2CutseedTime_, *ZmassChi2Cut_;
  TH1F* NZmassLess_;
  TH1F* NZmassMore_;
  TH1F* ZMassMoreSeedTime_;
  TH1F* ZMassLessSeedTime_;
  TH1F* numSclusters_;
  TH1F*chi2POOT_, *chi2MOOT_, *chi2IT_;
  TH1F*chi2Sig_, *chi2BkgL_, *chi2BkgH_;
  TH2F*chi2VsSeedTime_;
  TH2F*seedTime1Vsdi_eleMass_,*seedTime2Vsdi_eleMass_, *seedTime1VsseedTime2_;
  TH1F*di_eleMassSeedTimeMore_,*di_eleMassSeedTimeLess_, *di_eleMassSeedTimeIn_;
} theHists;


void HistSet::book(TFileDirectory subDir, const std::string& post) {

  nVertices_=subDir.make<TH1F>("num vertices","num vertices; num vertices",41,-0.5,40.5);
  mass_         =(TH1F*) subDir.make<TH1F>("mass","mass; m(ele,ele) [GeV]",80,50,130);
  //mass_         =(TH1F*) subDir.make<TH1F>("mass","mass; m(#gamma,#gamma) [GeV]",80,50,130);
  dZvertices_   =(TH1F*) subDir.make<TH1F>("dZvertices","dZvertices; #DeltaZ(ele_{1},ele_{2}) [cm]",250,0,25);
  Zvertices_    =(TH1F*) subDir.make<TH1F>("Zvertices","Zvertices; z vertex [cm]",250,-25,25);

  // Initialize histograms -- xtals
  chi2_                =(TH1F*) subDir.make<TH1F>("cluster chi2 ","cluster chi2 ; #chi^{2}",100,0,10);
  chi2POOT_             =(TH1F*) subDir.make<TH1F>("cluster chi2POOT ","cluster chi2POOT ; #chi^{2},(t > 3ns)",100,0,10);
  chi2MOOT_             =(TH1F*) subDir.make<TH1F>("cluster chi2MOOT ","cluster chi2MOOT ; #chi^{2} (t < -3ns)",100,0,10);
  chi2IT_               =(TH1F*) subDir.make<TH1F>("cluster chi2IT ","cluster chi2IT ; #chi^{2} (|t| < 2ns)",100,0,10);
  chi2Sig_              =(TH1F*) subDir.make<TH1F>("cluster chi2Sig ","cluster chi2Bkg ; #chi^{2} (76 < M_{ee} < 100 GeV)",100,0,10);
  chi2BkgL_             =(TH1F*) subDir.make<TH1F>("cluster chi2BkgL ","cluster chi2BkgL ; #chi^{2} (50 < M_{ee} < 76 GeV)",100,0,10);
  chi2BkgH_             =(TH1F*) subDir.make<TH1F>("cluster chi2BkgH ","cluster chi2BkgH ; #chi^{2} (100 < M_{ee} < 130 GeV)",100,0,10);
  chi2VsSeedTime_       =(TH2F*) subDir.make<TH2F>("Cluster chi2 Vs SeedTime","Cluster chi2 Vs ECAL Time; ECAL Time[ns]; #chi2",500, -25.0, 25.0, 100, 0.0, 10.0 );
  seedTime_            =(TH1F*) subDir.make<TH1F>("seed time","seed time; t_{seed} [ns]; num. seeds/0.05ns",binsTDistro_,-rangeTDistro_,rangeTDistro_);
  secondTime_          =(TH1F*) subDir.make<TH1F>("second time","second time; t_{second} [ns]; num. secs/0.05ns",binsTDistro_,-rangeTDistro_,rangeTDistro_);
  clusterTime_         =(TH1F*) subDir.make<TH1F>("cluster time","cluster time; t_{cluster} [ns]; num. clusters/0.05ns",binsTDistro_,-rangeTDistro_,rangeTDistro_);


  TOFcorrections_      = (TH1F*) subDir.make<TH1F>("TOF difference","TOF difference; #Delta TOF [ns]; num. seeds/0.05ns",binsTDistro_,-rangeTDistro_/2.,rangeTDistro_/2.);
  TOFcorrectionsVSdeltaEta_=(TH2F*) subDir.make<TH2F>("TOF corrections VS #Delta#eta","TOF corrections VS #Delta#eta; #Delta#eta; #Delta TOF [ns]; ",30,0,3.,binsTDistro_,-rangeTDistro_/2.,rangeTDistro_/2.);
  clusTimeDiffHistTOFVSdeltaEtaRightVertex_=(TH2F*) subDir.make<TH2F>("TOF-corr cluster time difference VS #Delta#eta RightVertex","TOF-corr cluster time difference VS #Delta#eta RightVertex; |#Delta#eta|; (t_{clus1} - t_{clus2}) TOF-corrected [ns]; ",30,0,3.,binsTDistro_,-rangeTDistro_/2.,rangeTDistro_/2.);
  clusTimeDiffHistTOFVSdeltaEtaWrongVertex_=(TH2F*) subDir.make<TH2F>("TOF-corr cluster time difference VS #Delta#eta WrongVertex","TOF-corr cluster time difference VS #Delta#eta WrongVertex; |#Delta#eta|;  (t_{clus1} - t_{clus2}) TOF-corrected [ns]; ",30,0,3.,binsTDistro_,-rangeTDistro_/2.,rangeTDistro_/2.);

  tCollVSdeltaEtaRightVertex_=(TH2F*) subDir.make<TH2F>("t_{coll} VS #Delta#eta RightVertex","t_{coll} VS #Delta#eta RightVertex; |#Delta#eta|; t_{coll} [ns]; ",30,0,3.,binsTDistro_,-rangeTDistro_/4.,rangeTDistro_/4.);
  tCollVStimeDiffHistTOF_=(TH2F*) subDir.make<TH2F>("TOF-corrected: (t_{clus1} + t_{clus2})/2  VS  (t_{clus1} - t_{clus2})/2","TOF-corrected: (t_{clus1} + t_{clus2})/2  VS  (t_{clus1} - t_{clus2})/2 [ns]; (t_{clus1} - t_{clus2})/2 ; (t_{clus1} + t_{clus2})/2 VS [ns]; ",binsTDistro_,-rangeTDistro_/2.,rangeTDistro_/2.,binsTDistro_,-rangeTDistro_/2.,rangeTDistro_/2.);


  seedTimeDiffHist_    =(TH1F*) subDir.make<TH1F>("time difference of seeds","seeds time difference; t_{seed1} - t_{seed2} [ns]; num. seed pairs/0.05ns",binsTDistro_,-rangeTDistro_,rangeTDistro_);
  seedTimeDiffHistTOF_ =(TH1F*) subDir.make<TH1F>("TOF-corr time difference of seeds","TOF-corr seed time difference; (t_{seed1} - t_{seed2}) TOF-corrected   [ns]; num. seed pairs/0.05ns",binsTDistro_,-rangeTDistro_,rangeTDistro_);

  secondTimeDiffHist_  =(TH1F*) subDir.make<TH1F>("time difference of seconds","second time difference; t_{second1} - t_{second2} [ns]; num. seed pairs/0.05ns",binsTDistro_,-rangeTDistro_,rangeTDistro_);//GF new
  secondTimeDiffHistTOF_ =(TH1F*) subDir.make<TH1F>("TOF-corr time difference of seconds","TOF-corr seconds time difference;  (t_{second1} - t_{second2}) TOF-corrected [ns]; num. sec pairs/0.05ns",binsTDistro_,-rangeTDistro_,rangeTDistro_);

  clusTimeDiffHist_    =(TH1F*) subDir.make<TH1F>("cluster time difference","cluster time difference;  t_{clus1} - t_{clus2} [ns]; num. cluster pairs/0.05ns",binsTDistro_,-rangeTDistro_,rangeTDistro_);
  clusTimeDiffHistTOF_ =(TH1F*) subDir.make<TH1F>("TOF-corr cluster time difference","TOF-corr cluster time difference; (t_{clus1} - t_{clus2}) TOF-corrected [ns]; num. cluster pairs/0.05ns",binsTDistro_,-rangeTDistro_,rangeTDistro_);
  //  clusTimeDiffHistTOFwrongVertex_ =(TH1F*) subDir.make<TH1F>("TOF-corr cluster time difference","TOF-corr cluster time difference; (t_{clus1} - t_{clus2}) TOF-corrected [ns]; num. cluster pairs/0.05ns",binsTDistro_,-rangeTDistro_,rangeTDistro_);
  clusTimeDiffHistTOFwrongVertex_=(TH1F*) subDir.make<TH1F>("TOF-corr cluster time difference wrong vertex","TOF-corr cluster time difference wronge vertex; (t_{clus1} - t_{clus2}) TOF-corrected [ns]; num. cluster pairs/0.05ns",binsTDistro_,-rangeTDistro_,rangeTDistro_);

  tColl_=(TH1F*) subDir.make<TH1F>("t_{coll}","t_{coll} [ns];  (t_{clus1} + t_{clus2})/2 [ns]",binsTDistro_,-rangeTDistro_/2.,rangeTDistro_/2.);


  numCryBC1            =(TH1F*) subDir.make<TH1F>("num cry in bc1","num cry in bc1; num cry",25,0.0,25.0);
  numCryBC2            =(TH1F*) subDir.make<TH1F>("num cry in bc2","num cry in bc2; num cry",25,0.0,25.0);
  timeVsEta_           =(TH2F*) subDir.make<TH2F>("timeVsEta","timeVsEta;#eta; SeedTime [ns]",100,-3.5,3.5,500,-25.0,25.0);
  timeVsEtaLead_       =(TH2F*) subDir.make<TH2F>("timeVsEtaLead","timeVsEtaLead;#eta_{lead}; SeedTime [ns]",50,-2.5,2.5,500,-25.0,25.0);
  timeVsEtaSub_        =(TH2F*) subDir.make<TH2F>("timeVsEtaSub","timeVsEtaSub; #eta_{sublead};SeedTime [ns]",50,-2.5,2.5,500,-25.0,25.0);
  outliersVsEtaPhi_    =(TH2F*) subDir.make<TH2F>("outliersVsEtaPhi","outliersVsEtaPhi; #eta; #phi",50,-2.5,2.5,72,-3.14,3.14);
  timeVsPhi_           =(TH2F*) subDir.make<TH2F>("timeVsPhi","timeVsPhi; #phi; SeedTime [ns]",100,-3.5,3.5,500,-25.0,25.0);
  seedAmpli_           =(TH1F*) subDir.make<TH1F>("E(seed)  ","E(seed) ; E [GeV]",130,0,130);
  timeVsSeedE_         =(TH2F*) subDir.make<TH2F>("timeVsSeedE","timeVsSeedE;seed E[GeV]; Seed Time [ns]",150,0.0,150,500,-25.0,25.0);
  secondAmpli_         =(TH1F*) subDir.make<TH1F>("E(second)  ","E(second) ; E [GeV]",130,0,130);
  diffSeedOther_       =(TH1F*) subDir.make<TH1F>("t_{seed}-t_{others}","t_{seed}-t_{others}; t_{seed}-t_{others} [ns]; num./0.05ns",binsTDistro_,-rangeTDistro_,rangeTDistro_);
  diffSeedOtherOverErr_ =(TH1F*) subDir.make<TH1F>("(t_{seed}-t_{others})/#sigma","(t_{seed}-t_{others})/#sigma; (t_{seed}-t_{others})/#sigma; num./0.05ns",binsTDistro_,-rangeTDistro_,rangeTDistro_);

  diffSeedSecond_      =(TH1F*) subDir.make<TH1F>("t_{seed}-t_{second}","t_{seed}-t_{second}; t_{seed}-t_{second} [ns]; num./0.05ns",binsTDistro_,-rangeTDistro_,rangeTDistro_);
  diffSeedSecondOverErr_ =(TH1F*) subDir.make<TH1F>("(t_{seed}-t_{second})/#sigma","(t_{seed}-t_{second})/#sigma; (t_{seed}-t_{second})/#sigma; num./0.05ns",binsTDistro_,-rangeTDistro_,rangeTDistro_);
  seedVSSecond_        =(TH2F*) subDir.make<TH2F>("t_{seed} VS t_{second}","t_{seed} VS t_{second}; t_{seed} [ns]; t_{second} [ns]",75,-1.5,1.5,75,-1.5,1.5);
  e2e9VsSx_        =(TH2F*) subDir.make<TH2F>("E2/E9 VS 1-E4/E1","E2/E9 VS 1-E4/E1; E2/E9; 1-E4/E1",50,0.0,5.0,50,-2.5,2.5);
  e2e9_            =(TH1F*) subDir.make<TH1F>("E2/E9 ","E2/E9; E2/E9",100,0.0,10.0);
  Sx_              =(TH1F*) subDir.make<TH1F>("1-E4/E1 ","1-E4/E2; 1-E4/E1",200,-10.0,10.0);
  gammaEt1_        =(TH1F*) subDir.make<TH1F>("Et_{1}^{e} ","Et_{1}^{e}; E_t_{1}^{e}",150,0.0,150.0);
  gammaEt2_        =(TH1F*) subDir.make<TH1F>("Et_{2}^{e} ","Et_{2}^{e}; E_t_{2}^{e}",150,0.0,150.0);
  NZmassseedTime_            =(TH1F*) subDir.make<TH1F>("NZMass seed time","NZMass seed time; t_{seed} [ns]; num. seeds/0.05ns",binsTDistro_,-rangeTDistro_,rangeTDistro_);
  ZmassBkgseedTime_            =(TH1F*) subDir.make<TH1F>("ZmassBkgseedTime","ZmassBkgseedTime; t_{seed} [ns]; num. seeds/0.05ns",binsTDistro_,-rangeTDistro_,rangeTDistro_);
  ZmassBkg_         =(TH1F*) subDir.make<TH1F>("ZmassBkg","ZmassBkg; m(e,e) [GeV]",80,50,130);
  ZmassChi2CutseedTime_            =(TH1F*) subDir.make<TH1F>("ZmassChi2CutseedTime","ZmassChi2CutseedTime; t_{seed} [ns]; num. seeds/0.05ns",binsTDistro_,-rangeTDistro_,rangeTDistro_);
  ZmassChi2Cut_         =(TH1F*) subDir.make<TH1F>("ZmassChi2Cut","ZmassChi2Cut; m(e,e) [GeV]",80,50,130);
  NZmass_         =(TH1F*) subDir.make<TH1F>("NZmass","NZmass; m(e,e) [GeV]",80,50,130);
  NZmassLess_         =(TH1F*) subDir.make<TH1F>("NZmassLess","Zmass with [t_ele1, t_ele2 ] < -3ns ; m(#gamma,#gamma) [GeV]",80,50,130);
  NZmassMore_         =(TH1F*) subDir.make<TH1F>("NZmassMore","Zmass With [t_ele1, t_ele2 ] >  3ns; m(#gamma,#gamma) [GeV]",80,50,130);
  
  ZMassLessSeedTime_            =(TH1F*) subDir.make<TH1F>("ZMassLessSeedTime","ZMassMoreSeedTime; t_{seed} [ns]; num. seeds/0.05ns",40.0,-20.0,0.0);
  ZMassMoreSeedTime_            =(TH1F*) subDir.make<TH1F>("ZMassMoreSeedTime","ZMassMoreSeedTime; t_{seed} [ns]; num. seeds/0.05ns",40.0,0.0,20.0);
  numSclusters_       =(TH1F*) subDir.make<TH1F>("num Of Photons"," Number of Photons; Number of Photons",50,0,50);

  di_eleMassSeedTimeIn_         =(TH1F*) subDir.make<TH1F>("di_eleMassSeedTimeIn_","di_eleMass[ |t_{e_{1,2}}| < 2ns ] ; m(e_{1},e_{2}) [GeV]",80,50.,130.);
  di_eleMassSeedTimeMore_         =(TH1F*) subDir.make<TH1F>("di_eleMassSeedTimeMore_","di_eleMass[ 3.0ns < t_{e_{1}} < 20.0ns or 3.0ns < t_{e_{2}} < 20ns] ; m(e_{1},e_{2}) [GeV]",80,50.,130.);
  di_eleMassSeedTimeLess_         =(TH1F*) subDir.make<TH1F>("di_eleMassSeedTimeLess_","di_eleMass[ -20.0 ns < t_{e_{1}} < -3ns or -20.0ns < t_{e_{2}} < -3ns]  ; m(e_{1},e_{2}) [GeV]",80,50.,130.);

  seedTime1Vsdi_eleMass_        =(TH2F*) subDir.make<TH2F>("seedTime1Vsdi_eleMass_","t_{e_{1}} Vs di_eleMass; t_{e_{1}} [ns];  m(e_{1},e_{2}) [GeV]",40,-10.,10.,80,50.0,130.0);
  seedTime2Vsdi_eleMass_        =(TH2F*) subDir.make<TH2F>("seedTime2Vsdi_eleMass_","t_{e_{2}} Vs di_eleMass; t_{e_{2}} [ns];  m(e_{1},e_{2}) [GeV]",40,-10.,10.,80,50.0,130.0);
  seedTime1VsseedTime2_        =(TH2F*) subDir.make<TH2F>("seedTime1Vsseedtime2_","t_{e_{1}} Vs t_{e_{2}}; t_{e_{1}} [ns];  t_{e_{2}} [ns]",40,-10.,10.,40,-10.0,10.0);

}
  
void HistSet::fill(int sc1, int sc2, int bc1, int bc2 ){

  float et1 = treeVars_.superClusterRawEnergy[sc1]/cosh( treeVars_.superClusterEta[sc1] );
  math::PtEtaPhiELorentzVectorD  el1(et1  ,
				     treeVars_.superClusterEta[sc1],
				     treeVars_.superClusterPhi[sc1],
				     treeVars_.superClusterRawEnergy[sc1] );  
  float et2 = treeVars_.superClusterRawEnergy[sc2]/cosh( treeVars_.superClusterEta[sc2] );
  math::PtEtaPhiELorentzVectorD  el2(et2 ,
				     treeVars_.superClusterEta[sc2],
				     treeVars_.superClusterPhi[sc2],
				     treeVars_.superClusterRawEnergy[sc2] );
  math::PtEtaPhiELorentzVectorD diEle = el1;
  diEle += el2;
  
  // ////////////////////////
  mass_      ->Fill(diEle.M());
  gammaEt1_->Fill(et1); gammaEt2_->Fill(et2);
  float dvertex = pow(treeVars_.superClusterVertexZ[sc1]-treeVars_.superClusterVertexZ[sc2],2);
  //dvertex       += pow(treeVars_.superClusterVertexY[sc1]-treeVars_.superClusterVertexY[sc2],2);
  //dvertex       += pow(treeVars_.superClusterVertexX[sc1]-treeVars_.superClusterVertexX[sc2],2);
  dvertex       = sqrt(dvertex);
  dZvertices_->Fill(dvertex);
  Zvertices_->Fill( (treeVars_.superClusterVertexZ[sc1]+treeVars_.superClusterVertexZ[sc2])/2 );
  nVertices_->Fill(treeVars_.nVertices);

  ClusterTime bcTime1 = timeAndUncertSingleCluster(bc1,treeVars_);
  ClusterTime bcTime2 = timeAndUncertSingleCluster(bc2,treeVars_);
  
  TOFcorrections_           -> Fill(extraTravelTime(sc2,treeVars_) - extraTravelTime(sc1,treeVars_) );
  TOFcorrectionsVSdeltaEta_ -> Fill(  fabs(treeVars_.superClusterEta[sc2]-treeVars_.superClusterEta[sc1]) , extraTravelTime(sc2,treeVars_) - extraTravelTime(sc1,treeVars_)  );

  int vtxOfThisEle=-99;
  // look for the vertex which electrons  are attached to: 
  for(int u=0; u<treeVars_.nVertices; u++){
    // matching done with 1mm tolerance
    if( fabs(treeVars_.superClusterVertexZ[sc2]-treeVars_.vtxZ[u]) < 0.1) {       vtxOfThisEle=u;     }
    //std::cout << u << "\t" << treeVars_.superClusterVertexZ[sc2] << "\t" << treeVars_.vtxZ[u] << std::endl;
  }
  //std::cout << "\n\tdebugging vertices" << std::endl;
  if(vtxOfThisEle >-99){

    for(int u=0; u<treeVars_.nVertices; u++){
      if(u==vtxOfThisEle) {
	clusTimeDiffHistTOFVSdeltaEtaRightVertex_ -> Fill(  fabs(treeVars_.superClusterEta[sc2]-treeVars_.superClusterEta[sc1]) , 
							    (bcTime1.time-extraTravelTime(sc1,u,treeVars_))  - (bcTime2.time-extraTravelTime(sc2,u,treeVars_))
							    //(bcTime1.time-extraTravelTime(sc1,treeVars_))  - (bcTime2.time-extraTravelTime(sc2,treeVars_))
							    );

	tCollVSdeltaEtaRightVertex_-> Fill(  fabs(treeVars_.superClusterEta[sc2]-treeVars_.superClusterEta[sc1]) , 
					     ( (bcTime1.time-extraTravelTime(sc1,u,treeVars_)) + (bcTime2.time-extraTravelTime(sc2,u,treeVars_)) ) /2.
					     );
	tCollVStimeDiffHistTOF_    -> Fill(  ( (bcTime1.time-extraTravelTime(sc1,u,treeVars_)) - (bcTime2.time-extraTravelTime(sc2,u,treeVars_)) ) /2. ,
					     ( (bcTime1.time-extraTravelTime(sc1,u,treeVars_)) + (bcTime2.time-extraTravelTime(sc2,u,treeVars_)) ) /2.
					     );
	tColl_                     -> Fill( ( (bcTime1.time-extraTravelTime(sc1,u,treeVars_)) + (bcTime2.time-extraTravelTime(sc2,u,treeVars_)) ) /2. );

      }// if correct vertex
      else   {
	clusTimeDiffHistTOFVSdeltaEtaWrongVertex_ -> Fill(  fabs(treeVars_.superClusterEta[sc2]-treeVars_.superClusterEta[sc1]) , 
							    (bcTime1.time-extraTravelTime(sc1,u,treeVars_))  - (bcTime2.time-extraTravelTime(sc2,u,treeVars_))
							    //(bcTime1.time-extraTravelTime(sc1,treeVars_))  - (bcTime2.time-extraTravelTime(sc2,treeVars_))
							    );
	clusTimeDiffHistTOFwrongVertex_           -> Fill( (bcTime1.time-extraTravelTime(sc1,u,treeVars_))  - (bcTime2.time-extraTravelTime(sc2,u,treeVars_) ));
      }// if wrong vertex
    }//loop on vertices

  } // if vertex matching succeeded
  //  else std::cout << "vertex was not found which matches electrons track... " << std::endl;

   /** Select events based on chi2 Cut **/
  //if(bcTime1.chi2 > 4.0 && bcTime2.chi2 > 4.0 ) continue; 
  
  chi2_->Fill(bcTime1.chi2);	  chi2_->Fill(bcTime2.chi2);
  chi2VsSeedTime_->Fill(bcTime1.seedtime, bcTime1.chi2); chi2VsSeedTime_ ->Fill(bcTime2.seedtime, bcTime2.chi2);
  // take care of the seeds
  // Fill  Tail mass Plots
  // if( fabs(diEle.M()  - 91) < 5 )
  seedTime1VsseedTime2_  ->Fill(bcTime1.seedtime, bcTime2.seedtime);
  // chi2 cuts
  if(bcTime1.chi2 < 4.0 && bcTime2.chi2 < 4.0)
     {
      ZmassChi2Cut_->Fill(diEle.M());
      ZmassChi2CutseedTime_->Fill(bcTime1.seedtime); ZmassChi2CutseedTime_->Fill(bcTime2.seedtime); 
       }
 
   // Both eles are in time |t| < 2ns
   if( fabs(bcTime1.seedtime) < 2.0 && fabs(bcTime2.seedtime) < 2.0)
     {
      di_eleMassSeedTimeIn_->Fill(diEle.M());
      }

  // either one ele is out of time t > 3ns
  if( (bcTime1.seedtime > 3.0 && bcTime1.seedtime < 20.0) ||  (bcTime2.seedtime > 3.0 && bcTime2.seedtime < 20.0) )
     {
      di_eleMassSeedTimeMore_->Fill(diEle.M());
     }
  // either one ele is out of time t < -3ns
  if( ( bcTime1.seedtime > -20.0 && bcTime1.seedtime < -3.0) ||  ( bcTime2.seedtime > -20.0 && bcTime2.seedtime < -3.0) )
     {
      di_eleMassSeedTimeLess_->Fill(diEle.M()); 
     }
  // Fill electron Time Vs di_eleMass
     seedTime1Vsdi_eleMass_ ->Fill(bcTime1.seedtime, diEle.M());
     seedTime2Vsdi_eleMass_ ->Fill(bcTime2.seedtime, diEle.M());

   // Fill t > 3, t < -3, |t| < 2  chi2 distribution/ZMass
  if(bcTime1.seedtime > 3.0 && bcTime2.seedtime > 3.0)
     {
         chi2POOT_              ->Fill(bcTime1.chi2);	  chi2POOT_->Fill(bcTime2.chi2);
         } 
  
  if(bcTime1.seedtime < -3.0 && bcTime2.seedtime < -3.0)
     {
         chi2MOOT_              ->Fill(bcTime1.chi2);	  chi2MOOT_->Fill(bcTime2.chi2);
         } 

  if(fabs(bcTime1.seedtime) < 2.0 && fabs(bcTime2.seedtime) < 2.0)
     {
         chi2IT_              ->Fill(bcTime1.chi2);	  chi2IT_->Fill(bcTime2.chi2);
         } 
   // signal==>
  if( diEle.M() > 80.0 && diEle.M() < 100.0 )
     {    
       chi2Sig_              ->Fill(bcTime1.chi2);	  chi2Sig_->Fill(bcTime2.chi2);
       NZmass_->Fill(diEle.M());
       NZmassseedTime_ ->Fill( bcTime1.seedtime) ; NZmassseedTime_->Fill( bcTime2.seedtime);
       } 

   // Backgraound ELectrons
  if( diEle.M()  < 76.0 || diEle.M() > 100.0 )
     {
       
          ZmassBkg_->Fill(diEle.M()); 
          ZmassBkgseedTime_ ->Fill( bcTime1.seedtime) ;  ZmassBkgseedTime_->Fill( bcTime2.seedtime);
       if( diEle.M() > 50.0 && diEle.M() < 76.0) 
         {
	   chi2BkgL_              ->Fill(bcTime1.chi2);	  chi2BkgL_ ->Fill(bcTime2.chi2);
         }

        if( diEle.M() > 100.0 && diEle.M() < 130.0 )
	   {
            chi2BkgH_              ->Fill(bcTime1.chi2);	  chi2BkgH_ ->Fill(bcTime2.chi2);
	    }

     }
     // Background==>
  /*
   if( diEle.M() < 80.0 ||  diEle.M() > 100.0 )
     { 
       NZmass_->Fill(diEle.M()); 
       NZmassseedTime_ ->Fill( bcTime1.seedtime) ; 
       NZmassseedTime_->Fill( bcTime2.seedtime);
       }
  */
  seedTime_            -> Fill(bcTime1.seedtime);  seedTime_->Fill(bcTime2.seedtime); 
  seedTimeDiffHist_    -> Fill( bcTime1.seedtime - bcTime2.seedtime );
  seedTimeDiffHistTOF_ -> Fill( (bcTime1.seedtime-extraTravelTime(sc1,treeVars_))  - (bcTime2.seedtime-extraTravelTime(sc2,treeVars_))  );
  
  // take care of the second-highest amplitude crystal
  if(bcTime1.second>-1) secondAmpli_->Fill(treeVars_.xtalInBCEnergy[bc1][bcTime1.second]);  // check that there's crystals beyond seed
  else                  secondAmpli_->Fill(0);  
  if(bcTime2.second>-1) secondAmpli_->Fill(treeVars_.xtalInBCEnergy[bc2][bcTime2.second]);
  else                  secondAmpli_->Fill(0);  
  if(bcTime1.second>-1 && bcTime2.second>-1)
    {
       secondTime_            -> Fill( bcTime1.secondtime);  secondTime_->Fill(bcTime2.secondtime); 
       secondTimeDiffHist_    -> Fill( bcTime1.secondtime - bcTime2.secondtime );
       secondTimeDiffHistTOF_ -> Fill( (bcTime1.secondtime-extraTravelTime(sc1,treeVars_))  - (bcTime2.secondtime-extraTravelTime(sc2,treeVars_))  );
     }
  
  
  clusterTime_         -> Fill(bcTime1.time);              clusterTime_ ->Fill(bcTime2.time);
  clusTimeDiffHist_    -> Fill(bcTime1.time - bcTime2.time );
  clusTimeDiffHistTOF_ -> Fill( (bcTime1.time - extraTravelTime(sc1,treeVars_) ) - (bcTime2.time -extraTravelTime(sc2,treeVars_)) );

  seedAmpli_->Fill(treeVars_.xtalInBCEnergy[bc1][bcTime1.seed]);
  seedAmpli_->Fill(treeVars_.xtalInBCEnergy[bc2][bcTime2.seed]);

  // Shape Study Variable
  e2e9VsSx_->Fill( treeVars_.clusterE2[sc1]/treeVars_.clusterE3x3[sc1], 1-treeVars_.clusterE2x2[sc1]/treeVars_.clusterE1[sc1] );
  e2e9VsSx_->Fill( treeVars_.clusterE2[sc2]/treeVars_.clusterE3x3[sc2], 1-treeVars_.clusterE2x2[sc2]/treeVars_.clusterE1[sc2] );
  e2e9_->Fill( treeVars_.clusterE2[sc1]/treeVars_.clusterE3x3[sc1] ); e2e9_->Fill( treeVars_.clusterE2[sc2]/treeVars_.clusterE3x3[sc2] );
  Sx_->Fill( 1-treeVars_.clusterE2x2[sc1]/treeVars_.clusterE1[sc1] ); Sx_->Fill( 1-treeVars_.clusterE2x2[sc2]/treeVars_.clusterE1[sc2] ); 
  timeVsSeedE_->Fill(treeVars_.xtalInBCEnergy[bc1][bcTime1.seed], bcTime1.seedtime); timeVsSeedE_->Fill(treeVars_.xtalInBCEnergy[bc2][bcTime2.seed], bcTime2.seedtime); timeVsEta_ -> Fill( treeVars_.superClusterEta[sc1] , (bcTime1.seedtime )); timeVsEta_ -> Fill( treeVars_.superClusterEta[sc2] , (bcTime2.seedtime )); //- extraTravelTime(sc1,treeVars_) ) );
  timeVsPhi_ -> Fill( treeVars_.superClusterPhi[sc1] , (bcTime1.seedtime )); timeVsPhi_ -> Fill( treeVars_.superClusterPhi[sc2] , (bcTime2.seedtime )); // - extraTravelTime(sc1,treeVars_) ) );
    numSclusters_->Fill(treeVars_.nSuperClusters);

    
   // SpitOut these delayed events:
  if( (bcTime1.seedtime > 3.0)  &&  (bcTime2.seedtime > 3.0) ) { NZmassMore_->Fill(diEle.M());    std::cout <<"Delayed Photon has [ SclusEt1, SclusEt2 ]=: " << "[" << et1 << " , " << et2 << "]" << " Z Mass=: " << diEle.M() << " GeV" <<" From Event Id= " << treeVars_.eventId << "..From Run Number =: " << treeVars_.runId << " And  SClus1 has [ Eta1, Eta2 ] : = "  << "[ "<< treeVars_.superClusterEta[sc1] << " , " <<  treeVars_.superClusterEta[sc2] <<  " ]"   << " And  SClus2 has [ Phi1, Phi2 ] : = " << "[ " <<  treeVars_.superClusterPhi[sc1]  << " ," <<  treeVars_.superClusterPhi[sc2] <<  " ]" << "with times [ timeEle1, timeEle2 ]: " << " [ " << bcTime1.seedtime << ", " << bcTime2.seedtime << "] " << std::endl; }  
  // std::cout << "otherstime:  " << bcTime1.otherstime << "\t" << bcTime2.otherstime << std::endl;
  // std::cout << "seedtime:  " << bcTime1.seedtime << "\t" << bcTime2.seedtime << std::endl;
 
  // SpitOut these delayed events:
  if( (bcTime1.seedtime < -3.0) && (bcTime2.seedtime < -3.0) ) {  NZmassLess_->Fill(diEle.M());  std::cout <<"Delayed Photon has [ SclusEt1, SclusEt2 ]=: " << "[" << et1 << " , " << et2 << "]" << " Z Mass=: " << diEle.M() << " GeV" <<" From Event Id= " << treeVars_.eventId << "..From Run Number =: " << treeVars_.runId << "And  SClus1 has [Eta1, Eta2 ] : = " << "[ "  << treeVars_.superClusterEta[sc1]  << " ," << treeVars_.superClusterEta[sc2]  <<  " ]"   << " And  SClus2 has [Phi1, Phi2 ] : = " << "[ " << treeVars_.superClusterPhi[sc1]  << " ," <<  treeVars_.superClusterPhi[sc2] <<  " ]" << "with times [ timeEle1, timeEle2 ]: " <<" [ " << bcTime1.seedtime << ", " << bcTime2.seedtime << "] " << std::endl; }  

if((diEle.M() > 76.0) && (diEle.M() < 100.0))
    {
    //Real Signal Positive Time
    if( (bcTime1.seedtime > 3.0 && bcTime1.seedtime < 13.0 )  ||  ( bcTime2.seedtime > 3.0  &&  bcTime2.seedtime < 13.0 ) )
       { ZMassMoreSeedTime_->Fill(bcTime1.seedtime); ZMassMoreSeedTime_->Fill(bcTime2.seedtime);
          std::cout<<"Got Nobel!" << std::endl;
	  }
    //Real Signal Negative Time
    if(( bcTime1.seedtime > -10.0 &&  bcTime1.seedtime < -3.0 )  ||  ( bcTime2.seedtime > -10.0 && bcTime2.seedtime < -3.0 )) 
       {
         ZMassLessSeedTime_->Fill(bcTime1.seedtime); ZMassLessSeedTime_->Fill(bcTime2.seedtime);
	  std::cout<<"Got Nobel!" << std::endl;
	   }
 
         }


 if(bcTime1.otherstime>-999) // check that there's crystals beyond seed
    {
      diffSeedOther_ -> Fill(bcTime1.seedtime-bcTime1.otherstime); 
      diffSeedOtherOverErr_->Fill( (bcTime1.seedtime-bcTime1.otherstime) / sqrt( pow(treeVars_.xtalInBCTimeErr[bc1][bcTime1.seed],2) -0.6*0.6+timingResParamConstEB*timingResParamConstEB + pow(bcTime1.otherstimeErr,2)) );
      diffSeedSecond_ -> Fill(bcTime1.seedtime-treeVars_.xtalInBCTime[bc1][bcTime1.second]); 
      seedVSSecond_ -> Fill(treeVars_.xtalInBCTime[bc1][bcTime1.second],bcTime1.seedtime); 
      diffSeedSecondOverErr_    -> Fill( (bcTime1.seedtime-treeVars_.xtalInBCTime[bc1][bcTime1.second]) 
					   / sqrt( pow(treeVars_.xtalInBCTimeErr[bc1][bcTime1.seed],2) 
						   +  pow(treeVars_.xtalInBCTimeErr[bc1][bcTime1.second],2)
						   - 2* 0.6*0.6 + 2*timingResParamConstEB*timingResParamConstEB 
						   )   
					   ); 
    }
  if(bcTime2.otherstime>-999) // check that there's crystals beyond seed
    {
      diffSeedOther_           -> Fill(bcTime2.seedtime-bcTime2.otherstime);
      diffSeedOtherOverErr_    ->Fill( (bcTime2.seedtime-bcTime2.otherstime) / sqrt( pow(treeVars_.xtalInBCTime[bc2][bcTime2.seed],2) -0.6*0.6+timingResParamConstEB*timingResParamConstEB + pow(bcTime2.otherstimeErr,2)) ); 
      diffSeedSecond_          -> Fill(bcTime2.seedtime-treeVars_.xtalInBCTime[bc2][bcTime2.second]); 
      diffSeedSecondOverErr_    -> Fill( (bcTime2.seedtime-treeVars_.xtalInBCTime[bc2][bcTime2.second]) 
					   / sqrt( pow(treeVars_.xtalInBCTimeErr[bc2][bcTime2.seed],2) 
						   +  pow(treeVars_.xtalInBCTimeErr[bc2][bcTime2.second],2)
						   - 2* 0.6*0.6 + 2*timingResParamConstEB*timingResParamConstEB 
						   )   
					   ); 
    }
  
  numCryBC1->Fill(bcTime1.numCry);
  numCryBC2->Fill(bcTime2.numCry);
  
  timeVsEtaLead_ -> Fill( treeVars_.superClusterEta[sc1] , (bcTime1.time - extraTravelTime(sc1,treeVars_) ) );
  timeVsEtaSub_  -> Fill( treeVars_.superClusterEta[sc2] , (bcTime2.time - extraTravelTime(sc2,treeVars_) ) ); 
  // catch location of time outliers
  if ( fabs(    (bcTime1.seedtime - extraTravelTime(sc1,treeVars_))    )>1.5 )  outliersVsEtaPhi_ -> Fill( treeVars_.superClusterEta[sc1] , treeVars_.superClusterPhi[sc1]); 
  if ( fabs(    (bcTime2.seedtime - extraTravelTime(sc2,treeVars_))    )>1.5 )  outliersVsEtaPhi_ -> Fill( treeVars_.superClusterEta[sc2] , treeVars_.superClusterPhi[sc2]); 

}
// end HistSet::fill


// ---------------------------------------------------------------------------------------
// ------------------ Function to initialize the histograms ------------------------------
void initializeHists(TFileDirectory subDir){

  mass_         = subDir.make<TH1F>("mass global","mass (global); m(ele,ele) [GeV]",80,50,130);
  dZvertices_   = subDir.make<TH1F>("dZvertices global","dZvertices (global); #DeltaZ(ele_{1},ele_{2}) [cm]",250,0,25);
  Zvertices_    = subDir.make<TH1F>("Zvertices global","Zvertices (global); z vertex [cm]",250,-25,25);
  nVertices_=subDir.make<TH1F>("num vertices global","num vertices (global); num vertices",41,-0.5,40.5);
  gammaEt1_ = subDir.make<TH1F>("Et1 global","Et1 (global); Et_{1}^{#gamma} [GeV]",150,0,150);
  gammaEt2_ = subDir.make<TH1F>("Et2 global","Et2 (global); Et_{2}^{#gamma} [GeV]",150,0,150);
  numSclusters      =(TH1F*) subDir.make<TH1F>("num Of Photons"," Number of Photons; Number of Photons",50,0,50);
}//end initializeHists


// ---------------------------------------------------------------------------------------
//! main program
int main (int argc, char** argv)
{
  // First parse arguments
  parseArguments(argc, argv);

  if (listOfFiles_.size()==0){
    std::cout << "\tno input file found" << std::endl;
    return(1);
  }
  else{
    std::cout << "\tfound " << listOfFiles_.size() << " input files: " << std::endl;
    for(std::vector<std::string>::const_iterator  file_itr=listOfFiles_.begin(); file_itr!=listOfFiles_.end(); file_itr++){
      std::cout << "\t" << (*file_itr) << std::endl;
    }
  }

  // Tree construction
  // FIX should turn this string into a configurable 
  TChain * chain = new TChain ("EcalTimeAnalysis") ;  // ntuple producer in CMSSW CVS
  //TChain * chain = new TChain ("EcalTimePi0Analysis") ;  // ntuple producer in UserCode/UMN space
  std::vector<std::string>::const_iterator file_itr;
  for(file_itr=listOfFiles_.begin(); file_itr!=listOfFiles_.end(); file_itr++){
    chain->Add( (*file_itr).c_str() );
  }
  int nEntries = chain->GetEntries () ;
  if (numEvents_==-1) numEvents_ = nEntries;
  std::cout << "\n\tFOUND "         <<  listOfFiles_.size() << " input files" << std::endl ;    
  std::cout << "\n\tFOUND "         <<  nEntries << " events" << std::endl ;    
  std::cout << "\tWILL run on: "    <<  numEvents_ << " events" << std::endl;
  std::cout << "\tOutput file: "    <<  outputRootName_ << std::endl;
  std::cout << "\tminAOverSigma: "  <<  minAmpliOverSigma_ << std::endl;
  std::cout << "\teTGammaMinEB: "   <<  eTGammaMinEB_ << std::endl;
  std::cout << "\ts4s9GammaMinEB: " <<  s4s9GammaMinEB_ << std::endl;
  std::cout << "\teTPi0MinEB: "     <<  eTPi0MinEB_ << std::endl;
  std::cout << "\teTGammaMinEE: "   <<  eTGammaMinEE_ << std::endl;
  std::cout << "\ts4s9GammaMinEE: " <<  s4s9GammaMinEE_ << std::endl;
  std::cout << "\teTPi0MinEE: "     <<  eTPi0MinEE_ << std::endl;
  std::cout << "\tminRun: "         <<  minRun_ << std::endl;
  std::cout << "\tmaxRun: "         <<  maxRun_ << std::endl;
  std::cout << "\tminLS: "          <<  minLS_ << std::endl;
  std::cout << "\tmaxLS: "          <<  maxLS_ << std::endl;
	
  setBranchAddresses (chain, treeVars_);
  
  // setting up the TFileService in the ServiceRegistry;
  edmplugin::PluginManager::Config config;
  edmplugin::PluginManager::configure(edmplugin::standard::config());
  std::vector<edm::ParameterSet> psets;
  edm::ParameterSet pSet;
  pSet.addParameter("@service_type",std::string("TFileService"));
  pSet.addParameter("fileName",std::string("TimePerf-plots.root")); // this is the file TFileService will write into
  psets.push_back(pSet);
  static edm::ServiceToken services(edm::ServiceRegistry::createSet(psets));
  static edm::ServiceRegistry::Operate operate(services);
  edm::Service<TFileService> fs;

  TFileDirectory subDirECALECAL=fs->mkdir("ECALECAL");  
  HistSet plotsECALECAL;
  plotsECALECAL.book(subDirECALECAL,std::string("ECALECAL"));

  TFileDirectory subDirEBEB=fs->mkdir("EBEB");  
  HistSet plotsEBEB;
  plotsEBEB.book(subDirEBEB,std::string("EBEB"));
  
  TFileDirectory subDirEEEE=fs->mkdir("EEEE");  
  HistSet plotsEEEE;
  plotsEEEE.book(subDirEEEE,std::string("EEEE"));
    
  TFileDirectory subDirEBEE=fs->mkdir("EBEE");  
  HistSet plotsEBEE;
  plotsEBEE.book(subDirEBEE,std::string("EBEE"));
    
  TFileDirectory subDirEBEBequalShare=fs->mkdir("EBEBequalShare");  
  HistSet plotsEBEBequalShare;
  plotsEBEBequalShare.book(subDirEBEBequalShare,std::string("EBEBequalShare"));
    
  TFileDirectory subDirEBEBunevenShare=fs->mkdir("EBEBunevenShare");  
  HistSet plotsEBEBunevenShare;
  plotsEBEBunevenShare.book(subDirEBEBunevenShare,std::string("EBEBunevenShare"));
  
  timeCorrector theCorr;
  //std::cout << "\ncreated object theCorr to be used for timeVsAmpliCorrections" << std::endl;
  //std::cout << "\ninitializing theCorr" << std::endl;
  //theCorr.initEB( std::string("EBmod4") );
  //theCorr.initEE( std::string("EElow") );
  theCorr.initEB( "EB" );
  theCorr.initEE( "EE" );


  //Initialize output root file
  //saving_ = new TFile(outputRootName_.c_str (),"recreate");

  // Initialize the histograms
  TFileDirectory subDirGeneral=fs->mkdir("General");  
  initializeHists(subDirGeneral);

  int eventCounter = 0;
  /////////////////////////////////////////////////////
  // Main loop over entries
  for (int entry = 0 ; (entry < nEntries && eventCounter < numEvents_); ++entry)
  {
    chain->GetEntry (entry) ;
    // Keep the event?
    bool keepEvent = includeEvent(treeVars_.l1ActiveTriggers,
        treeVars_.l1NActiveTriggers,trigIncludeVector,trigExcludeVector)
            && includeEvent(treeVars_.l1ActiveTechTriggers,
                treeVars_.l1NActiveTechTriggers,ttrigIncludeVector,ttrigExcludeVector);
    if(!keepEvent)
      continue;

    
    // do analysis if the run is in the desired range  
    if( treeVars_.runId<minRun_  || maxRun_<treeVars_.runId) continue;
    
    // do analysis if the LS is in the desired range  
    if( treeVars_.lumiSection<minLS_  || maxLS_<treeVars_.lumiSection) continue;
    
    //bool verticesAreOnlyNextToNominalIP;
    int  count=0;
    
    for(int v=0; v<treeVars_.nVertices; v++  )
	{        if (fabs(treeVars_.vtxZ[0])<15) count++; }
    
   // if ( treeVars_.nVertices >0 && count==treeVars_.nVertices ) verticesAreOnlyNextToNominalIP = true;
   // else                                                        verticesAreOnlyNextToNominalIP = false;
    
    //    --vertex: require vertex@IP (1), veto it (2) or either (0, or unset)
   // if (flagOneVertex_ ==1 && (!verticesAreOnlyNextToNominalIP) ) continue;
   // if (flagOneVertex_ ==2 && (verticesAreOnlyNextToNominalIP) )  continue;
    
    // if evet being actually processed, increment counter of analyzed events
    eventCounter++;
    
    speak_=false;
    if (entry<10 || entry%10000==0) speak_=true;

    if (speak_)  std::cout << "\n\n------> reading entry " << entry << "\tLS: " << treeVars_.lumiSection << " <------\n" ; 
    if (speak_)  std::cout << "  found " << treeVars_.nSuperClusters << " superclusters" << std::endl ;
    if (speak_)  std::cout << "  found " << treeVars_.nClusters << " basic clusters" << std::endl ;

    // skipe evenet with 3 or more electrons
   // if ( treeVars_.nSuperClusters >= 3 ) continue; 
    ///////////////////////////////////////////////////////////////////////
    // outer loop on supercluster
    numSclusters->Fill(treeVars_.nSuperClusters);
    for (int sc1=0; sc1<treeVars_.nSuperClusters; sc1++){

      float et1 = treeVars_.superClusterRawEnergy[sc1]/cosh( treeVars_.superClusterEta[sc1] );
      //if (et1<20) continue;
      if (et1 < 20.0) continue;

      math::PtEtaPhiELorentzVectorD  el1(et1  ,
					 treeVars_.superClusterEta[sc1],
					 treeVars_.superClusterPhi[sc1],
					 treeVars_.superClusterRawEnergy[sc1] );

      ///////////////////////////////////////////////////////////////////////
      // inner loop on supercluster
      for (int sc2=(sc1+1); sc2<treeVars_.nSuperClusters; sc2++){

	float et2 = treeVars_.superClusterRawEnergy[sc2]/cosh( treeVars_.superClusterEta[sc2] );
	///if (et2<20) continue;
	if (et2 < 20.0) continue;
	
	math::PtEtaPhiELorentzVectorD  el2(et2 ,
					   treeVars_.superClusterEta[sc2],
					   treeVars_.superClusterPhi[sc2],
					   treeVars_.superClusterRawEnergy[sc2] );
      
	// there seems to be a problem with vertexing - since nearly none of the electrons have the same vertex... CHECK!
	float dvertex = pow(treeVars_.superClusterVertexZ[sc1]-treeVars_.superClusterVertexZ[sc2],2);
	//dvertex       += pow(treeVars_.superClusterVertexY[sc1]-treeVars_.superClusterVertexY[sc2],2);
	//dvertex       += pow(treeVars_.superClusterVertexX[sc1]-treeVars_.superClusterVertexX[sc2],2);
	dvertex       = sqrt(dvertex);
	
	math::PtEtaPhiELorentzVectorD diEle = el1;
	diEle += el2;

	// ////////////////////////
	mass_      ->Fill(diEle.M());
	dZvertices_->Fill(dvertex);
	Zvertices_->Fill( (treeVars_.superClusterVertexZ[sc1]-treeVars_.superClusterVertexZ[sc2])/2 );
	nVertices_->Fill(treeVars_.nVertices);
        gammaEt1_->Fill(et1);
        gammaEt2_->Fill(et2);


	// require invariant mass
	//if( fabs( diEle.M() - 91 ) > 20 ) continue;
	if( fabs( diEle.M() - 91 ) > 61 ) continue;
	// require two electrons from the same vertex
	//if ( dvertex > 0.01 )             continue; 
	//if ( dvertex > 0.1 )             continue;
	//if ( dvertex > 15 )             continue;

	if(0) std::cout << "di-electron system mass: " << diEle.M() << " vertex distance: " << dvertex << std::endl;

	// at this stage I have a suitable di-electron system for time studies

	float tmpEne=-9999;
	// loop on BC and match to sc1  ===============
	int bc1=-1;
	for (int bc=0; bc<treeVars_.nClusters; bc++){
	  if ( (pow(treeVars_.superClusterEta[sc1]-treeVars_.clusterEta[bc],2)+ pow(treeVars_.superClusterPhi[sc1]-treeVars_.clusterPhi[bc],2) ) < 0.02 
	       && treeVars_.clusterEnergy[bc]>tmpEne) {
	    tmpEne=treeVars_.clusterEnergy[bc];
	    bc1=bc;
	  }// end - if good bc candidate
	}// end - loop over BC

	
	tmpEne=-9999;
	// loop on BC and match to sc2 ==============
	int bc2=-1;
	for (int bc=0; bc<treeVars_.nClusters; bc++){
	  if ( pow(treeVars_.superClusterEta[sc2]-treeVars_.clusterEta[bc],2)+ pow(treeVars_.superClusterPhi[sc2]-treeVars_.clusterPhi[bc],2) < 0.02 
	       && treeVars_.clusterEnergy[bc]>tmpEne) {
	    tmpEne=treeVars_.clusterEnergy[bc];
	    bc2=bc;
	  }// end - if good bc candidate
	}// end - loop over BC
	
	// protect in case of no matching
	if(bc1==-1 || bc2==-1) continue;
	if(0) {
	std::cout << "\n\nsc1 : " << treeVars_.superClusterEta[sc1] << " " << treeVars_.superClusterPhi[sc1] << " " << treeVars_.superClusterRawEnergy[sc1] << std::endl;
	std::cout << "bc1 : " << treeVars_.clusterEta[bc1] << " " << treeVars_.clusterPhi[bc1] << " " << treeVars_.clusterEnergy[bc1] << "\n"<< std::endl;
	std::cout << "sc2 : " << treeVars_.superClusterEta[sc2] << " " << treeVars_.superClusterPhi[sc2] << " " << treeVars_.superClusterRawEnergy[sc2] << std::endl;
	std::cout << "bc2 : " << treeVars_.clusterEta[bc2] << " " << treeVars_.clusterPhi[bc2] << " " << treeVars_.clusterEnergy[bc2] << std::endl;
	}
	
	ClusterTime bcTime1 = timeAndUncertSingleCluster(bc1,treeVars_);
	ClusterTime bcTime2 = timeAndUncertSingleCluster(bc2,treeVars_);

	if(! (bcTime1.isvalid && bcTime2.isvalid) ) continue;

	// fill the structures which hold all the plots
	plotsECALECAL.fill(sc1,sc2, bc1,bc2);
	if      ( fabs(treeVars_.clusterEta[bc1])<1.479    &&  fabs(treeVars_.clusterEta[bc2])<1.479 ){
 	  plotsEBEB.fill(sc1,sc2, bc1,bc2);

	  float energyRatio1 = treeVars_.xtalInBCEnergy[bc1][bcTime1.seed];
	  if(bcTime1.second>-1) {energyRatio1 /= treeVars_.xtalInBCEnergy[bc1][bcTime1.second]; }
	  else { energyRatio1 /= 99999; }
	  float energyRatio2 = treeVars_.xtalInBCEnergy[bc2][bcTime2.seed];
	  if(bcTime2.second>-1) {energyRatio2 /= treeVars_.xtalInBCEnergy[bc2][bcTime2.second]; }
	  else { energyRatio2 /= 99999; }

	  float minRatio = 0.7; float maxRatio = 1.3;
	  if(minRatio<energyRatio1 && minRatio<energyRatio2 && energyRatio1<maxRatio && energyRatio2<maxRatio) 	  plotsEBEBequalShare.fill(sc1,sc2, bc1,bc2);  

	  minRatio = 2; maxRatio = 10;
	  if(minRatio<energyRatio1 && minRatio<energyRatio2 && energyRatio1<maxRatio && energyRatio2<maxRatio) 	  plotsEBEBunevenShare.fill(sc1,sc2, bc1,bc2);  
	  
	}// if EBEB, and subcases
	else if ( fabs(treeVars_.clusterEta[bc1])>1.5    &&  fabs(treeVars_.clusterEta[bc2])>1.5 ) 	  plotsEEEE.fill(sc1,sc2, bc1,bc2);
	else if ( (fabs(treeVars_.clusterEta[bc1])<1.4 && fabs(treeVars_.clusterEta[bc2])>1.5) ||
		  (fabs(treeVars_.clusterEta[bc1])>1.5 && fabs(treeVars_.clusterEta[bc2])<1.4)    ) 	plotsEBEE.fill(sc1,sc2, bc1,bc2);

	// if I've found a pair of supercluster, bail out of loop to repeat using twice the same supercluster
	break;	
	
      }// end loop sc2
    }// end loop sc1
    
  }   // end of loop over entries
  

  delete chain ;
  
  return 0 ;
}
