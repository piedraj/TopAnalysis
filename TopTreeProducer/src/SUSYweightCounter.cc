// 
// Package:    SUSYweightCounter
// Class:      SUSYweightCounter
// 
/**\class SUSYweightCounter SUSYweightCounter.cc TopAnalysis/SUSYweightCounter/src/SUSYweightCounter.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
*/
//
// Original Author:  
//         Created:  Fri Jan  9 11:55:23 CET 2009
// $Id$
//
// 

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
//#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Utilities/interface/InputTag.h"


#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"


#include "TH1F.h"
#include "TH1D.h"
#include "TFile.h"
#include <math.h> 
#include <sstream>
#include <string>
#include <stdlib.h>
#include <string.h>
#include "TMath.h"


using namespace std;
using namespace reco;
using namespace edm;


//
// class declaration
//



class SUSYweightCounter : public edm::EDAnalyzer {
public:
  explicit SUSYweightCounter(const edm::ParameterSet&);
  ~SUSYweightCounter();
  
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
 
  string theHistosFileName;
  TFile* theHistosFile;

 TH1D * h_eventCount ;
 TH1D * h_hdampWeight;
 TH1D * h_pdfWeight;

 bool doLHE;
 bool doHdamp;
 bool doPdf;
 int  nPdf;
 bool isaMCatNLO;

};

//
// constructors and destructor
//
SUSYweightCounter::SUSYweightCounter(const edm::ParameterSet& iConfig):
  theHistosFileName(iConfig.getUntrackedParameter<string>("histosFileName")) ,
  doLHE(iConfig.getUntrackedParameter<bool>("doLHE",false)),
  doHdamp(iConfig.getUntrackedParameter<bool>("doHdamp",false)),
  doPdf(iConfig.getUntrackedParameter<bool>("doPdf",false)),
  nPdf(iConfig.getUntrackedParameter<int>("nPdf",102)),
  isaMCatNLO(iConfig.getUntrackedParameter<bool>("isaMCatNLO",true))
{

}


SUSYweightCounter::~SUSYweightCounter()
{}



void
SUSYweightCounter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
// In the future, this should be checked at AOD level in the  Handle<LHEEventProduct> comments;
//   iEvent.getByLabel("source", comments);	# model T2tt_mstop_mLSP cross-section uncertainty

                edm::Handle<LHEEventProduct> lheEventHandle;
       if(doLHE){
                    iEvent.getByLabel("externalLHEProducer",lheEventHandle);
    //
    //                //<weightgroup name='scale_variation' combine='envelope'>
    //                //<weight id='1001'> muR=1 muF=1 hdamp=mt=172.5 </weight>
    //                //<weight id='1002'> muR=1 muF=2 hdamp=mt=172.5 </weight>
    //                //<weight id='1003'> muR=1 muF=0.5 hdamp=mt=172.5 </weight>
    //                //<weight id='1004'> muR=2 muF=1 hdamp=mt=172.5 </weight>
    //                //<weight id='1005'> muR=2 muF=2 hdamp=mt=172.5 </weight>
    //                //<weight id='1006'> muR=2 muF=0.5 hdamp=mt=172.5 </weight>
    //                //<weight id='1007'> muR=0.5 muF=1 hdamp=mt=172.5 </weight>
    //                //<weight id='1008'> muR=0.5 muF=2 hdamp=mt=172.5 </weight>
    //                //<weight id='1009'> muR=0.5 muF=0.5 hdamp=mt=172.5 </weight>
    //

       for(int w=0;w<9;w++){
        const LHEEventProduct::WGT& wgt = lheEventHandle->weights().at(w);
            h_hdampWeight->Fill(w,wgt.wgt);
	}
}
//<weightgroup name='hdamp_variation' combine='envelope'>
////<weight id='3001'> muR=1 muF=1 hdamp=0 </weight>
////<weight id='3002'> muR=1 muF=2 hdamp=0 </weight>
////<weight id='3003'> muR=1 muF=0.5 hdamp=0 </weight>
////<weight id='3004'> muR=2 muF=1 hdamp=0 </weight>
////<weight id='3005'> muR=2 muF=2 hdamp=0 </weight>
////<weight id='3006'> muR=2 muF=0.5 hdamp=0 </weight>
////<weight id='3007'> muR=0.5 muF=1 hdamp=0 </weight>
////<weight id='3008'> muR=0.5 muF=2 hdamp=0 </weight>
////<weight id='3009'> muR=0.5 muF=0.5 hdamp=0 </weight>
////<weight id='3010'> muR=1 muF=1 hdamp=86.25 </weight>
////<weight id='3011'> muR=1 muF=2 hdamp=86.25 </weight>
////<weight id='3012'> muR=1 muF=0.5 hdamp=86.25 </weight>
////<weight id='3013'> muR=2 muF=1 hdamp=86.25 </weight>
////<weight id='3014'> muR=2 muF=2 hdamp=86.25 </weight>
////<weight id='3015'> muR=2 muF=0.5 hdamp=86.25 </weight>
////<weight id='3016'> muR=0.5 muF=1 hdamp=86.25 </weight>
////<weight id='3017'> muR=0.5 muF=2 hdamp=86.25 </weight>
////<weight id='3018'> muR=0.5 muF=0.5 hdamp=86.25 </weight>
////<weight id='3019'> muR=1 muF=1 hdamp=350 </weight>
////<weight id='3020'> muR=1 muF=2 hdamp=350 </weight>
////<weight id='3021'> muR=1 muF=0.5 hdamp=350 </weight>
////<weight id='3022'> muR=2 muF=1 hdamp=350 </weight>
////<weight id='3023'> muR=2 muF=2 hdamp=350 </weight>
////<weight id='3024'> muR=2 muF=0.5 hdamp=350 </weight>
////<weight id='3025'> muR=0.5 muF=1 hdamp=350 </weight>
////<weight id='3026'> muR=0.5 muF=2 hdamp=350 </weight>
////<weight id='3027'> muR=0.5 muF=0.5 hdamp=350 </weight>
//
if(doHdamp){  
    for(int w=9+nPdf;w<9+nPdf+27;w++){
        const LHEEventProduct::WGT& wgt = lheEventHandle->weights().at(w);
            h_hdampWeight->Fill(w-nPdf,wgt.wgt);
	}
}
//    <weightgroup type='PDF_variation' combine='gaussian'>
//      <weight id='2001'> pdfset=292201 </weight>
//      <weight id='2002'> pdfset=292202 </weight>
//      <weight id='2003'> pdfset=292203 </weight>
// ...
//      <weight id='2100'> pdfset=292300 </weight>
//     <weight id='2101'> pdfset=292301 </weight>
//     <weight id='2102'> pdfset=292302 </weight>/
//    </weightgroup>
/* POWHEGV2
  100 NNPDF3.0 PDF set 260001-260100
  NNPDF 3.0 alphas=0.117 variation PDF set 265000
  NNPDF 3.0 alphas=0.119 variation PDF set 266000
  52+1 CT10 PDF variations PDF set 11000-11052
  CT10 alphas=0.117 variation PDF set 11067
  CT10 alphas=0.119 variation PDF set 11069
  50+1 MMHT2014nlo68clas118 PDF variations PDF set 25200-25250
  5 MMHT2014nlo68cl 5 alphas variations PDF set 25260-25264
*/

/*FxFx 102 weight variations
     <weightgroup type='PDF_variation' combine='gaussian'>
      <weight id='2001'> pdfset=292201-292302>
	102 variations on NNPDF30_nlo_nf_5_pdfas
*/

/* MLM 446-9=437 weight variations
   <weightgroup type="NNPDF30_lo_as_0130.LHgrid" combine="hessian"> 100 variations. weight id="10" - 110
   <weightgroup type="NNPDF30_lo_as_0130_nf_4.LHgrid" combine="hessian"> 100 variations. weight id="111-211
  <weightgroup type="NNPDF30_lo_as_0118.LHgrid" combine="hessian">
    <weight id="212">Member 0</weight>
  <weightgroup type="NNPDF23_lo_as_0130_qed.LHgrid" combine="hessian"> 100 variations.  <weight id="213"-313
  <weightgroup type="NNPDF23_lo_as_0119_qed.LHgrid" combine="hessian">
    <weight id="314">Member 0</weight>
  <weightgroup type="cteq6l1.LHgrid" combine="hessian">
    <weight id="315">Member 0</weight>
  <weightgroup type="MMHT2014lo68cl.LHgrid" combine="hessian"> 50 variations.
    <weight id="316-366
  <weightgroup type="MMHT2014lo_asmzsmallrange.LHgrid" combine="hessian">
    <weight id="367">Member 0</weight>
    <weight id="368">Member 1</weight>
    <weight id="369">Member 2</weight>
  <weightgroup type="HERAPDF15LO_EIG.LHgrid" combine="hessian"> 20 variations
    <weight id="370-390"
  <weightgroup type="NNPDF30_nlo_as_0118.LHgrid" combine="hessian">
    <weight id="391">Member 0</weight>
  <weightgroup type="NNPDF23_nlo_as_0119.LHgrid" combine="hessian">
    <weight id="392">Member 0</weight>
  <weightgroup type="CT10nlo.LHgrid" combine="hessian"> 52 variations.    <weight id="393-445
  <weightgroup type="MMHT2014nlo68cl.LHgrid" combine="hessian">
    <weight id="446">Member 0</weight>
*/


if(doPdf){
    for(int w=9;w<9+nPdf;w++){
        const LHEEventProduct::WGT& wgt = lheEventHandle->weights().at(w);
            h_pdfWeight->Fill(w-8,wgt.wgt);
                }
}

if (isaMCatNLO){
const LHEEventProduct::WGT& wgt = lheEventHandle->weights().at(0);
double w= wgt.wgt;
h_eventCount->Fill(TMath::Sign((double)1.0, w));}
else h_eventCount->Fill(1);
}




// ------------ method called once each job just before starting event loop  ------------
void 
SUSYweightCounter::beginJob()
{
  
  theHistosFile = new TFile(theHistosFileName.c_str(), "RECREATE");
  theHistosFile->cd();

  h_eventCount = new TH1D("h_eventCount", ";events;", 2, -2.0, 2.0); 
  h_hdampWeight = new TH1D("h_hdampWeight","weight;sum;",36,-0.5,35.5);
  h_pdfWeight =  new TH1D("h_pdfWeight","weight;sum;",nPdf,0.5,nPdf+0.5);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SUSYweightCounter::endJob() {
  
  theHistosFile->cd();
  h_eventCount->Write();
  h_hdampWeight->Write();
  h_pdfWeight->Write();
  theHistosFile->Write();
  theHistosFile->Close();
  
}


//define this as a plug-in
DEFINE_FWK_MODULE(SUSYweightCounter);

