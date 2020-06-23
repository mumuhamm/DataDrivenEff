/// \author  - Muhammad Alibordi
// Test of RooKeyPDF has ability to discriminate the background and

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TH1D.h"
#include "TRandom.h"
//#include "/afs/cern.ch/work/m/mumuhamm/private/CMSSW_10_2_7/src/Bplus/RooJohnsonLocal.cxx"

using namespace RooFit ;
using namespace std;


void bp_ctEff_MC(string filename ){



  auto ctform = new TF1("ctEffFn","expo(0)*ROOT::Math::Chebyshev4(x,[2],[3],[4],[5],[6])",0.007,0.3);
  ctform->SetParameter(1,-1);
  ctform->FixParameter(2,1);


    double ltbin[] = {0.00655738,0.0131148,0.0196721,0.0262295,0.0327869,0.0393443,0.0459016,0.052459,0.0590164,0.0655738,0.0721311,0.0786885,0.0852459,0.0918033,0.0983607,0.104918,0.111475,0.118033,0.12459,0.131148,0.137705,0.144262,0.15082,0.157377,0.163934,0.170492,0.177049,0.183607,0.190164,0.196721,0.203279,0.209836,0.216393,0.222951,0.229508,0.236066,0.242623,0.24918,0.255738,0.262295,0.268852,0.27541,0.281967,0.288525,0.295082};//,0.301639,0.308197,0.314754,0.321311,0.327869,0.334426,0.340984,0.347541,0.354098,0.360656,0.367213,0.37377,0.380328,0.386885,0.393443,0.4};
    
  Int_t  binnum = sizeof(ltbin)/sizeof(Double_t)-1;


  auto resobulk = new TH1D("resobulk", "Resobulk ; From B_{d}^{0}; Events ",  binnum, -0.004, 0.004);
  auto hctaucut = new TH1D("hctaucut","ct eff from MC   ;ct[cm];#epsilon [a.u.]",  binnum, ltbin);
  auto hctaunocut = new TH1D("hctaunocut","Gen lifetime ;ct[cm];#epsilon [a.u.]",  binnum, ltbin);

  Float_t ctreco , npv, cterr;
  Double_t kappa_val;
  Double_t par1, par2, par3, par4, par5, par6, par7;
    Double_t time = 4.911e-02;//5.0095320e-02;//2018 pdg Bplus 491.1micron, B0 = 455.7 micron, 500.95320
  TTree          *fChain, *copy;
  TFile* f = new TFile(filename.c_str());
  fChain = new TTree;
  fChain=(TTree*)f->Get("treeFit");
  if (!fChain){
    cout << "No TTree found in input file, returning" << endl;
    return;
  }
  TRandom3 *decay = new TRandom3(456);
  fChain->SetBranchAddress("BsCt2DMC",&ctreco);
  fChain->SetBranchAddress("BsCt2DMCErr",&cterr);
  Long64_t nentries = fChain->GetEntries();
  cout << "Start Processing " << nentries << " events" <<"\n";
  for (int jentry=0; jentry<nentries;jentry++) {

    fChain->GetEntry(jentry);
    if (cterr !=cterr)continue;
    if (ctreco !=ctreco)continue;
    hctaucut->Fill(ctreco);
    
  }

  TFile *fb0 = new TFile("/Users/md/Documents/BsToJpsiPhi/BPB0_3D/fittree_ntuBdMC2017.root");
    
  TTree* b0file = (TTree*)fb0->Get("treeFit");
  Int_t n_entriesb0 = b0file->GetEntries();
  Float_t ctrecob0, ctgenb0, cterrb0;
  b0file->SetBranchAddress("BsCt2DMC",&ctrecob0);
  b0file->SetBranchAddress("BsCt2DMC_GEN",&ctgenb0);
  b0file->SetBranchAddress("BsCt2DMCErr",&cterrb0);
  for( int k = 0 ; k < nentries ; k++ )
    {
      b0file->GetEntry(k);
      resobulk->Fill(ctgenb0-ctrecob0);
      //pullbulk->Fill((ctgenb0-ctrecob0)/cterrb0);
    }
    
  for( int i = 0 ; i < nentries ; i++ )
    {
      Double_t ctaugen1 =  decay->Exp(time);
      Double_t actualctgen1 =  ctaugen1;// + cterror->GetRandom();//resobulk->GetRandom();
      hctaunocut->Fill(actualctgen1);
    }
    
   
 
  TCanvas *ceff= new TCanvas("ceff", "ceff",0,0,600,600);
  auto hdivideI = (TH1D*)hctaucut->Clone("hdivideI");
  hdivideI->Sumw2();
  hdivideI->Divide(hctaunocut);
  hdivideI->Draw("colz");    //colz ep
 /* hdivideI->Fit(ctform);
  TFitResultPtr r =    hdivideI->Fit(ctform, "S");
  TMatrixDSym cov = r->GetCovarianceMatrix();   //  to access the covariance matrix
  Double_t chi2   = r->Chi2();                  // to retrieve the fit chi2
  Double_t par0   = r->Parameter(0);            // retrieve the value for the parameter 0
  Double_t err0   = r->ParError(0);             // retrieve the error for the parameter 0
  r->Print("V");                                // print full information of fit including covariance matrix
  */


}
