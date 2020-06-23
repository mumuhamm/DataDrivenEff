//Auhor - Muhammad Alibordi
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TF1.h>
#include <TH1D.h>
#include "TGraphAsymmErrors.h"
#include "TVirtualFFT.h"
#include "TBinomialEfficiencyFitter.h"
#include "TVectorF.h"
#include "TPaveText.h"
#include <vector>
#include "TTree.h"
#include "TFile.h"
#include "TEventList.h"
#include "Riostream.h"
#include "string.h"
#include "TList.h"
#include "TDirectory.h"
#include "TCut.h"
#include "TChain.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TEfficiency.h"
#include "TLorentzVector.h"
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooAbsPdf.h"
#include "RooFitResult.h"
#include "RooMCStudy.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooPolynomial.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "TH1.h"
#include "TH2.h"
#include "TPad.h"
#include "TGraphAsymmErrors.h"
#include "RooDataHist.h"
/*You measure the b+ events in each bin of CT fitting the mass
Then you get a plot with the b+ events in CT (binned plot)
Which is background free, like the reco plot you got with the signal MC sample
Then you proceed as usual*/

using namespace std;

void mass_ct_bin(string filename)// 2017/18 MC sample is the filename
{
    //Int_t nbins=10;
    Double_t bins[] = {0.007,0.0073,0.0076,0.0079, 0.008, 0.009, 0.01,0.011,0.012, 0.013,0.014,0.015,0.016,0.017,0.018,0.019,  0.02165, 0.02865,0.03565,0.04265, 0.04965, 0.05665,  0.06365,0.07065,0.07765,0.08465,0.09165,0.09865,0.10565,0.11265,0.11965,0.12665,0.13365,0.14065,0.14765,0.15465,0.16165,0.16865,0.17565,0.18265,0.18965,0.19665,0.20365,0.21065,0.21765,0.22465,0.23165,0.23865,0.24565,0.25265,0.25965,0.26665,0.27365,0.28065,0.28765,0.29465,0.3, 0.32,0.34,0.36,0.38,0.4};
    Int_t  nbins = sizeof(bins)/sizeof(Double_t) - 1;
    TH1D *mass_h[nbins];
    
    for (int l = 0; l<nbins; l++){
        mass_h[l] = new TH1D(Form("BPMass%dth_Bin", l),Form("BPMass%dth_Bin; M_{B^{+}} (GeV); Events",l),nbins,5.18,5.4);
        
    }
    
    Float_t mass  , ct;
    
    TTree          *fChain, *copy;   //!pointer to the analyzed TTree or TChain
    TFile* f = new TFile(filename.c_str());
    fChain = new TTree;
    fChain=(TTree*)f->Get("treeFit");
    fChain->SetBranchAddress("BsCt2DMC",&ct);
    fChain->SetBranchAddress("svmass",&mass);
    
    Long64_t nentries = fChain->GetEntries();
    cout << "Start Processing " << nentries << " events" <<"\n";
    for (Long64_t jentry=0; jentry<nentries;jentry++)
    {
        fChain->GetEntry(jentry);
        if (jentry%100000==0) cout << "processing event " << jentry << "/" << nentries <<"\n";
    for (int j = 0 ; j < nbins; j++)
    {
        if (ct > ((0.4/nbins)*j) && ct < ((0.4/nbins)*(j+1)) )
        {
            
            mass_h[j]->Fill(mass);
            
                //std::cout<<mass<<"\t"<<j<<"\n";
            
        }}}
    for ( int p = 0; p< nbins ; p++){ std::cout<<((0.4/nbins)*(p+1))<<","<<endl;}
    TFile *outf = TFile::Open("bpm_ctbin18data.root","RECREATE");
    
    TCanvas *cx = new TCanvas("cx", "cx",0,0,800,600);
    
    for ( int m =0 ; m <nbins; m++)
    {
            mass_h[m]->Draw("colz");
            cx->SaveAs(Form("/Users/md/Documents/Data_MC_Sample/BsDataMC/massPlot_bin/massBP%d_th_Bin.png",m));
            mass_h[m]->Write(Form("bm18_%d_bin",m));
        
        }
    outf->Close();
    
  
}
