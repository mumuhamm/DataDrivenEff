#include "TStyle.h"
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
#include "TFile.h"
#include "TTree.h"
#include "TPad.h"
#include "TCanvas.h"


#include "RooDataHist.h"
#include "TH1D.h"
using namespace RooFit ;

Double_t ctFunction(Double_t *x,Double_t *par)
   {
      Double_t func;
        Double_t expoinv;
        Double_t lin, bolic;
        Double_t xx = x[0];
        lin         = par[0] + par[3]*xx;
        expoinv     = par[1]*(1./(1.+exp(-xx/par[2])));
        bolic            = par[3]*xx*xx;
        Double_t  poly2          = par[0] + par[1]*xx + par[2]*xx*xx + par[3]*xx*xx*xx + par[4]*xx*xx*xx*xx ;
        Double_t  sigmoid     = par[5]*(1./(1.+exp(-xx/par[6])));
        func            = poly2+sigmoid; //-bolic;
        Double_t g= lin+expoinv;
        //return g ;
        return func;
    }


void BsCt_Toy(string filename)
{

                /*auto ctform = new TF1("ctform",ctFunction,0.007,0.3,7);
         ctform->SetParameter(0,0.05);  
         ctform->SetParLimits(0,-1.,1.0);
         ctform->SetParName(0,"p0");
         ctform->SetParameter(1,0.017); 
         ctform->SetParLimits(1,-1,1);
         ctform->SetParName(1,"p1");
         ctform->SetParameter(2,-1.7);  
         ctform->SetParLimits(2,-2.0,2.0);
         ctform->SetParName(2,"p2");
         ctform->SetParameter(3,8.3);  
         ctform->SetParLimits(3,-5.0,10.0);
         ctform->SetParName(3,"p3");
         ctform->SetParameter(4,-13.1);   
         ctform->SetParLimits(4,-20.0,20.0);
         ctform->SetParName(4,"p4");
         ctform->SetParameter(5,0.014);
         ctform->SetParLimits(5,0.0,1.);
         ctform->SetParName(5,"p5");
         ctform->SetParameter(6,1.22); 
         ctform->SetParLimits(6,-2.0,2.);
         ctform->SetParName(6,"p6");*/
         
    TF1 *ctform = new TF1("ctEffFn","expo(0)*ROOT::Math::Chebyshev4(x,[2],[3],[4],[5],[6])",0.007,0.5);
    ctform->SetParameter(1,-1);
    ctform->FixParameter(2,1);


    
    


  Int_t nbins;
  std::cout<<" Give the value of bin number "<<"\n";
  cin>>nbins;

   Double_t ctpar0,ctpar1,ctpar2,ctpar3,ctpar4, ctpar5, ctpar6;
          Double_t kappa_val;

    Double_t bins[] = {0.007, 0.009, 0.011, 0.013, 0.015, 0.017, 0.019, 0.021, 0.023, 0.025,
        0.030, 0.035, 0.040, 0.045, 0.050, 0.055, 0.060, 0.065, 0.070, 0.075, 0.080, 0.085, 0.090, 0.095, 0.100,
        0.125, 0.150, 0.175, 0.200, 0.225, 0.250, 0.275, 0.300};
          Int_t  binnum = sizeof(bins)/sizeof(Double_t) - 1;

      auto hctaucut = new TH1D("hctaucut","ct(cm);ct(cm);#epsilon (ct) (a.u.)", binnum, bins);
      auto hctaunocut = new TH1D("hctaunocut","ct(cm);ct(cm);#epsilon(ct) (a.u.)", binnum, bins);
      auto resobulk = new TH1D("resobulk", "Resobulk ;Resolution ; Events ",  nbins, -0.004, 0.004);
      auto pullbulk = new TH1D("pullbulk", "pullbulk ; ; Events ",  nbins, -5.0, 5.0);

   TTree          *fChain, *copy;   //!pointer to the analyzed TTree or TChain
   TFile* f = new TFile(("/Users/ab/Documents/BsToJpsiPhi/FourTree/"+filename).c_str());
   fChain = new TTree;
   fChain=(TTree*)f->Get("treeFit");
   if (!fChain){
    cout << "No TTree found in input file, returning" << endl;
    return;
  }



  Long64_t n_entries = fChain->GetEntries();
        cout << "Start Processing " << n_entries << " events" <<"\n";
        Float_t ctreco, ctgencut, cterr;

        fChain->SetBranchAddress("BsCt2DMC",&ctreco);
        fChain->SetBranchAddress("BsCt2DMC_GEN",&ctgencut);
        fChain->SetBranchAddress("BsCt2DMCErr",&cterr);
        for (Int_t i=0;i<n_entries;i++) {
                                                           fChain->GetEntry(i);
                                                           hctaucut->Fill(ctreco);
                                                           resobulk->Fill(ctgencut-ctreco);
                                                           pullbulk->Fill((ctgencut-ctreco)/cterr);
                                                      }



  TFile *fIn2 = new TFile("/Users/ab/Documents/BsToJpsiPhi/FourTree/ntuBsGEN.root");
    if (!fIn2){return;}
    TTree* ffgen = (TTree*)fIn2->Get("OutTree");
    Int_t ff_entries = ffgen->GetEntries();
    Float_t  ctaugen;

    ffgen->SetBranchAddress("ctau_GEN",&ctaugen);
    for (Int_t j=0;j<ff_entries;j++) {

        ffgen->GetEntry(j);

        Double_t actualctgen =  ctaugen +  resobulk->GetRandom();
        hctaunocut->Fill(actualctgen);

    }

    TCanvas *ckap= new TCanvas("ckap", "ckap",0,0,800,600);
    pullbulk->Draw();
    kappa_val = pullbulk->GetRMS();
    std::cout<<"Kappa_value"<<kappa_val<<"\n";

    TCanvas *creso= new TCanvas("creso", "creso",0,0,800,600);
    resobulk->Draw();




    auto c4 = new TCanvas("c4", "c4", 0, 0, 1000, 800);
    auto hdivideIV = (TH1D*)hctaucut->Clone("hdivideIV");
    hdivideIV->Sumw2();
    hdivideIV->Divide(hctaunocut);
    hdivideIV->Draw("colz");    //colz ep
    hdivideIV->Fit(ctform,"M");
                     TFitResultPtr taure =    hdivideIV->Fit(ctform, "S");
                     TMatrixDSym cov_tau = taure->GetCovarianceMatrix();
                     Double_t chi2_tau   = taure->Chi2();
                     Double_t par0_tau   = taure->Parameter(0);
                     Double_t err0_tau   = taure->ParError(0);
                     taure->Print("V");
                     
                     
    ctpar0 = ctform->GetParameter(0);
    ctpar1 = ctform->GetParameter(1);
    ctpar2 = ctform->GetParameter(2);
    ctpar3 = ctform->GetParameter(3);
    ctpar4 = ctform->GetParameter(4);
    ctpar5 = ctform->GetParameter(5);
    ctpar6 = ctform->GetParameter(6);
      
    std::cout<<"par1"<<ctpar0<<"n";
   RooRealVar *ctp0 = new RooRealVar("ctp0","ctp0", ctpar0);
   RooRealVar *ctp1 = new RooRealVar("ctp1","ctp1", ctpar1);
   RooRealVar *ctp2 = new RooRealVar("ctp2","ctp2", ctpar2);
   RooRealVar *ctp3 = new RooRealVar("ctp3","ctp3", ctpar3);
   RooRealVar *ctp4 = new RooRealVar("ctp4","ctp4", ctpar4);  
   RooRealVar *ctp5 = new RooRealVar("ctp5","ctp5", ctpar5);  
   RooRealVar *ctp6 = new RooRealVar("ctp6","ctp6", ctpar6);  
                     
                     
    RooRealVar *BsCt2DMC = new RooRealVar ("BsCt2DMC","Bs ct", 0.007,0.3,"cm");//
    RooRealVar *BsCt2DMCErr = new RooRealVar ("BsCt2DMCErr","Bs ct Err", 0.0002,0.005,"cm");//
      RooRealVar *tau=new RooRealVar ("tau","tau",0.0441,0.03,0.055,"ps");//,0.043//origunal value in MC       
      
      
      
          
     RooRealVar *kappa = new RooRealVar("kappa", "kappa", kappa_val);
    RooFormulaVar *kapreso = new RooFormulaVar("kapreso", "kappa*BsCt2DMCErr", RooArgList(*kappa, *BsCt2DMCErr));
    RooRealVar *bias = new RooRealVar("bias","bias",0);
    RooRealVar *sigma = new RooRealVar("sigma","per-event error scale factor",1);
    RooGaussModel *gm= new RooGaussModel("gm","gauss model scaled bt per-event error", *BsCt2DMC, *bias, *sigma, *kapreso);                
                     
      RooGenericPdf *cteffFunc = new RooGenericPdf("cteffFunc", "exp(ctp0+ctp1*BsCt2DMC)*ROOT::Math::Chebyshev4(BsCt2DMC,ctp2,ctp3,ctp4,ctp5,ctp6)",
                                                 RooArgList(*ctp0, *ctp1, *ctp2, *ctp3, *ctp4, *ctp5, *ctp6, *BsCt2DMC));
  
    RooDataSet* expDataDterr1 = cteffFunc->generate(*BsCt2DMC,10000) ;
    RooDataSet* expDataDterr2 = gm->generate(*BsCt2DMCErr,10000) ;

        RooDecay *decay_gm= new RooDecay("decay_gm","decay",*BsCt2DMC, *tau, *gm, RooDecay::SingleSided) ;
        RooEffProd *modelEff=new RooEffProd("modelEffbfrconv","model with efficiency", *decay_gm, *cteffFunc) ;
        RooMCStudy* mcstudy = new RooMCStudy(*modelEff,RooArgSet(*BsCt2DMC,*BsCt2DMCErr),Binned(kTRUE),Silence(),Extended(), FitOptions(Save(kTRUE),PrintEvalErrors(0))) ;
        mcstudy->generateAndFit(1000,50000, kTRUE) ;
        
        
  // Make plots of the distributions of mean, the error on mean and the pull of mean
  RooPlot* frame1 = mcstudy->plotParam(*tau,Bins(40), Title("Decay length(ct) (cm)")) ;
  RooPlot* frame2 = mcstudy->plotError(*tau,Bins(40),Title("#Delta ct from Toy")) ;
  RooPlot* frame3 = mcstudy->plotPull(*tau,Bins(40),Title("Pull distribution of ct"),FitGauss(kTRUE)) ;
         // Plot distribution of minimized likelihood
  RooPlot* frame4 = mcstudy->plotNLL(Title("-Log(likelihood)"),Bins(40)) ;

  

  // Access some of the saved fit results from individual toys
  TH2* corrHist000 = mcstudy->fitResult(562)->correlationHist("c562") ;
  TH2* corrHist127 = mcstudy->fitResult(839)->correlationHist("c839") ;




  // Draw all plots on a canvas
  gStyle->SetPalette(1) ;
  gStyle->SetOptStat(0) ;
  TCanvas* c = new TCanvas("c","c",900,900) ;
  c->Divide(2,2) ;
  c->cd(1) ; gPad->SetLeftMargin(0.15) ; frame1->GetYaxis()->SetTitleOffset(1.4) ; frame1->Draw() ;
  c->cd(2) ; gPad->SetLeftMargin(0.15) ; frame2->GetYaxis()->SetTitleOffset(1.4) ; frame2->Draw() ;
  c->cd(3) ; gPad->SetLeftMargin(0.15) ; frame3->GetYaxis()->SetTitleOffset(1.4) ; frame3->Draw() ;
  c->cd(4) ; gPad->SetLeftMargin(0.15) ; frame4->GetYaxis()->SetTitleOffset(1.4) ; frame4->Draw() ;
        
     
              
     
}
