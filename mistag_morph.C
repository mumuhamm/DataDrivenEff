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
#include "RooJohnsonLocal.cxx"
using namespace RooFit ;
using namespace std;



void mistag_morph()
{
    
  Int_t nbins = 100;
  //std::cout<<"Give the bin numbers"<<"\n";
  //cin>>nbins;
    
    
  TChain* chain_data = new TChain("treeFit");
  chain_data->Add("/Users/ab/Documents/Data_MC_Sample/fittree_ntuBsMC2017.root");
  Int_t nevt = (int)chain_data->GetEntries();
  std::cout<<"Number of total events"<<nevt<<"\n";
  //Creating a data set which we are going to fit with the variables defined above
   
  RooRealVar *svmass= new RooRealVar("svmass", "M_{B_{s}} GeV/c^{2}",5.25,5.49);
  RooRealVar *mistag = new RooRealVar("mistag","Mistag fraction of original B and Bbar",0.04,0.70);
  RooRealVar *morph_par_sig = new RooRealVar("morph_par_sig", "morph_par_sig", 0.0,1);
  RooRealVar *morph_par_bkg = new RooRealVar("morph_par_bkg", "morph_par_bkg", 0.0,0.5);
    
  RooDataSet *data = new RooDataSet("data", "data", RooArgSet(*svmass, *mistag), Import(*chain_data));
     
  
  RooRealVar *par4 = new RooRealVar("par4", "par4", 0.4);
  RooExponential *bg_mass_model = new RooExponential("bg_mass_model","bg_mass_model",*svmass,*par4);
  //RooGenericPdf *bg_mass_model = new RooGenericPdf("bg_mass_model","mass bg formula", "par1+par2*@3+par3*@3*@3", RooArgSet(*par1,*par2,*par3,*svmass));
    
  RooRealVar *mean_mistag = new RooRealVar("mean_mistag", "mean_mistag", 0.3);
  RooRealVar *sigma_mistag = new RooRealVar("sigma_mistag", "sigma_mistag", 0.2);
  RooGaussian *bg_mistag_model = new RooGaussian("bg_mistag_model","bg_mistag_model",*mistag,*mean_mistag, *sigma_mistag);// RooArgSet(*par11,*par22,*par33,*mistag));
  RooProdPdf *bgpdfsample = new RooProdPdf("bgpdfsample", "bgpdfsample", RooArgList(*bg_mass_model, *bg_mistag_model));
    
    
  RooDataSet* data_bkg = bgpdfsample->generate(RooArgSet(*svmass,*mistag),98895);
    
  RooDataSet *dataPBG = new RooDataSet("dataPBG", "dataPBG", RooArgSet(*svmass, *mistag));
  dataPBG->append(*data);    
  dataPBG->append(*data_bkg);
  dataPBG->Print("v");
 
  
    
    
    
  const double SB1_L=5.24;
  const double SB1_H=5.28;
    
  const double SR_L=5.33;
  const double SR_H=5.40;
    
  const double SB2_L=5.45;
  const double SB2_H=5.49;
    
  svmass->setRange("sbleft",SB1_L,SB1_H);//SB1
  svmass->setRange("sbright",SB2_L,SB2_H);  //SB2
  svmass->setRange("signalcent",SR_L,SR_H);//Signal
  //svmass->setRange("fullragne", 5.24, 5.49);//FullRange
    
  TCut signalregion= Form(" svmass>%f && svmass<%f",SR_L,SR_H);
  TCut sidebandregion= Form(" (svmass>%f && svmass<%f) || (svmass>%f && svmass<%f)",SB1_L,SB1_H,SB2_L,SB2_H);
    
  RooAbsReal* integral_mass1 = bg_mass_model->createIntegral(*svmass,NormSet(*svmass),Range("signalcent")) ;
  RooAbsReal* integral_mass2 = bg_mass_model->createIntegral(*svmass,NormSet(*svmass),Range("sbleft,sbright")) ;
    
  cout<<"============="<<"\n";
  Double_t Integral_SR =integral_mass1->getVal();
  Double_t Integral_SB =integral_mass2->getVal();
  std::cout<<"Side Band region integral value: "<<Integral_SB<<"\n";
  std::cout<<"Signal region integral value: "<<Integral_SR<<"\n";
  std::cout<<"The Ratio SF: "<<Integral_SB/Integral_SR<<"\n";
  cout<<"============="<<"\n";

  RooDataSet *data_SigReg = (RooDataSet*)dataPBG->reduce(signalregion); 
  RooDataSet *data_SBReg = (RooDataSet*)dataPBG->reduce(sidebandregion);
     
  RooDataSet *data_TEST = (RooDataSet*)data->reduce(sidebandregion);
     
  TH1F * mthistSR = (TH1F*)data_SigReg->createHistogram("mistag",nbins);
  TH1F * mthistSB = (TH1F*)data_SBReg->createHistogram("mistag",nbins);
  
  TH1F *SignalMinusBG=(TH1F*)mthistSR->Clone();
  mthistSB->Scale(Integral_SB/Integral_SR);//1835.64/1605.95);//28241.7/365090
  SignalMinusBG->Add(mthistSB,-1);
  SignalMinusBG->Sumw2();     
  for (int itera=0; itera<nbins; itera++)
    {
      if (SignalMinusBG->GetBinContent(itera)<0) SignalMinusBG->SetBinContent(itera,0);
    }
          
 
  TCanvas *lq = new TCanvas();
  SignalMinusBG->Draw();
 
      
   
    
   
   
  RooDataHist * bkgdatahist = new RooDataHist("bkgdatahist", "bkgdatahist", *mistag, Import(*mthistSB));
  RooHistPdf * bkghistpdf = new RooHistPdf("bkghistpdf", "bkghistpdf", *mistag, *bkgdatahist,0);
  RooDataHist * sigdatahist = new RooDataHist("sigdatahist", "sigdatahist", *mistag, Import(*SignalMinusBG));
  RooHistPdf * sighistpdf = new RooHistPdf("sighistpdf", "sighistpdf", *mistag, *sigdatahist,0);
  
    TVectorT<double> paramVec = TVectorD(1);
    RooArgList pdf_sig_mistag;
    pdf_sig_mistag.add(*sighistpdf);
    pdf_sig_mistag.Print();
    RooArgList pdf_bkg_mistag;
    pdf_bkg_mistag.add(*bkghistpdf);
    pdf_bkg_mistag.Print();
    
    RooArgList varlist;
    varlist.add(*mistag);
    
    RooRealVar *nSig = new RooRealVar("nSig", "Number of Signal Events in SIGNAL MC",1400,0,(int)chain_data->GetEntries());
    RooRealVar *nBkg = new RooRealVar("nBkg", "Number of Backgound Events in produced  MC",800,0,988950);//(int)chain_data->GetEntries());
    RooMomentMorph *morph_signal = new RooMomentMorph("morph_signal","morph_signal",*morph_par_sig,varlist,pdf_sig_mistag, paramVec,RooMomentMorph::Linear);
    morph_signal->Print("v");
    
    RooMomentMorph *morph_background = new RooMomentMorph("morph_background","morph_background",*morph_par_bkg,varlist,pdf_bkg_mistag, paramVec,RooMomentMorph::Linear);
    morph_background->Print("v");
    
    RooAddPdf *morph = new RooAddPdf("morph","morph",RooArgList(*morph_signal,*morph_background),RooArgList(*nSig, *nBkg));
   
    RooFitResult* fitRes = morph->fitTo(*dataPBG,Save(), Extended(1));//data_SigReg
    fitRes->Print("v");
    gStyle->SetOptStat(0) ;
    gStyle->SetPalette(1) ;
    TH2* hcorr = fitRes->correlationHist() ;
    TCanvas* c = new TCanvas("Moment_Morph","MM correaltion matrix",800,400) ;
    gPad->SetLeftMargin(0.15) ; hcorr->GetYaxis()->SetTitleOffset(1.4) ; hcorr->Draw("colz") ;
    
    RooPlot* misTagpl = mistag->frame(Title("Moment Morphing of mistag"),Bins(nbins));
    dataPBG->plotOn(misTagpl);
    morph->paramOn(misTagpl);
    morph->plotOn(misTagpl,  Components(*morph_signal), LineColor(3), LineWidth(2), LineStyle(4));
    morph->plotOn(misTagpl, Components(*morph_background), LineColor(2), LineWidth(2), LineStyle(5));
    
    
    TCanvas* mis = new TCanvas("mis","kernelestimation",600,600) ;
    gPad->SetLeftMargin(0.15) ; misTagpl->GetYaxis()->SetTitleOffset(1.4) ;misTagpl->Draw() ;
    
}
    
    
    
