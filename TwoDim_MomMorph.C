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



void TwoDim_MomMorph()
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
    RooRealVar *morph_par_sig = new RooRealVar("morph_par_sig", "morph_par_sig", -1,1);
    RooRealVar *morph_par_bkg = new RooRealVar("morph_par_bkg", "morph_par_bkg", -1,1);
    
    RooDataSet *data = new RooDataSet("data", "data", RooArgSet(*svmass, *mistag), Import(*chain_data));
     
    RooRealVar *par1 = new RooRealVar("par1", "par1", 0.4);
    RooRealVar *par2 = new RooRealVar("par2", "par2", 0.6);
    RooRealVar *par3 = new RooRealVar("par3", "par3", 0.3);
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
 
    RooPlot* bsmass1 = svmass->frame(Title("M_{B_{s}} (GeV/c^{2}) with BG"),Bins(nbins));
    dataPBG->plotOn(bsmass1,DataError(RooAbsData::SumW2));
    RooPlot* misTag = mistag->frame(Title("The mistag distribution with BG"),Bins(nbins));
    dataPBG->plotOn(misTag);
    TCanvas *lz = new TCanvas();
    lz->Divide(2,1);
    lz->cd(1);bsmass1->Draw();lz->cd(2);misTag->Draw();
    
    
    
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
    std::cout<<"The Ratio SF: "<<Integral_SB/Integral_SR<<"\n";//1835.64/1605.95
    cout<<"============="<<"\n";
//================================================================================Subtract background events from signal
     RooDataSet *data_SigReg = (RooDataSet*)dataPBG->reduce(signalregion); 
     RooDataSet *data_SBReg = (RooDataSet*)dataPBG->reduce(sidebandregion);
     
     RooDataSet *data_TEST = (RooDataSet*)data->reduce(sidebandregion);
    
     TH1F * mthistSR = (TH1F*)data_SigReg->createHistogram("mistag",nbins);
     TH1F * mthistSB = (TH1F*)data_SBReg->createHistogram("mistag",nbins);
     TH1F * mthist_TEST= (TH1F*)data_TEST->createHistogram("mistag",nbins);
     TH1F *SignalMinusBG=(TH1F*)mthistSR->Clone();
     mthistSB->Scale(Integral_SB/Integral_SR);//1835.64/1605.95);//28241.7/365090
     SignalMinusBG->Add(mthistSB,-1);
     SignalMinusBG->Sumw2();     
     for (int itera=0; itera<nbins; itera++)
      {
         if (SignalMinusBG->GetBinContent(itera)<0) SignalMinusBG->SetBinContent(itera,0);
      }
    
    //================================================================================Subtract signal events from background
    
    TH1F *BGMinusSig=(TH1F*)mthistSB->Clone();
    mthistSR->Scale(0.044);
    BGMinusSig->Add(mthistSR,-1);
    BGMinusSig->Sumw2();
    for (int itera1=0; itera1<nbins; itera1++)
    {
        if (BGMinusSig->GetBinContent(itera1)<0) BGMinusSig->SetBinContent(itera1,0);
    }
    
    
    
      TCanvas *lh = new TCanvas("lh","lh",600,600);
      lh->Divide(2,1);lh->cd(1); mthistSR->Draw();lh->cd(2); mthistSB->Draw();
      TCanvas *lq = new TCanvas();
      SignalMinusBG->Draw();
      std::cout<<" Number of events in  mistag signal histogram from signal data applying side-band selection: "<<mthist_TEST->GetEntries()<<"\n";
      
   
    
   
   
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
    
    
    RooMomentMorph *morph_signal = new RooMomentMorph("morph_signal","morph_signal",*morph_par_sig,varlist,pdf_sig_mistag, paramVec,RooMomentMorph::Linear);
    morph_signal->Print("v");
    
    RooMomentMorph *morph_background = new RooMomentMorph("morph_background","morph_background",*morph_par_bkg,varlist,pdf_bkg_mistag, paramVec,RooMomentMorph::Linear);
    morph_background->Print("v");
    
    
    
    
    RooRealVar *conts = new RooRealVar("conts","conts", 1.75079,0.0, 5);
    RooExponential *expoBg = new RooExponential("expoBg","expoBG",*svmass,*conts);
    
    RooRealVar *m1par = new RooRealVar("m1par", "m1par", 0.3,0.6);
    RooRealVar *m2par = new RooRealVar("m2par", "m2par", 0.5, 0.7);
    RooRealVar *m3par = new RooRealVar("m3par", "m3par", 0.1, 0.5);
    RooGenericPdf *polyBg = new RooGenericPdf("ployBg", "@0+@1*@3+@2*@3*@3", RooArgSet(*m1par,*m2par,*m3par,*svmass));
    

    RooRealVar *mean = new RooRealVar("mean","mean",5.36679, 5.35, 5.37) ;
    RooRealVar *sigma1 = new RooRealVar("sigma1","sigma1",0.03,0.,0.5) ;
    RooRealVar *sigma2 = new RooRealVar("sigma2","sigma2",0.018,0.,0.5) ;
    RooRealVar *sigma3 = new RooRealVar("sigma3","sigma3",0.022,0.,0.5) ;
   
     RooRealVar *mu= new RooRealVar("mu", "mu", 5.36679, 5.35, 5.37);
    RooRealVar *lambda = new RooRealVar("lambda", "lambda", 0.5, 0, 5);
    RooRealVar *gamma = new RooRealVar("gamma", "gamma", 2.29776e-02, 0, 0.02);
    RooRealVar *delta = new RooRealVar("delta", "delta", 1., 0, 10);
    
    
    
    RooJohnsonLocal *john = new RooJohnsonLocal("john", "john", *svmass, *mu, *lambda, *gamma, *delta);
    //===============================================================================================
    RooGaussian *gauss1 = new RooGaussian("gauss1","gauss1",*svmass,*mean,*sigma1) ;
    RooGaussian *gauss2 = new RooGaussian("gauss2","gauss2",*svmass,*mean,*sigma2) ;
    RooGaussian *gauss3 = new RooGaussian("gauss3","gauss3",*svmass,*mean,*sigma3) ;
   

    RooRealVar *frac1 = new RooRealVar("frac1","frac1",0.1681,0.1, 0.4);
    RooRealVar *frac2 = new RooRealVar("frac2","frac2",0.3,0.1,0.60);
    RooRealVar *frac3 = new RooRealVar("frac3","frac3",0.2,0.1,1);
    
    // Create adaptive kernel estimation pdf. In this configuration the input data, is mirrored over the boundaries to minimize edge effects in distribution, that do not fall to zero towards the edges,An adaptive kernel estimation pdf on the same data without mirroring option
    //===========================================================================================================kernel estimation pdf
    //RooKeysPdf *kerestisig= new RooKeysPdf("kerestisig","kerestisig",*mistag,*data,RooKeysPdf::NoMirror);
    //RooKeysPdf *kerestiMTBG= new RooKeysPdf("kerestiMTBG","kerestiMBG",*mistag,*data_bkg,RooKeysPdf::NoMirror);
    
    //RooKeysPdf *kerestiMBG= new RooKeysPdf("kerestiMBG","kerestiMBG",*svmass,*data2,RooKeysPdf::NoMirror);
    
    
    RooAddPdf *T_gaus = new RooAddPdf("T_gaus","T_gaus",RooArgList(*gauss1,*gauss2,*gauss3),RooArgList(*frac1,*frac2));

    RooRealVar *nSig = new RooRealVar("nSig", "Number of Signal Events in SIGNAL MC",1400,0,(int)chain_data->GetEntries());
    RooRealVar *nBkg = new RooRealVar("nBkg", "Number of Backgound Events in produced  MC",800,0,988950);//(int)chain_data->GetEntries());
   
    //RooProdPdf *sigpdf = new RooProdPdf("sigpdf", "mass*mistag",RooArgList(*john,*sighistpdf));
    RooProdPdf *sigpdf = new RooProdPdf("sigpdf", "mass*mistag",RooArgList(*john,*morph_signal));
    //RooProdPdf *sigpdf = new RooProdPdf("sigpdf", "mass*mistag",RooArgList(*john,*kerestisig));
   
    
    RooProdPdf *bkgpdf = new RooProdPdf("bkgpdf", "mass*mistag",RooArgList(*expoBg,*morph_background));
    //RooProdPdf *bkgpdf = new RooProdPdf("bkgpdf", "mass*mistag",RooArgList(*expoBg,*kerestiMTBG));
   
    
    

    RooAddPdf *MTpdf = new RooAddPdf("MTpdf","MTpdf",RooArgList(*sigpdf,*bkgpdf),RooArgList(*nSig, *nBkg));
    
    
    
    
    
    RooFitResult* fitRes = MTpdf->fitTo(*dataPBG,Save(), Extended(1));//data_SigReg
    fitRes->Print("v");
    gStyle->SetOptStat(0) ;
    gStyle->SetPalette(1) ;
    TH2* hcorr = fitRes->correlationHist() ;
    TCanvas* cor= new TCanvas("Moment_Morph","MM 2D-fit correaltion matrix",800,400) ;
    gPad->SetLeftMargin(0.15) ; hcorr->GetYaxis()->SetTitleOffset(1.4) ; hcorr->Draw("colz") ;
    std::cout<<"=========================================================================="<<"\n";
    
    RooAbsReal* iMT_sig = sigpdf->createIntegral(*svmass,NormSet(*svmass),Range("signalcent")) ;
    RooAbsReal* iMT_bkg = bkgpdf->createIntegral(*svmass,NormSet(*svmass),Range("signalcent")) ;
    Double_t P_SIG_SR = iMT_sig->getVal();
    Double_t N_SIG_SR =  (iMT_sig->getVal())*(nSig->getValV()) ;
    Double_t N_SIGErr_SR = nSig->getError()*iMT_sig->getVal();
    std::cout << "Probability of getting signal events in signal region ="<< P_SIG_SR<< "\n" ;
    std::cout << "Number of signal events in the signal region  S ="<<N_SIG_SR<<"\n" ;
    std::cout<< " Error in signal count :  Serr = "<<N_SIGErr_SR<<"\n";
    
    std::cout<<"========================"<<"\n";
    Double_t P_BKG_SR =  iMT_bkg->getVal();
    Double_t N_BKG_SR =  (iMT_bkg->getVal())*(nBkg->getValV());
    Double_t N_BKGErr_SR = nBkg->getError()*iMT_bkg->getVal();
    std::cout << "Probability of getting background events in signal region ="<<P_BKG_SR << "\n" ;
    std::cout << "Number of background events in the signal region B ="<<N_BKG_SR<<"\n" ;
    std::cout<<" Error in background count : Berr = "<<N_BKGErr_SR<<"\n";
    
    
    std::cout<<"========================="<<"\n";
    RooAbsReal* iMT_sig_sb = sigpdf->createIntegral(*svmass,NormSet(*svmass),Range("sbleft,sbright")) ;
    Double_t P_SIG_SB = iMT_sig_sb->getVal();
    Double_t N_SIG_SB = (iMT_sig_sb->getVal())*(nSig->getValV());
    Double_t N_SIGErr_SB = nSig->getError()*iMT_sig_sb->getVal();
    std::cout << "Probability of getting signal events in sideband-region  ="<<P_SIG_SB<< "\n" ;
    std::cout << "Number of signal events in the side-band region  ="<<N_SIG_SB<<"\n" ;
    std::cout<<" Error in signal count in the side-band:  "<<N_SIGErr_SB<<"\n";
       
    std::cout<<"========================="<<"\n";
    RooAbsReal* iMT_bkg_sb = bkgpdf->createIntegral(*svmass,NormSet(*svmass),Range("sbleft,sbright")) ;
    Double_t P_BKG_SB =   iMT_bkg_sb->getVal();
    std::cout << "Probability of getting BKG events in sideband-region  ="<<P_BKG_SB<< "\n" ;
    std::cout << "Number of BKG events in the side-band region  ="<< (iMT_bkg_sb->getVal())*(nBkg->getValV()) <<"\n" ;
    std::cout<<" Error in BKG count in the side-band: "<<nBkg->getError()*iMT_bkg_sb->getVal()<<"\n";
   
    Double_t foM =  N_SIG_SR/(sqrt(N_SIG_SR+N_BKG_SR));
    std::cout<<" Total signal pdf integral value  I_Sig ="<<(P_SIG_SR+P_SIG_SB)<<"\n";
    std::cout<<" Total background pdf integral value  I_Bkg ="<<(P_BKG_SR+P_BKG_SB)<<"\n";
    std::cout<<" figure of merit in the signal region ="<<foM<<"\n";
    
    std::cout<<" Scale Factor  SF = " << (iMT_bkg_sb->getVal())/(iMT_bkg->getVal())<<"\n";
    std::cout<<"=========================================================================="<<"\n";
   
    RooPlot* bsmass = svmass->frame(Title("M_{B_{s}} (GeV/c^{2})"),Bins(nbins));
    dataPBG->plotOn(bsmass,DataError(RooAbsData::SumW2));
    MTpdf->plotOn(bsmass) ;
    MTpdf->paramOn(bsmass);
    //MTpdf->plotOn(bsmass, LineColor(kBlue), LineWidth(1));    
    
    RooPlot* pullframe = svmass->frame(RooFit::Title("Mass pull"));
    RooHist* hpull1 = bsmass->pullHist();
    
    pullframe->addPlotable(hpull1,"P0") ;
    pullframe->SetMinimum(-3) ;
    pullframe->SetMaximum(+3) ;
    pullframe->SetYTitle("pull");
    pullframe->SetMarkerStyle(20);
    pullframe->SetNdivisions(10);
    TCanvas * pull_Can = new TCanvas("pull_Can", "mass Pull", 800, 200);
    pullframe->Draw();
    Double_t chisquare_mass = bsmass->chiSquare();
    cout<<"Chi square of mass fit is :"<< chisquare_mass<< endl;
    MTpdf->plotOn(bsmass, Components(*john), LineColor(3), LineWidth(1), LineStyle(4));
    //MTpdf->plotOn(bsmass, Range("signalcent"),NormRange("signalcent"),Components(*john), LineColor(3), LineWidth(1), LineStyle(4));
    //MTpdf->plotOn(bsmass,Components(*gauss1), LineColor(3), LineWidth(1), LineStyle(4));
    //MTpdf->plotOn(bsmass, Components(*gauss2), LineColor(2), LineWidth(1), LineStyle(3));
    //MTpdf->plotOn(bsmass, Components(*gauss3), LineColor(6), LineWidth(1), LineStyle(5));
    MTpdf->plotOn(bsmass,Components(*expoBg), LineColor(46), LineWidth(1), LineStyle(6));
    
    RooPlot* misTagpl = mistag->frame(Title("A-kernel-estimation for mistag w/o mirroring"),Bins(nbins));
    dataPBG->plotOn(misTagpl);
    MTpdf->paramOn(misTagpl);
    //MTpdf->plotOn(misTagpl, Components(*kerestisig), LineColor(3), LineWidth(1), LineStyle(4));
    //MTpdf->plotOn(misTagpl, Components(*kerestiMTBG), LineColor(2), LineWidth(1), LineStyle(5));
    
   // MTpdf->plotOn(misTagpl, Components(*mtSigPdf), LineColor(3), LineWidth(1), LineStyle(4));
   // MTpdf->plotOn(misTagpl, Components(*mtBkgPdf), LineColor(2), LineWidth(1), LineStyle(5));
    MTpdf->plotOn(misTagpl,  Components(*morph_signal), LineColor(3), LineWidth(2), LineStyle(4));
    MTpdf->plotOn(misTagpl, Components(*morph_background), LineColor(2), LineWidth(2), LineStyle(5));
    
    TCanvas *c = new TCanvas("c", "c",0,0,600,600);
    TPad *pad1 = new TPad("pad1","pad1",0,0.33,1,1);
    TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.33);
    pad1->SetBottomMargin(0.00001);
    pad1->SetBorderMode(0);
    pad2->SetTopMargin(0.00001);
    pad2->SetBottomMargin(0.1);
    pad2->SetBorderMode(0);
    pad1->Draw();
    pad2->Draw();
    pad1->cd();
    gStyle->SetOptTitle(0);
    c->SetFillColor(0);
    c->SetBorderSize(2);
    c->SetLeftMargin(0.1422222);
    c->SetRightMargin(0.04444445);
    bsmass->SetStats(0);
    bsmass->Draw();
    auto cms1 = new TLatex(5.25, 32000, "#bf{CMS} #it{Simulations} 2017, #sqrt{s} = 13 TeV");
    cms1->SetNDC(false);
    cms1->SetTextColor(12);
    cms1->SetTextFont(42);
    cms1->SetTextSize(0.055);
    cms1-> Draw();
    pad2->cd();
    pullframe->SetStats(0);
    pullframe->Draw();
    c->cd();
    
    TCanvas* mis = new TCanvas("mis","kernelestimation",600,600) ;
    gPad->SetLeftMargin(0.15) ; misTagpl->GetYaxis()->SetTitleOffset(1.4) ;misTagpl->Draw() ;
   
   
   
   
   
   
   
   
   
   
    
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
    /*RooRealVar *mul = new RooRealVar("mul","mul",0.3);
     RooRealVar *sigma = new RooRealVar("sigma","sigma",0.4);
     RooRealVar *excon = new RooRealVar("excon","excon",0.3);
     RooExponential *p1 = new RooExponential("p1","p1",*svmass,*excon);
     RooGaussian *p2 = new RooGaussian("p2","p2",*mistag,*mul,*sigma);
     
     
     RooProdPdf * bgpdfsample = new RooProdPdf("bgpdfsample", "bgpdfsample", RooArgList(*p1, *p2));*/
   
   
   
   
   
   
   
   
   
   
   
   /*
    //TLegend* legend1 = new TLegend();
    //legend1->AddEntry(mthistSB, "BG_Combined");
    // legend1->AddEntry(mthist_TEST, "SB_Signal");
    
    
    //std::cout<<" Number of events in  mistag BG histogram from combined data : "<<mthistSB->GetEntries()<<"\n";
   
   
    //mthistSB ->SetMarkerColor(2);mthistSB->SetMarkerStyle(20);mthistSB->SetMarkerSize(1.3);mthistSB->Draw("colz");
    //mthist_TEST->SetMarkerStyle(23);mthist_TEST->SetMarkerColor(4);mthist_TEST->SetMarkerSize(1.3);mthist_TEST->Draw("same" "colz");
  
   const double SB1_L=5.24;
const double SB1_H=5.28;

const double SR_L=5.33;
const double SR_H=5.40;

const double SB2_L=5.45;
const double SB2_H=5.49;

TCut signalregion= Form("BsCt2DPVCosTheta>0.02 && BsCt2DPVCosTheta<0.3 && BsFitM>%f && BsFitM<%f",SR_L,SR_H);
TCut sidebandregion= Form("BsCt2DPVCosTheta>0.02 && BsCt2DPVCosTheta<0.3 && (BsFitM>%f && BsFitM<%f) || (BsFitM>%f && BsFitM<%f)",SB1_L,SB1_H,SB2_L,SB2_H);

tree->Draw("BsCtErr2DCostheta>>histotimerrSignal(100,0.0006,0.01)",signalregion,"goff");
tree->Draw("BsCtErr2DCostheta>>histotimerrBG(100,0.0006,0.01)",sidebandregion,"goff");
TH1F *histotimerrSignal = (TH1F*)gDirectory->Get("histotimerrSignal");
TH1F *histotimerrBG = (TH1F*)gDirectory->Get("histotimerrBG");
TH1F *SignalMinusBG=(TH1F*)histotimerrSignal->Clone();
histotimerrBG->Scale(intPoly->getVal()/intPoly1->getVal());
SignalMinusBG->Add(histotimerrBG,-1);

SignalMinusBG->Sumw2();

for (int itera=0; itera<100; itera++) {
if (SignalMinusBG->GetBinContent(itera)<0) SignalMinusBG->SetBinContent(itera,0);
}

RooDataHist errdataS("errdataS","errdataS", RooArgList(*w->var("BsCtErr2DCostheta")), RooFit::Import(*SignalMinusBG,kTRUE));
RooDataHist errdataSPlot("errdataSPlot","errdataSPlot", RooArgList(*w->var("BsCtErr2DCostheta")), RooFit::Import(*SignalMinusBG,kFALSE));
RooHistPdf errhistoS("errhistoS","errhistoS" ,RooArgList(*w->var("BsCtErr2DCostheta")),errdataS);
RooDataHist errdataBG("errdataBG","errdataBG", RooArgList(*w->var("BsCtErr2DCostheta")), RooFit::Import(*histotimerrBG,kTRUE));
RooHistPdf errhistoBG("errhistoBG","errhistoBG" ,RooArgList(*w->var("BsCtErr2DCostheta")),errdataBG);
   
   
   
  */
  
   
   /* RooDataSet *data_SigReg2 = (RooDataSet*)dataPBG->reduce(c1);
     RooRealVar *w= new RooRealVar("w","w",-1.0,1.0) ; 
     w->setVal(1.0) ; data_SigReg->addColumn(*w, kFALSE) ;
     w->setVal(-0.16) ; data_SigReg2->addColumn(*w, kFALSE) ;
     //RooDataSet *data_SigReg = new RooDataSet("data_SigReg", "data_SigReg", RooArgSet(*svmass, *mistag));
     data_SigReg->append(*data_SigReg2);    
      //data_SigReg->append(*data_SigReg2); 
      
      
   
   RooRealVar *w = new RooRealVar("w","w",0.2) ;
      RooDataSet *bkgdata = (RooDataSet*)dataPBG->reduce(c2);    
      bkgdata->addColumn(*w) ;
       RooDataSet *wdata = new RooDataSet(bkgdata->GetName(),bkgdata->GetTitle(),bkgdata,*bkgdata->get(),0,w->GetName()) ;
      wdata->Print("v");
      RooDataSet *sigdata = (RooDataSet*)dataPBG->reduce(c1);
      sigdata->Print("v");
  */
    
    }
