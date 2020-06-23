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
//#include "RooJohnsonLocal.cxx"
using namespace RooFit ;
using namespace std;



void mass_fit_bin(){
    gStyle->SetOptStat(0) ;
      gStyle->SetPalette(1) ;
    
    TH1D *m_bin1[60];
    RooDataHist *bindata[60];
    RooRealVar *nSig[60];
    RooRealVar *nBkg_Combi[60];
    RooRealVar *nBkg_Prompt[60];
    RooFitResult* fitRes[60]; TH2* hcorr[60];
    RooPlot* bsmass[60];
    TLatex *cms1[60];TLegend *leg[60];
    RooAddPdf *model[60];
    TRandom3 *decay = new TRandom3(456);
    Double_t time = 4.9110e-02;//2018 pdg Bplus 491.1micron
    Int_t nbins = 60;
    int ltbnum = 60; double ltbinerr[ltbnum], sigevterrnorm[ltbnum], sigevtnorm[ltbnum];
    Double_t histobins[60], histobinserr[60];
    
    TFile *myf = TFile::Open("/Users/md/Documents/Data_MC_Sample/BsDataMC/bpm_ctbin17data.root","READ");
    //TCanvas *cc = new TCanvas("cc","cc",0,0,600,600);
    
    for(int k =0 ; k < 60; k++){
       m_bin1[k] = (TH1D*)myf->Get(Form("bm17_%d_bin", k));
        //m_bin1[k]->Draw();
        //cc->SaveAs(Form("/Users/md/Documents/Data_MC_Sample/BsDataMC/massPlot_bin/pltfrmfit/mplot_%d.pdf",k));
   }
   
  RooRealVar *massvar = new RooRealVar("massvar", "massvar", 5.18,5.40);
  for(int k1 =0 ; k1 < 60; k1++){
  bindata[k1] = new RooDataHist("bindata", "bindata", *massvar, Import(*m_bin1[k1]));
  }
    
  RooRealVar *conts = new RooRealVar("conts","conts", -3.0, 3.0);
  RooExponential *expoBg = new RooExponential("expoBg","expoBG",*massvar,*conts);
   
    RooRealVar *mean= new RooRealVar("mean", "mean", 5.27929, 5.2, 5.35);
    RooRealVar *sigma = new RooRealVar("sigma", "sigma",0.03,0.,2.5);
    RooGaussian *gauss= new RooGaussian("gauss","gauss",*massvar,*mean,*sigma) ;
    
    
  RooRealVar *prompt_p1 =new RooRealVar("prompt_p1","prompt background mass polynomial", 0.1, -0.2, 0.3);
  RooPolynomial *prompt_mass =new RooPolynomial("prompt_mass", " prompt background mass",*massvar, RooArgList(*prompt_p1));
    
    
  RooRealVar *mu= new RooRealVar("mu", "mu", 5.27929, 5.2, 5.35);
  RooRealVar *lambda = new RooRealVar("lambda", "lambda",  0.5, 0, 1);
  RooRealVar *gamma = new RooRealVar("gamma", "gamma",  0., -1, 1);
  RooRealVar *delta = new RooRealVar("delta", "delta", 1., 0, 10);
  RooJohnson *john = new RooJohnson("john", "john", *massvar, *mu, *lambda, *gamma, *delta);

    
    
     for(int k3 =0 ; k3 < 60; k3++){
       nSig[k3] = new RooRealVar("nSig", "Number of Signal Events in SIGNAL MC",7,(int)m_bin1[k3]->GetEntries());//e05 for all era 3 for 17D
       nBkg_Combi[k3] = new RooRealVar("nBkg_Combi", "Number of Backgound Events in produced  MC", 1 ,(int)m_bin1[k3]->GetEntries());//1000 for all other era
       nBkg_Prompt[k3] = new RooRealVar("nBkg_Prompt", "Number of Backgound Events in produced  MC", 1 ,(int)m_bin1[k3]->GetEntries());//500 for all other era
     
 // RooAddPdf *model = new RooAddPdf("model","model",RooArgList(*john,*expoBg, *prompt_mass), RooArgList(*nSig,*nBkg_Combi,*nBkg_Prompt));
    
    
        model[k3] = new RooAddPdf("model","model",RooArgList(*john,*expoBg), RooArgList(*nSig[k3],*nBkg_Combi[k3]));
    }
    
    
      /* double ltbin[62]={0.007 ,0.0073,0.0076,0.0079, 0.008, 0.009, 0.01,0.011,0.012, 0.013,0.014,0.015,0.016,0.017,0.018,0.019,  0.02165, 0.02865,0.03565,0.04265, 0.04965, 0.05665,  0.06365,0.07065,0.07765,0.08465,0.09165,0.09865,0.10565,0.11265,0.11965,0.12665,0.13365,0.14065,0.14765,0.15465,0.16165,0.16865,0.17565,0.18265,0.18965,0.19665,0.20365,0.21065,0.21765,0.22465,0.23165,0.23865,0.24565,0.25265,0.25965,0.26665,0.27365,0.28065,0.28765,0.29465,0.3,0.32,0.34,0.36,0.38, 0.4};*/
    
    
    double ltbin[] = {0.00655738,0.0131148,0.0196721,0.0262295,0.0327869,0.0393443,0.0459016,0.052459,0.0590164,0.0655738,0.0721311,0.0786885,0.0852459,0.0918033,0.0983607,0.104918,0.111475,0.118033,0.12459,0.131148,0.137705,0.144262,0.15082,0.157377,0.163934,0.170492,0.177049,0.183607,0.190164,0.196721,0.203279,0.209836,0.216393,0.222951,0.229508,0.236066,0.242623,0.24918,0.255738,0.262295,0.268852,0.27541,0.281967,0.288525,0.295082,0.301639,0.308197,0.314754,0.321311,0.327869,0.334426,0.340984,0.347541,0.354098,0.360656,0.367213,0.37377,0.380328,0.386885,0.393443,0.4};
    
    Int_t  binnum = sizeof(ltbin)/sizeof(Double_t)-1;
    
    TCanvas *c = new TCanvas("c", "c",0,0,600,600);
    auto recotau = new TH1D("recotau","Reco ct(BG removed) ;ct[cm];#epsilon [a.u.]",  binnum, ltbin);
    auto gentau = new TH1D("gentau","Coupling independent Gen ct ;ct[cm];#epsilon [a.u.]",  binnum, ltbin);
    
   for(int k2 =0 ; k2 < 60; k2++){
       fitRes[k2] = model[k2]->fitTo(*bindata[k2],Save(),NumCPU(8));//,ConditionalObservables(*BsCt2DMCErr));//data_SigReg
       fitRes[k2]->Print("v");
       RooRealVar* par1_fitresult = (RooRealVar*)fitRes[k2]->floatParsFinal().find("nSig");
       sigevtnorm[k2] = par1_fitresult->getVal()/(ltbin[k2+1]-ltbin[k2]);
       sigevterrnorm[k2] = par1_fitresult->getAsymErrorHi()/(ltbin[k2+1]-ltbin[k2]);
       recotau->SetBinContent(k2,sigevtnorm[k2]);
       recotau->SetBinError(k2, sigevterrnorm[k2]);
       ltbinerr[k2] = 0;
       std::cout<<par1_fitresult->getVal()<<"\t"<<par1_fitresult->getAsymErrorHi()<<"\n";
       bsmass[k2] = massvar->frame(Title("M_{B^{+}} (GeV/c^{2})"),Bins(40));
       bindata[k2]->plotOn(bsmass[k2],DataError(RooAbsData::SumW2));
       model[k2]->plotOn(bsmass[k2]) ;
       model[k2]->paramOn(bsmass[k2]);
       model[k2]->plotOn(bsmass[k2], RooFit::LineColor(kGreen),RooFit::Components("john"), RooFit::Name("signal"), LineWidth(2), LineStyle(4));
       model[k2]->plotOn(bsmass[k2],RooFit::LineColor(kRed),RooFit::Components("expoBg"), RooFit::Name("combinatorial"), LineWidth(2), LineStyle(6));
           leg[k2] = new TLegend(0.7,0.7,0.9,0.9);
           leg[k2]->AddEntry(bsmass[k2]->findObject("signal"),"B^{+}#rightarrow J/#psi K^{+}","l");
           leg[k2]->AddEntry(bsmass[k2]->findObject("combinatorial"),"Combinatorial","l");
           bsmass[k2]->SetStats(0);
           bsmass[k2]->Draw();
           leg[k2]->Draw("same");
           cms1[k2] = new TLatex(5.18, 5300, "#bf{CMS} #it{Preliminary} 2017, M_{B^{+}} in #it{ct-bin}, #sqrt{s} = 13 TeV");
           cms1[k2]->SetNDC(false);
           cms1[k2]->SetTextColor(12);
           cms1[k2]->SetTextFont(42);
           cms1[k2]->SetTextSize(0.055);
           cms1[k2]-> Draw();
           c->Update();
           c->SaveAs(Form("/Users/md/Documents/Data_MC_Sample/BsDataMC/massPlot_bin/pltfrmfit/fittedM_%d.pdf",k2));
       
   }
    
    
       TCanvas *cmlt= new TCanvas("cmlt", "cmlt",0,0,800,600);
       cmlt->SetLogy();
       auto* gr_mct = new TGraphErrors(ltbnum, ltbin,sigevtnorm,ltbinerr,sigevterrnorm);
       gr_mct->GetXaxis()->SetTitle("ct");
       gr_mct->GetYaxis()->SetTitle("Signal events of M_{B^{+}} in ct bin ");
       gr_mct->SetMarkerColor(kBlue);
       gr_mct->SetMarkerStyle(20);
       gr_mct->SetLineWidth(2);
       gr_mct->SetLineColor(kBlue);
       gr_mct->Draw("AP");
    
    
    double ltbin_my[] = {0.00655738,0.0131148,0.0196721,0.0262295,0.0327869,0.0393443,0.0459016,0.052459,0.0590164,0.0655738,0.0721311,0.0786885,0.0852459,0.0918033,0.0983607,0.104918,0.111475,0.118033,0.12459,0.131148,0.137705,0.144262,0.15082,0.157377,0.163934,0.170492,0.177049,0.183607,0.190164,0.196721,0.203279,0.209836,0.216393,0.222951,0.229508,0.236066,0.242623,0.24918,0.255738,0.262295,0.268852,0.27541,0.281967,0.288525,0.295082,0.301639,0.308197,0.314754,0.321311,0.327869,0.334426,0.340984,0.347541,0.354098,0.360656,0.367213,0.37377,0.380328,0.386885,0.393443,0.4};
    
    
    std::cout<<binnum<<"\n";
    auto bdecaytime = new TH1D("bdecaytime","Gen lifetime ;ct[cm];#epsilon [a.u.]",  binnum, ltbin);
    TFile *fb0 = new TFile("/Users/md/Documents/BsToJpsiPhi/Angle_Eff/fittree_ntuBdMC2018.root");
    auto resobulk = new TH1D("resobulk", "Resobulk ; From B_{d}^{0}; Events ",  nbins, -0.004, 0.004);
    auto pullbulk = new TH1D("pullbulk", "pullbulk ; From B_{d}^{0}; Events ",  nbins, -5.0, 5.0);
    
    TTree* b0file = (TTree*)fb0->Get("treeFit");
    Int_t n_entriesb0 = b0file->GetEntries();
    Float_t ctrecob0, ctgenb0, cterrb0;
    b0file->SetBranchAddress("BsCt2DMC",&ctrecob0);
    b0file->SetBranchAddress("BsCt2DMC_GEN",&ctgenb0);
    b0file->SetBranchAddress("BsCt2DMCErr",&cterrb0);
    for( int k = 0 ; k < n_entriesb0 ; k++ )
    {
    b0file->GetEntry(k);
        resobulk->Fill(ctgenb0-ctrecob0);
        pullbulk->Fill((ctgenb0-ctrecob0)/cterrb0);
    }
    
    for( int i = 0 ; i < 41941700 ; i++ )
    {
        Double_t ctaugen1 =  decay->Exp(time);
        Double_t actualctgen1 =  ctaugen1 + resobulk->GetRandom();
        bdecaytime->Fill(actualctgen1);
    }
    
   for(int k5 =0 ; k5 < binnum; k5++)
   {
       //std::cout<<k5<<"\n";
       histobins[k5] = (bdecaytime->GetBinContent(k5));///(ltbin[k5+1]-ltbin[k5]);
       histobinserr[k5] =(bdecaytime->GetBinError(k5));///(ltbin[k5+1]-ltbin[k5]);
       gentau->SetBinContent(k5,histobins[k5]);
       gentau->SetBinError(k5, histobinserr[k5]);
       std::cout<<k5<<"\t"<<histobins[k5]<<"\t"<<histobinserr[k5]<<"\t"<<ltbin[k5]<<"\n";//<<ltbinerr[k5]<<"\n";
   }
    
    TCanvas *cmlt_gen= new TCanvas("cmlt_gen", "cmlt_gen",0,0,800,600);
    cmlt_gen->SetLogy();
    auto* gr_mct_gen = new TGraphErrors(binnum,ltbin,histobins, ltbinerr, histobinserr);
    gr_mct_gen->GetXaxis()->SetTitle("ct");
    gr_mct_gen->GetYaxis()->SetTitle("Signal events of M_{B^{+}} in ct bin ");
    gr_mct_gen->SetMarkerColor(kRed);
    gr_mct_gen->SetMarkerStyle(21);
    gr_mct_gen->SetLineWidth(2);
    gr_mct_gen->SetLineColor(kRed);
    gr_mct_gen->Draw("AP");
    gr_mct_gen->Print();
    TCanvas *cm = new TCanvas("cm","Multigraph",200,10,700,700);
    TMultiGraph *mg = new TMultiGraph();
    mg->Add(gr_mct);
    mg->Add(gr_mct_gen);
    mg->Draw("AC*");
    
    TLegend *legend = new TLegend(0.1,0.7,0.48,0.9);
    //legend->SetHeader("Header","C");
    legend->AddEntry(gr_mct,"From Data");
    legend->AddEntry(gr_mct_gen,"From TRandom");
    legend->Draw();
    cm->Update();
    
    TCanvas *cm1 = new TCanvas("cm1","cm1",0,0,900,700);
    cm1->Divide(2,1);
    cm1->cd(1);recotau->Draw();
    cm1->cd(2);gentau->Draw();
    
    
    auto ctform = new TF1("ctEffFn","expo(0)*ROOT::Math::Chebyshev4(x,[2],[3],[4],[5],[6])",0.007,0.4);
    ctform->SetParameter(1,-1);
    ctform->FixParameter(2,1);
    
    TCanvas *ceff= new TCanvas("ceff", "ceff",0,0,600,600);
    auto hdivideI = (TH1D*)recotau->Clone("hdivideI");
    hdivideI->Sumw2();
    hdivideI->Divide(gentau);
    hdivideI->Draw("colz");    //colz ep
   /* hdivideI->Fit(ctform);
    TFitResultPtr r =    hdivideI->Fit(ctform, "S");
       TMatrixDSym cov = r->GetCovarianceMatrix();   //  to access the covariance matrix
       Double_t chi2   = r->Chi2();                  // to retrieve the fit chi2
       Double_t par0   = r->Parameter(0);            // retrieve the value for the parameter 0
       Double_t err0   = r->ParError(0);             // retrieve the error for the parameter 0
       r->Print("V");                                // print full information of fit including covariance matrix
    
  
    TCanvas *ccor = new TCanvas("Corr Matrix","M/ct/cterr correaltion matrix",800,400) ;
    hcorr[k2] = fitRes[k2]->correlationHist() ;
          gPad->SetLeftMargin(0.15) ; hcorr[k2]->GetYaxis()->SetTitleOffset(1.4) ; hcorr[k2]->Draw("colz") ;
          ccor->SaveAs(Form("/Users/md/Documents/Data_MC_Sample/BsDataMC/massPlot_bin/pltfrmfit/CorrelationMatrixBinFit_%d.pdf",k2));
    RooPlot* pullframe = massvar->frame(RooFit::Title("Mass pull"));
    RooHist* hpull1 = bsmass->pullHist();
    pullframe->addPlotable(hpull1,"P0") ;
    pullframe->SetMinimum(-3) ;
    pullframe->SetMaximum(+3) ;
    pullframe->SetYTitle("pull");
    pullframe->SetMarkerStyle(20);
    pullframe->SetNdivisions(10);
    Double_t chisquare_mass = bsmass->chiSquare();
    
    model->plotOn(bsmass, RooFit::LineColor(kGreen),RooFit::Components("john"), RooFit::Name("signal"), LineWidth(2), LineStyle(4));
    model->plotOn(bsmass,RooFit::LineColor(kRed),RooFit::Components("expoBg"), RooFit::Name("combinatorial"), LineWidth(2), LineStyle(6));
    //model->plotOn(bsmass,RooFit::LineColor(kMagenta),RooFit::Components("prompt_mass"), RooFit::Name("prompt"), LineWidth(2), LineStyle(6));
    TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);
    leg->AddEntry(bsmass->findObject("signal"),"B^{+}#rightarrow J/#psi K^{+}","l");
    leg->AddEntry(bsmass->findObject("combinatorial"),"Combinatorial","l");
    //leg->AddEntry(bsmass->findObject("prompt"),"B^{+}#rightarrow J/#psi X","l");
    
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
    leg->Draw("same");
    auto cms1 = new TLatex(5.18, 5300, "#bf{CMS} #it{Preliminary} 2017, M_{B^{+}} in #it{ct-bin}, #sqrt{s} = 13 TeV");
    cms1->SetNDC(false);
    cms1->SetTextColor(12);
    cms1->SetTextFont(42);
    cms1->SetTextSize(0.055);
    cms1-> Draw();
    pad2->cd();
    pullframe->SetStats(0);
    pullframe->Draw();
    c->cd();
    cout<<"Chi square of mass fit is :"<< chisquare_mass<< endl;
    
    
    */
    

}
