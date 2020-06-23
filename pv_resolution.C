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



Double_t kappafit(Double_t *x,Double_t *par)
{
    Double_t func;
    Double_t expoinv;
    Double_t lin, bolic;
    Double_t xx = x[0];
    lin         = par[0] + par[1]*xx;//+ par[2]*xx*xx*xx;// + par[5]*xx*xx*xx*xx;
    expoinv     = par[1]*(1./(1.+exp(-xx/par[2])));
    bolic            = par[3]*xx*xx;
    Double_t  poly2          = par[0] + par[1]*xx + par[2]*xx*xx + par[3]*xx*xx*xx + par[4]*xx*xx*xx*xx ;
    Double_t  sigmoid     = par[0]*(1./(1.+exp(-xx/par[1])));
    func            = poly2+sigmoid; //-bolic;
    return sigmoid;
}


void pv_resolution(string filename)// 2017/18 MC sample is the filename
{
    auto myfunction = new TF1("myfunction",kappafit,5.0,30.0,2);
    myfunction->SetParameter(0,1.00);  // par[0] : offset
    myfunction->SetParLimits(0,-1000.,1000.0);
    myfunction->SetParName(0,"p0");
    myfunction->SetParameter(1,1.5); // par[1] : width at f(x) minimum
    myfunction->SetParLimits(1,-1000.0,1005.);
    myfunction->SetParName(1,"p1");
   
    
    Int_t nbins=10;
    //std::cout<<" Give the value of bin number "<<"\n";
    //cin>>nbins;
    Double_t  ct_time[nbins], rmsreso[nbins], kappa_value[nbins], npvarr[nbins], kappa_valueError[nbins], pverror[nbins];
    TH1D *resolution[nbins];
    TH1D *pull[nbins];
    TH1D *kappa_PV[nbins];
    TH1D *NPV_ar[nbins];
    
    for (int l = 0; l<nbins; l++){
        resolution[l] = new TH1D(Form("Resolution%dth_Bin", l),Form("Resolution %dth_Bin; ct_{GEN}-ct_{RECO}; Events",l),nbins,-0.004,0.004);
        kappa_PV[l] = new TH1D(Form("KappaVsPVN%dth_Bin", l),Form("Kappa_In_Bins_Of_PVMultiplicity%dth_Bin; ct_{GEN}-ct_{RECO}/#Delta ct; Events",l),  nbins, -5., 5.);
        pull[l] = new TH1D(Form("Pull%dth_Bin", l),Form("Pull%dth_Bin; ct_{GEN}-ct_{RECO}/#Delta ct; Events",l),  nbins, -5., 5.);
        NPV_ar[l] = new TH1D(Form("npv%dth_Bin", l),Form("npv%dth_Bin; Multiplicity of PV; Events",l),  nbins, 5, 30);
    }
    auto lifetimegen = new TH1D("lifetimegen", "Decay Length GEN; ct (cm); Events ",  nbins, 0.007, 0.3);
    auto resobulk = new TH1D("resobulk", "Resobulk ; ; Events ",  nbins, -0.004, 0.004);
    auto pullbulk = new TH1D("pullbulk", "pullbulk ; ; Events ",  nbins, -5.0, 5.0);
    auto primaryvtx = new TH1D("primaryvtx", "primaryvtx; PV multiplicity (); events", nbins, 0, 100);
    auto pv_vs_kappa = new TProfile("pv_vs_kappa", "#kappa(N_{PV}); pv multiplicity ();#kappa ()",20,0,30, -2.2,2.2);
    Float_t ctreco  , ctgen,  cterr, ctgencut;
    Int_t npv;
    Double_t p_vtx, pull_measure;
    Double_t kappa_val;
    TTree          *fChain, *copy;   //!pointer to the analyzed TTree or TChain
    TFile* f = new TFile(filename.c_str());
    fChain = new TTree;
    fChain=(TTree*)f->Get("treeFit");
    fChain->SetBranchAddress("BsCt2DMC",&ctreco);
    fChain->SetBranchAddress("BsCt2DMC_GEN",&ctgen);
    fChain->SetBranchAddress("BsCt2DMCErr",&cterr);
    fChain->SetBranchAddress("Bs_NPV",&npv);
    Long64_t nentries = fChain->GetEntries();
    cout << "Start Processing " << nentries << " events" <<"\n";
    for (Long64_t jentry=0; jentry<nentries;jentry++)
    {
        fChain->GetEntry(jentry);
        if (jentry%100000==0) cout << "processing event " << jentry << "/" << nentries <<"\n";
        if (cterr !=cterr)continue;
        if (ctreco !=ctreco)continue;
        p_vtx = npv;
        pull_measure = ((ctgen-ctreco)/cterr*1.0);
        
        resobulk->Fill(ctgen-ctreco);
        pullbulk->Fill((ctgen-ctreco)/cterr);
        
        primaryvtx->Fill(npv);
        for (int j1 = 0 ; j1 < nbins; j1++){
            if (npv > ((30/nbins)*j1) && npv < ((30/nbins)*(j1+1)) )
            {
               
                kappa_PV[j1]->Fill((ctgen-ctreco)/cterr);
                
            }
        }
        for (int j = 0 ; j < nbins; j++)
        {
            if (ctreco > ((0.3/nbins)*j) && ctreco < ((0.3/nbins)*(j+1)) )
            {
                resolution[j]->Fill(ctgen-ctreco);
                pull[j]->Fill((ctgen-ctreco)/cterr);
                
            }}}
    
    TCanvas *cx = new TCanvas("cx", "cx",0,0,800,600);
    for ( int m =0 ; m <nbins; m++)
    {
        std::cout<<"Rms values in each bin for Resolution"<<m<<"\t"<<resolution[m]->GetRMS()<<"\n";
        //std::cout<<"Rms values in each bin for Resolution"<<m<<"\t"<<resolution[m]->GetRMS()<<"\n";
        std::cout<<"Rms values in each bin for ct Pull"<<m<<"\t"<<pull[m]->GetRMS()<< "Kappa factor :"<<(pull[m]->GetRMS())<<"\n";
        //resolution[m]->Draw("colz");
        // cx->SaveAs(Form("/eos/user/m/mumuhamm/codes/AllEffResolution/BsResolution%d_th_Bin.png",m));
        kappa_PV[m]->Draw("colz");
        kappa_PV[m]->Fit("gaus");
        cx->SaveAs(Form("/Users/md/Documents/BsToJpsiPhi/Angle_Eff/Plots/kappa_pv/Kappa_In_Bins_Of_PVMultiplicity%d_th_Bin.png",m));
    }
    for(int z = 0; z<nbins; z++)
    {
       
        rmsreso[z]= resolution[z]->GetRMS();
        ct_time[z]= (0.3/nbins) + (0.3/nbins)*z;
    }
    for(int z1 = 0; z1 < nbins; z1++){
        kappa_value[z1] = (kappa_PV[z1]->GetRMS()) ;
        kappa_valueError[z1]= (kappa_PV[z1]->GetRMSError());
        std::cout<<"print the kappa value "<<kappa_value[z1]<<"Bin error"<<kappa_valueError[z1]<<"\n";
        npvarr[z1]= (30/nbins) + (30/nbins)*z1;
        std::cout<<"print pv number"<<npvarr[z1]<<"\n";
           pverror[z1]= 0;
      
    }
 
   // TCanvas *ce= new TCanvas("ce", "ce",0,0,800,600);
    //TGraph* gr_resoct = new TGraph(nbins, ct_time,rmsreso);
    //gr_resoct->Draw("AC");
    TCanvas *c1= new TCanvas("c1", "c1",0,0,800,600);
    auto* gr_kappapv = new TGraphErrors(nbins, npvarr,kappa_value, pverror, kappa_valueError);
    gr_kappapv->Fit(myfunction);
    TFitResultPtr r =    gr_kappapv->Fit(myfunction, "S");
    TMatrixDSym cov = r->GetCovarianceMatrix();   //  to access the covariance matrix
    Double_t chi2   = r->Chi2();                  // to retrieve the fit chi2
    Double_t par0   = r->Parameter(0);            // retrieve the value for the parameter 0
    Double_t err0   = r->ParError(0);             // retrieve the error for the parameter 0
    r->Print("V");                                // print full information of fit including covariance matrix
     gr_kappapv->GetXaxis()->SetTitle("PV multiplicity");
    gr_kappapv->GetYaxis()->SetTitle("#kappa = #frac{#delta(ct)}{#Delta(ct)}");
    gr_kappapv->Draw("AP");
    auto cms1 = new TLatex(0.0, 1.25, "#bf{CMS} #it{Simulations} 2017, #sqrt{#bf{s}} = #bf{13 TeV}");
    cms1->SetNDC(false);
    cms1->SetTextColor(12);
    cms1->SetTextFont(42);
    cms1->SetTextSize(0.055);
    cms1-> Draw();
    c1->cd();
    TCanvas *c3= new TCanvas("c3", "c3",0,0,800,600);
    pullbulk->Draw();
    kappa_val = pullbulk->GetRMS();
    std::cout<<" kappa value for main pull"<<pullbulk->GetRMS()<<"\n";
    TCanvas *c4= new TCanvas("c4", "c4",0,0,800,600);
    resobulk->Draw();
  
}
