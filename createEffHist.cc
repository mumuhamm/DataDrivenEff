//Author Md Alibordi

 #include "RooRealVar.h"
 #include "RooDataSet.h"
 #include "TCanvas.h"
 #include "TAxis.h"
 #include "TMath.h"
 #include "RooPlot.h"

using namespace RooFit;
using namespace std;

// want to produce and save plots of distributions and efficiency?
bool plot = true;

static const int nBins = 9;
float binBorders [nBins+1] = { 1, 2, 4.3, 6, 8.68, 10.09, 12.86, 14.18, 16, 19};

TCanvas* c   [2*nBins];
TCanvas* cx1 [2*nBins];
TCanvas* cy1 [2*nBins];
TCanvas* cz1 [2*nBins];



void createEffHist(int tagFlag=1, int xbins=70, int ybins = 70, int zbins = 30)
{
    
    if ( ybins<1 ) ybins = xbins;
    if ( zbins<1 ) zbins = xbins;
    if ( xbins<1 ) return;
    
    RooRealVar  costhetaMC("costhetaMC","cos(#theta_{T})",-1,1);
    RooRealVar cospsiMC("cospsiMC","cos(#psi_{T})",-1,1);
    RooRealVar phiMC("phiMC","#phi",-TMath::Pi(),TMath::Pi());
    RooArgSet vars (costhetaMC, cospsiMC, phiMC);

  //string shortString = Form(tagFlag?"b%ict":"b%iwt",q2Bin);
  //string longString  = Form(tagFlag?"q2 bin %i correct-tag":"q2 bin %i wrong-tag",q2Bin);
  //int confIndex = (tagFlag?q2Bin:q2Bin+nBins);

  // Load ntuples
  TChain* t_den = new TChain();
  TChain* t_num = new TChain();
  t_den->Add("/Users/ab/Documents/Data_MC_Sample/ntuBsGEN.root/OutTree");
  t_num->Add("/Users/ab/Documents/Data_MC_Sample/fittree_ntuBsMC2017.root/treeFit");
  Int_t denEntries = t_den->GetEntries();
  Int_t numEntries = t_num->GetEntries();
  Int_t counter;

 Float_t costhetagen, cospsigen, phigen;
  Float_t recoCosTheta, recoCosPsi, recoPhi;
  Int_t recotag, gentag;
  t_den->SetBranchAddress( "angle_costheta_GEN" , &costhetagen );
  t_den->SetBranchAddress( "angle_cospsi_GEN" , &cospsigen);
  t_den->SetBranchAddress( "angle_phi_GEN", &phigen);
  t_den->SetBranchAddress( "Tag", &gentag);
  t_num->SetBranchAddress( "BscosthetaMC" , &recoCosTheta );
  t_num->SetBranchAddress( "BscospsiMC" , &recoCosPsi );
  t_num->SetBranchAddress( "BsphiMC", &recoPhi );
  t_num->SetBranchAddress( "tag", &recotag);

  RooDataSet* recodata    = new RooDataSet( "recodata"   , "GEN distribution", vars );
  RooDataSet* numData = new RooDataSet( "numData", "RECO distribution after selections", vars ); 

  // Prepare denominator datasets
  cout<<"Starting denominator dataset filling..."<<endl;
  counter=0;
  for (int iCand=0; iCand<denEntries; ++iCand) {
    t_den->GetEntry(iCand);
    
    // status display
    if ( iCand > 1.0*counter*denEntries/100 ) {
      cout<<counter<<"%"<<endl;
      counter += 10;
    }
    // fill
    costhetaMC.setVal(costhetagen);
    cospsiMC.setVal(cospsigen);
    phiMC.setVal(phigen);
    recodata->add(vars);
  }

  // Prepare numerator datasets
  cout<<"Starting numerator dataset filling..."<<endl;
  counter=0;
  for (int iCand=0; iCand<numEntries; ++iCand) {
    t_num->GetEntry(iCand);
   
    // status display
    if ( iCand > 1.0*counter*numEntries/100 ) {
      cout<<counter<<"%"<<endl;
      counter += 10;
    }
    // fill
    costhetaMC.setVal(recoCosTheta);
    cospsiMC.setVal(recoCosPsi);
    phiMC.setVal(recoPhi);
    numData->add(vars);
  }
  cout<<"Dataset prepared"<<endl;

  // compute and print average efficiency
  double avgEff = numData->sumEntries()/recodata->sumEntries();
  cout<<"Average efficiency = "<<avgEff<<endl;

  /*if (plot) {
    // Plot 1D distributions of numerator and denominator datasets
    double rescFac = 0.03;
    auto cnew = new TCanvas("cnew","Num and Den 1D projections ",1200,700);
    TLegend* leg = new TLegend(0.4,0.8,0.9,0.9);
    RooPlot* xframe = costhetaMC.frame(Title(" cos(#theta_{T}) distributions"));
    RooPlot* yframe = cospsiMC.frame(Title("cos(#psi_{T}) distributions"));
    RooPlot* zframe = phiMC.frame(Title("#phi distributions"));
    recodata->plotOn(xframe,Rescale(rescFac),MarkerColor(kRed+1),LineColor(kRed+1),Binning(30),Name("plDenDist"));
    recodata->plotOn(yframe,Rescale(rescFac),MarkerColor(kRed+1),LineColor(kRed+1),Binning(30));
    recodata->plotOn(zframe,Rescale(rescFac),MarkerColor(kRed+1),LineColor(kRed+1),Binning(30));
    numData->plotOn(xframe,Binning(30),Name("plNumDist"));
    numData->plotOn(yframe,Binning(30));
    numData->plotOn(zframe,Binning(30));
    xframe->GetYaxis()->SetTitleOffset(1.6);
    yframe->GetYaxis()->SetTitleOffset(1.6);
    zframe->GetYaxis()->SetTitleOffset(1.6);
    xframe->SetMaximum(xframe->GetMaximum()*rescFac*1.15);
    yframe->SetMaximum(yframe->GetMaximum()*rescFac*1.15);
    zframe->SetMaximum(zframe->GetMaximum()*rescFac*1.15);
    leg->AddEntry(xframe->findObject("plDenDist"),"Generator-level distributions","lep");
    leg->AddEntry(xframe->findObject("plNumDist"),"Post-selection RECO distributions","lep");
    cnew->Divide(3,1);
    cnew->cd(1);
    gPad->SetLeftMargin(0.17); 
    xframe->Draw();
    leg->Draw("same");
    cnew->cd(2);
    gPad->SetLeftMargin(0.17); 
    yframe->Draw();
    leg->Draw("same");
    cnew->cd(3);
    gPad->SetLeftMargin(0.17); 
    zframe->Draw();
    leg->Draw("same");

    cnew->SaveAs("dist_GEN_RECO.pdf" );
  }*/

  // create numerator and denominator histograms
  TH3D* denHist = (TH3D*)recodata->createHistogram( "denHist",
						   costhetaMC,   Binning(xbins,-1,1) ,
						   YVar(cospsiMC,Binning(ybins,-1,1)),
						   ZVar(phiMC,Binning(zbins,-TMath::Pi(),TMath::Pi())) );
  TH3D* numHist = (TH3D*)numData->createHistogram("numHist",
						   costhetaMC,     Binning(xbins,-1,1) ,
						   YVar(cospsiMC,Binning(ybins,-1,1)),
						   ZVar(phiMC,Binning(zbins,-TMath::Pi(),TMath::Pi())) );

  // save histograms
  TFile* fout = new TFile( Form("effHist_%i_%i_%i.root",xbins,ybins,zbins),"RECREATE");
  denHist->Write();
  numHist->Write();
  fout->Write();
  fout->Close();

  if (plot) {
    // plot 1D slices of the efficiency
    vector <TEfficiency*> effHistsX; 
    vector <TEfficiency*> effHistsY;
    vector <TEfficiency*> effHistsZ;
    auto cx1new= new TCanvas("cx1new","Binned efficiency  cos(theta_{T}) slices",1200,1000) ;
    auto cy1new = new TCanvas("cy1new","Binned efficiency cos(psi_{T}) slices",1200,1000) ;
    auto cz1new = new TCanvas("cz1new","Binned efficiency  phi slices" ,1200,1000) ;
    cx1new->Divide(4,6);
    cy1new->Divide(4,6);
    cz1new->Divide(4,6);

    // width of the slices in the hidden variables ("border" is half of it)
    double border = 0.075;

    // loop over slice grid
    for (int i=0; i<5; ++i)
    {
    for (int j=0; j<5; ++j) {

	// central values and borders of the slices in the hidden variables
	double centA = -0.8 + 1.6*i/4;
	double centB = -0.8 + 1.6*j/4;
	double lowA  = TMath::Max( centA - border,  1e-4-1 );
	double lowB  = TMath::Max( centB - border,  1e-4-1 );
	double highA = TMath::Min( centA + border, -1e-4+1 );
	double highB = TMath::Min( centB + border, -1e-4+1 );
    
	// slicing num and den distributions
	auto numProjX = numHist->ProjectionX("numProjX", 
					     numHist->GetYaxis()->FindBin(lowA), numHist->GetYaxis()->FindBin(highA),
					     numHist->GetZaxis()->FindBin(lowB*TMath::Pi()), numHist->GetZaxis()->FindBin(highB*TMath::Pi()),"e");
	auto numProjY = numHist->ProjectionY("numProjY", 
					     numHist->GetXaxis()->FindBin(lowA), numHist->GetXaxis()->FindBin(highA),
					     numHist->GetZaxis()->FindBin(lowB*TMath::Pi()), numHist->GetZaxis()->FindBin(highB*TMath::Pi()),"e");
	auto numProjZ = numHist->ProjectionZ("numProjZ", 
					     numHist->GetXaxis()->FindBin(lowA), numHist->GetXaxis()->FindBin(highA ),
					     numHist->GetYaxis()->FindBin(lowB), numHist->GetYaxis()->FindBin(highB),"e");
	auto denProjX = denHist->ProjectionX("denProjX", 
					     denHist->GetYaxis()->FindBin(lowA), denHist->GetYaxis()->FindBin(highA),
					     denHist->GetZaxis()->FindBin(lowB*TMath::Pi()), denHist->GetZaxis()->FindBin(highB*TMath::Pi()),"e");
	auto denProjY = denHist->ProjectionY("denProjY", 
					     denHist->GetXaxis()->FindBin(lowA), denHist->GetXaxis()->FindBin(highA),
					     denHist->GetZaxis()->FindBin(lowB*TMath::Pi()), denHist->GetZaxis()->FindBin(highB*TMath::Pi()),"e");
	auto denProjZ = denHist->ProjectionZ("denProjZ", 
					     denHist->GetXaxis()->FindBin(lowA), denHist->GetXaxis()->FindBin(highA),
					     denHist->GetYaxis()->FindBin(lowB), denHist->GetYaxis()->FindBin(highB),"e");

	// producing 1D efficiencies from the slices
	effHistsX.push_back( new TEfficiency(*numProjX,*denProjX) );
	effHistsX.back()->SetName( Form("effHistX_%i_%i",i,j) );
	effHistsX.back()->SetTitle( Form("Efficiency slice cospsi=%1.2f phi=%1.2f;cos(#theta_{T});Efficiency", centA, centB*TMath::Pi()));
    
	effHistsY.push_back( new TEfficiency(*numProjY,*denProjY) );
	effHistsY.back()->SetName( Form("effHistY_%i_%i",i,j) );
	effHistsY.back()->SetTitle( Form("Efficiency slice costheta=%1.2f phi=%1.2f;cos(#psi_{T});Efficiency", centA, centB*TMath::Pi()));

	effHistsZ.push_back( new TEfficiency(*numProjZ,*denProjZ) );
	effHistsZ.back()->SetName( Form("effHistZ_%i_%i",i,j) );
	effHistsZ.back()->SetTitle( Form("Efficiency slice costheta=%1.2f cospsi=%1.2f;#phi;Efficiency",centA,centB));

	// plot in canvas
	cx1new->cd(5*j+i+1);
	gPad->SetLeftMargin(0.18); 
	effHistsX.back()->Draw();
	cx1new->cd(5*j+i+1)->Update(); 
	auto graphx = effHistsX.back()->GetPaintedGraph(); 
	graphx->SetMinimum(0);
	// graphx->SetMaximum(0.2);
	graphx->GetYaxis()->SetTitleOffset(1.7);
	cx1new->cd(5*j+i+1)->Update();

	cy1new->cd(5*j+i+1);
	gPad->SetLeftMargin(0.18); 
	effHistsY.back()->Draw();
	cy1new->cd(5*j+i+1)->Update(); 
	auto graphy = effHistsY.back()->GetPaintedGraph(); 
	graphy->SetMinimum(0);
	// graphy->SetMaximum(0.2); 
	graphy->GetYaxis()->SetTitleOffset(1.7);
	cy1new->cd(5*j+i+1)->Update();

	cz1new->cd(5*j+i+1);
	gPad->SetLeftMargin(0.18); 
	effHistsZ.back()->Draw();
	cz1new->cd(5*j+i+1)->Update(); 
	auto graphz = effHistsZ.back()->GetPaintedGraph(); 
	graphz->SetMinimum(0);
	// graphz->SetMaximum(0.2); 
	graphz->GetYaxis()->SetTitleOffset(1.7);
	cz1new->cd(5*j+i+1)->Update();

      }}    

    cx1new->SaveAs( Form("effHist_%i_%i_%i_CosTheta_slices_dp%i.pdf",xbins,ybins,zbins,(int)(border*200)));
    cy1new->SaveAs( Form("effHist_%i_%i_%i_CosPsi_slices_dp%i.pdf",xbins,ybins,zbins,(int)(border*200)));
    cz1new->SaveAs( Form("effHist_%i_%i_%i_Phi_slices_dp%i.pdf",xbins,ybins,zbins,(int)(border*200)));
  }
    
}
