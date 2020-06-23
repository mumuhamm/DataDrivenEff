#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooWorkspace.h"
#include "RooCategory.h"
#include "TKey.h"
#include "RooExponential.h"
#include <map>
#include "TCut.h"
#include "RooHistPdf.h"
#include "RooHist.h"
#include "TVirtualPad.h"
#include "RooDataHist.h"
#include <string>
#include "TEventList.h"
#include "RooFit.h"
#include "RooRealVar.h"
#include "TFile.h"
#include "RooDataSet.h"
#include "TTree.h"
#include "TH2D.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooGaussian.h"
#include "RooGenericPdf.h"
#include "RooFitResult.h"
#include "RooAddPdf.h"
#include "RooAddition.h"
#include "RooAbsPdf.h"
#include "RooPolynomial.h"
#include "RooGaussModel.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooDecay.h"
#include "RooDataHist.h"
#include "RooAddModel.h"
#include "RooProdPdf.h"
#include "RooProduct.h"
#include "RooPlot.h"
#include "TH1D.h"
#include "TRandom.h"
#include "RooMinuit.h"
#include "RooExtendPdf.h"
#include "RooChi2Var.h"
#include "Math/Functor.h"
#include "TRandom3.h"
#include "Math/DistFunc.h"
#include "RooClassFactory.h"
#include "RooFitResult.h"
#include "RooDataSet.h"
#include "RooRealConstant.h"
#include "RooConstVar.h"
#include "Roo1DTable.h"
#include "RooBDecay.h"
#include "RooFormulaVar.h"
#include "RooTruthModel.h"
#include "RooRealSumPdf.h"
#include "Math/SpecFunc.h"
#include "RooBMixDecay.h"
#include "RooBCPEffDecay.h"
#include "Riostream.h"
#include "RooRandom.h"
#include "TMath.h"
#include "RooFun1TRPdf.h"
#include "RooFun2TRPdf.h"
#include "RooFun3TRPdf.h"
#include "RooFun4TRPdf.h"
#include "RooFun5TRPdf.h"
#include "RooFun6TRPdf.h"
#include "RooEffcthPdf.h"
#include "RooMCStudy.h"
#include "RooArgSet.h"
#include "RooLegendre.h"
#include "RooSpHarmonic.h"
#include "RooBifurGauss.h"

using namespace RooFit; 

void modelbuilderMC2(string target) {
  
    //gROOT->SetStyle("Plain");
    
    //#####################################################################
    //                                                                    #
    //   Definition of the PDF parameters for Signal and Background       #
    //                                                                    #
    //#####################################################################
    
    
      //                    Parameters non fixed
    //############################################################################
    
    // Declare observables  (M(Bs)[BdFitM],Ct(Bs)[BdCt2DBS],Bdcostheta,Bdcospsi,Bdphi)
    //___________________________________________________________________________
    RooRealVar *BsCt2DMC= new RooRealVar ("BsCt2DMC","Bs ct", 0.,0.3,"cm");//
    
    
    RooRealVar *BscosthetaMC= new RooRealVar ("BscosthetaMC","cos(#theta_{T})", -1,1);
    RooRealVar *BscospsiMC= new RooRealVar ("BscospsiMC","cos(#psi_{T})", -1,1);
    RooRealVar *BsphiMC= new RooRealVar ("BsphiMC","#phi_{T}", -TMath::Pi(),TMath::Pi(),"rad");//cosdelta2
    
    
    //Physical parameters  included in the signal PDF
    // --------------------------------------------------
    //Time dependent functions, Amplitudes |A_0|,|A_perp|,|A_par|
    //Time Amplitudes variables |A_0|,|A_perp|,|A_par|
    
/*    A_0=new RooRealVar ("A_0","|A_0|^2|",0.53,0.1,0.9);
    A_1=new RooRealVar ("A_1","|A_1|^2",0.34,0.,1.); 
    A_2=new RooFormulaVar ("A_2","|A_2|*2","(1-@0)*@1",RooArgList(*A_0,*A_1)); 
    A_3=new RooFormulaVar ("A_3","|A_3|*2","(1-@0)*(1-@1)",RooArgList(*A_0,*A_1)); 
    A_4=new RooFormulaVar ("A_4","A_4","sqrt(@0)*sqrt(@1)",RooArgList(*A_3,*A_2));
    A_5=new RooFormulaVar ("A_5","|A_0||A_pa|","sqrt(@0)*sqrt(@1)",RooArgList(*A_0,*A_3));
    A_6=new RooFormulaVar ("A_6","A_6","sqrt(@0)*sqrt(@1)",RooArgList(*A_0,*A_2));*/
   
    RooRealVar *A_0=new RooRealVar ("A_0","|A_0|^2|",0.51,0.30,0.7);
    RooRealVar *A_S=new RooRealVar ("A_S","|A_0|^2|",0.05,0.,0.1);
    RooRealVar *A_pe=new RooRealVar ("A_pe","|A_pe|^2",0.24,0.1,0.6);
    RooFormulaVar *A_pa=new RooFormulaVar ("A_pa","|A_pa|*2","1-@0-@1",RooArgList(*A_0,*A_pe));
    //RooRealVar *A_pa=new RooRealVar ("A_pa","|A_pa|^2",0.24,0.2,0.3);
    //RooFormulaVar *A_pe=new RooFormulaVar ("A_pe","|A_pe|*2","1-@0-@1",RooArgList(*A_0,*A_pa));
    RooFormulaVar *A_4=new RooFormulaVar ("A_4","A_4","sqrt(@0)*sqrt(@1)",RooArgList(*A_pa,*A_pe));
    RooFormulaVar *A_5=new RooFormulaVar ("A_5","|A_0||A_pa|","sqrt(@0)*sqrt(@1)",RooArgList(*A_0,*A_pa));
    RooFormulaVar *A_6=new RooFormulaVar ("A_6","A_6","sqrt(@0)*sqrt(@1)",RooArgList(*A_0,*A_pe));
    RooFormulaVar *A_8=new RooFormulaVar ("A_8","A_8","sqrt(@0)*sqrt(@1)",RooArgList(*A_S,*A_pa));
    RooFormulaVar *A_9=new RooFormulaVar ("A_9","|A_9","sqrt(@0)*sqrt(@1)",RooArgList(*A_S,*A_pe));
    RooFormulaVar *A_10=new RooFormulaVar ("A_10","A_10","sqrt(@0)*sqrt(@1)",RooArgList(*A_S,*A_0)); 
    //physical parameters, Phi_s, delta_1, delta_2, dm=Mass_L-Mass_H, DGamma=(Gamma_L-Gamma_H)/2
    
    RooRealVar *phi_s=new RooRealVar ("phi_s","phi_s",-0.04,-1,1);//-0.04,-1.,0.);//-0.04
    RooRealVar *deltaPa=new RooRealVar("deltaPa","#deltaPa",3.4,0.,4);//,2.5);   
    RooRealVar *deltaPe=new RooRealVar("deltaPe","#deltaPe",2.3,0.,4);//,-1,4); // =0.17  
    RooRealVar *deltaSPe=new RooRealVar("deltaSPe","#deltaPe",-0.028,-3.,3.);//,-1,4); // =0.17  
    RooRealVar *dm=new RooRealVar ("dm","dm",589.7,500,700); //590.1// in um*c, PDG value 17.69 in 1/ps
    RooRealVar *dGam=new RooRealVar ("dGam","dGam",3.21,1.,6.);//2.3//original value in MC
    RooRealVar *tau=new RooRealVar ("tau","tau",0.0447,0.03,0.055,"ps");//,0.043//origunal value in MC

    double tagmis = 1;
    RooRealVar *mistag = new RooRealVar("mistag","Mistag fraction of original B and Bbar",tagmis);
    RooCategory *tag = new RooCategory("tag","Flavour tag of the B meson");
    tag->defineType("Bsbar",1);
    tag->defineType("Bs",-1);
    tag->defineType("untag",0);
    
    
 //   RooDataSet *Alldataset= readDataset(AllFilenames, "Alldataset");
   
 //    int nentries = Alldataset->sumEntries();
 //   std::cout<< nentries << std::endl;

    RooTruthModel truth("truth","truth",*BsCt2DMC); 
/*
    RooFormulaVar fsinh1("fsinh1","fsinh1","-cos(@0)",RooArgList(*phi_s));
    RooFormulaVar fsin1("fsin1","fsin1","-@1*(1-2*@2)*sin(@0)",RooArgList(*phi_s,*tag,*mistag));
    RooBDecay *Amp0=new  RooBDecay("Amp0","Amp0",*BsCt2DMC,*tau,*dGam,RooConst(1),fsinh1,RooConst(0),fsin1,*dm,truth,RooBDecay::SingleSided);

    RooFormulaVar fsinh2("fsinh2","fsinh2","cos(@0)",RooArgList(*phi_s));
    RooFormulaVar fsin2("fsin2","fsin2","@1*(1-2*@2)*sin(@0)",RooArgList(*phi_s,*tag,*mistag));
    RooBDecay *Ampe=new  RooBDecay("Ampe","Ampe",*BsCt2DMC,*tau,*dGam,RooConst(1),fsinh2,RooConst(0),fsin2,*dm,truth,RooBDecay::SingleSided);


    RooFormulaVar fsinh4("fsinh4","fsinh4","-sin(@0)*cos(@1-@2)",RooArgList(*phi_s,*deltaPe,*deltaPa));
    RooFormulaVar fsin4("fsin4","fsin4","@1*(1-2*@2)*cos(@0)*cos(@3-@4)",RooArgList(*phi_s,*tag,*mistag,*deltaPe,*deltaPa));
    RooFormulaVar fcos4("fcos4","fcos4","-@0*(1-2*@1)*sin(@2-@3)",RooArgList(*tag,*mistag,*deltaPe,*deltaPa));
    RooBDecay *Ampa4=new RooBDecay("Ampa4","Ampa4",*BsCt2DMC,*tau,*dGam,RooRealConstant::value(0),fsinh4,fcos4,fsin4,*dm,truth,RooBDecay::SingleSided);

    RooFormulaVar fsinh6("fsinh6","fsinh6","-sin(@0)*cos(@1)",RooArgList(*phi_s,*deltaPe));
    RooFormulaVar fsin6("fsin6","fsin6","@1*(1-2*@2)*cos(@0)*cos(@3)",RooArgList(*phi_s,*tag,*mistag,*deltaPe));
    RooFormulaVar fcos6("fcos6","fcos6","-@0*(1-2*@1)*sin(@2)",RooArgList(*tag,*mistag,*deltaPe));
    RooBDecay *Ampa6=new RooBDecay("Ampa6","Ampa6",*BsCt2DMC,*tau,*dGam,RooRealConstant::value(0),fsinh6,fcos6,fsin6,*dm,truth,RooBDecay::SingleSided);


    RooFormulaVar coshGBasis("coshGBasis","exp(-@0/@1)*cosh(@0*@2/2)", RooArgList(*BsCt2DMC,*tau,*dGam));
    RooFormulaVar sinhGBasis("sinhGBasis","exp(-@0/@1)*sinh(@0*@2/2)", RooArgList(*BsCt2DMC,*tau,*dGam));
    RooFormulaVar cosGBasis("cosGBasis","exp(-@0/@1)*cos(@0*@2)",RooArgList(*BsCt2DMC,*tau,*dm));
    RooFormulaVar sinGBasis("sinGBasis","exp(-@0/@1)*sin(@0*@2)",RooArgList(*BsCt2DMC,*tau,*dm));
    RooAbsReal* coshGConv = truth.convolution(&coshGBasis,BsCt2DMC);
    RooAbsReal* sinhGConv = truth.convolution(&sinhGBasis,BsCt2DMC);
    RooAbsReal* cosGConv = truth.convolution(&cosGBasis,BsCt2DMC);
    RooAbsReal* sinGConv = truth.convolution(&sinGBasis,BsCt2DMC);
    RooAddition *myAmp0=new RooAddition("myAmp0","myAmp0",RooArgList(*sinhGConv,*coshGConv,*sinGConv),RooArgList(fsinh1,RooConst(1),fsin1));
    RooAddition *myAmpe=new RooAddition("myAmpe","myAmpe",RooArgList(*sinhGConv,*coshGConv,*sinGConv),RooArgList(fsinh2,RooConst(1),fsin2));
    RooAddition *myAmpa4=new RooAddition("myAmpa4","myAmpa4",RooArgList(*sinhGConv,*cosGConv,*sinGConv),RooArgList(fsinh4,fcos4,fsin4));
    RooAddition *myAmpa6=new RooAddition("myAmpa6","myAmpa6",RooArgList(*sinhGConv,*cosGConv,*sinGConv),RooArgList(fsinh6,fcos6,fsin6));
*/
/*    RooFormulaVar fsinh1("fsinh1","fsinh1","-cos(@0)",RooArgList(*phi_s));
    RooFormulaVar fsin1bar("fsin1bar","fsin1bar","-sin(@0)",RooArgList(*phi_s));
    RooFormulaVar fsin1("fsin1","fsin1","sin(@0)",RooArgList(*phi_s));
    RooBDecay *Amp0=new  RooBDecay("Amp0","Amp0",*BsCt2DMC,*tau,*dGam,RooConst(1),fsinh1,RooConst(0),fsin1,*dm,truth,RooBDecay::SingleSided); 

    RooFormulaVar fsinh2("fsinh2","fsinh2","cos(@0)",RooArgList(*phi_s));
    RooFormulaVar fsin2bar("fsin2bar","fsin2bar","sin(@0)",RooArgList(*phi_s));
    RooFormulaVar fsin2("fsin2","fsin2","-sin(@0)",RooArgList(*phi_s));
    RooBDecay *Ampe=new  RooBDecay("Ampe","Ampe",*BsCt2DMC,*tau,*dGam,RooConst(1),fsinh2,RooConst(0),fsin2,*dm,truth,RooBDecay::SingleSided);


    RooFormulaVar fsinh4("fsinh4","fsinh4","-sin(@0)*cos(@1-@2)",RooArgList(*phi_s,*deltaPe,*deltaPa));
    RooFormulaVar fsin4bar("fsin4bar","fsin4bar","cos(@0)*cos(@1-@2)",RooArgList(*phi_s,*deltaPe,*deltaPa));
    RooFormulaVar fcos4bar("fcos4bar","fcos4bar","-sin(@0-@1)",RooArgList(*deltaPe,*deltaPa));
    RooFormulaVar fsin4("fsin4","fsin4","-cos(@0)*cos(@1-@2)",RooArgList(*phi_s,*deltaPe,*deltaPa));
    RooFormulaVar fcos4("fcos4","fcos4","sin(@0-@1)",RooArgList(*deltaPe,*deltaPa));
    RooBDecay *Ampa4=new RooBDecay("Ampa4","Ampa4",*BsCt2DMC,*tau,*dGam,RooRealConstant::value(0),fsinh4,fcos4,fsin4,*dm,truth,RooBDecay::SingleSided);    

    RooFormulaVar fsinh6("fsinh6","fsinh6","-sin(@0)*cos(@1)",RooArgList(*phi_s,*deltaPe));
    RooFormulaVar fsin6bar("fsin6bar","fsin6bar","cos(@0)*cos(@1)",RooArgList(*phi_s,*deltaPe));
    RooFormulaVar fcos6bar("fcos6bar","fcos6bar","-sin(@0)",RooArgList(*deltaPe));
    RooFormulaVar fsin6("fsin6","fsin6","-cos(@0)*cos(@1)",RooArgList(*phi_s,*deltaPe));
    RooFormulaVar fcos6("fcos6","fcos6","sin(@0)",RooArgList(*deltaPe));
    RooBDecay *Ampa6=new RooBDecay("Ampa6","Ampa6",*BsCt2DMC,*tau,*dGam,RooRealConstant::value(0),fsinh6,fcos6,fsin6,*dm,truth,RooBDecay::SingleSided);    
*/
    RooFormulaVar fsinh1("fsinh1","fsinh1","-cos(@0)",RooArgList(*phi_s));
    RooFormulaVar fsin1bar("fsin1bar","fsin1bar","-sin(@0)",RooArgList(*phi_s));
    RooFormulaVar fsin1("fsin1","fsin1","sin(@0)",RooArgList(*phi_s));
    RooFormulaVar fsin1tag("fsin1tag","fsin1tag","-(1-2*@2)*@1*sin(@0)",RooArgList(*phi_s,*tag,*mistag));
    
    RooFormulaVar fsinh2("fsinh2","fsinh2","cos(@0)",RooArgList(*phi_s));
    RooFormulaVar fsin2bar("fsin2bar","fsin2bar","sin(@0)",RooArgList(*phi_s));
    RooFormulaVar fsin2("fsin2","fsin2","-sin(@0)",RooArgList(*phi_s));
    RooFormulaVar fsin2tag("fsin2tag","fsin2tag","(1-2*@2)*@1*sin(@0)",RooArgList(*phi_s,*tag,*mistag));

    RooFormulaVar fsinh4("fsinh4","fsinh4","-sin(@0)*cos(@1-@2)/sin(@1-@2)",RooArgList(*phi_s,*deltaPe,*deltaPa));
    RooFormulaVar fsin4bar("fsin4bar","fsin4bar","cos(@0)*cos(@1-@2)",RooArgList(*phi_s,*deltaPe,*deltaPa));
    RooFormulaVar fcos4bar("fcos4bar","fcos4bar","-sin(@0-@1)",RooArgList(*deltaPe,*deltaPa));
    RooFormulaVar fsin4("fsin4","fsin4","-cos(@0)*cos(@1-@2)",RooArgList(*phi_s,*deltaPe,*deltaPa));
    RooFormulaVar fcos4("fcos4","fcos4","sin(@0-@1)",RooArgList(*deltaPe,*deltaPa));
    RooFormulaVar fsin4tag("fsin4tag","fsin4tag","(1-2*@4)*@3*cos(@0)*cos(@1-@2)/sin(@1-@2)",RooArgList(*phi_s,*deltaPe,*deltaPa,*tag,*mistag));
    RooFormulaVar fcos4tag("fcos4tag","fcos4tag","-(1-2*@1)*@0",RooArgList(*tag,*mistag));
    
    RooFormulaVar fsinh5("fsinh5","fsinh5","-cos(@0)",RooArgList(*phi_s));
    RooFormulaVar fcosh5("fcosh5","fcosh5","cos(@0)",RooArgList(*deltaPa));
    RooFormulaVar fsin5bar("fsin5bar","fsin5bar","-sin(@0)*cos(@1)",RooArgList(*phi_s,*deltaPa));
    RooFormulaVar fsin5("fsin5","fsin5","sin(@0)*cos(@1)",RooArgList(*phi_s,*deltaPa));
    RooFormulaVar fsin5tag("fsin5tag","fsin5tag","-(1-2*@2)*@1*sin(@0)",RooArgList(*phi_s,*tag,*mistag));

    RooFormulaVar fsinh6("fsinh6","fsinh6","-sin(@0)*cos(@1)/sin(@1)",RooArgList(*phi_s,*deltaPe));
    RooFormulaVar fsin6bar("fsin6bar","fsin6bar","cos(@0)*cos(@1)",RooArgList(*phi_s,*deltaPe));
    RooFormulaVar fcos6bar("fcos6bar","fcos6bar","-sin(@0)",RooArgList(*deltaPe));
    RooFormulaVar fsin6("fsin6","fsin6","-cos(@0)*cos(@1)",RooArgList(*phi_s,*deltaPe));
    RooFormulaVar fcos6("fcos6","fcos6","sin(@0)",RooArgList(*deltaPe));
    RooFormulaVar fsin6tag("fsin6tag","fsin6tag","(1-2*@3)*@2*cos(@0)*cos(@1)/sin(@1)",RooArgList(*phi_s,*deltaPe,*tag,*mistag));
    RooFormulaVar fcos6tag("fcos6tag","fcos6tag","-(1-2*@1)*@0",RooArgList(*tag,*mistag));

    RooFormulaVar coshGBasis("coshGBasis","exp(-@0/@1)*cosh(@0*@2/2)", RooArgList(*BsCt2DMC,*tau,*dGam));
    RooFormulaVar sinhGBasis("sinhGBasis","exp(-@0/@1)*sinh(@0*@2/2)", RooArgList(*BsCt2DMC,*tau,*dGam));
    RooFormulaVar cosGBasis("cosGBasis","exp(-@0/@1)*cos(@0*@2)",RooArgList(*BsCt2DMC,*tau,*dm));
    RooFormulaVar sinGBasis("sinGBasis","exp(-@0/@1)*sin(@0*@2)",RooArgList(*BsCt2DMC,*tau,*dm));
    RooAbsReal* coshGConv = truth.convolution(&coshGBasis,BsCt2DMC);
    RooAbsReal* sinhGConv = truth.convolution(&sinhGBasis,BsCt2DMC);
    RooAbsReal* cosGConv = truth.convolution(&cosGBasis,BsCt2DMC);
    RooAbsReal* sinGConv = truth.convolution(&sinGBasis,BsCt2DMC);
    RooAddition *myAmp0bar=new RooAddition("myAmp0bar","myAmp0bar",RooArgList(*sinhGConv,*coshGConv,*sinGConv),RooArgList(fsinh1,RooConst(1),fsin1bar));
    RooAddition *myAmpebar=new RooAddition("myAmpebar","myAmpebar",RooArgList(*sinhGConv,*coshGConv,*sinGConv),RooArgList(fsinh2,RooConst(1),fsin2bar));
    RooAddition *myAmpa4bar=new RooAddition("myAmpa4bar","myAmpa4bar",RooArgList(*sinhGConv,*cosGConv,*sinGConv),RooArgList(fsinh4,fcos4bar,fsin4bar));
    RooAddition *myAmpa6bar=new RooAddition("myAmpa6bar","myAmpa6bar",RooArgList(*sinhGConv,*cosGConv,*sinGConv),RooArgList(fsinh6,fcos6bar,fsin6bar));
    RooAddition *myAmp0=new RooAddition("myAmp0","myAmp0",RooArgList(*sinhGConv,*coshGConv,*sinGConv),RooArgList(fsinh1,RooConst(1),fsin1));
    RooAddition *myAmpe=new RooAddition("myAmpe","myAmpe",RooArgList(*sinhGConv,*coshGConv,*sinGConv),RooArgList(fsinh2,RooConst(1),fsin2));
    RooAddition *myAmpa4=new RooAddition("myAmpa4","myAmpa4",RooArgList(*sinhGConv,*cosGConv,*sinGConv),RooArgList(fsinh4,fcos4,fsin4));
    RooAddition *myAmpa6=new RooAddition("myAmpa6","myAmpa6",RooArgList(*sinhGConv,*cosGConv,*sinGConv),RooArgList(fsinh6,fcos6,fsin6));
    RooAddition *myAmp0tag=new RooAddition("myAmp0tag","myAmp0tag",RooArgList(*sinhGConv,*coshGConv,*sinGConv),RooArgList(fsinh1,RooConst(1),fsin1tag));
    RooAddition *myAmpetag=new RooAddition("myAmpetag","myAmpetag",RooArgList(*sinhGConv,*coshGConv,*sinGConv),RooArgList(fsinh2,RooConst(1),fsin2tag));
    RooAddition *myAmpa4tag=new RooAddition("myAmpa4tag","myAmpa4tag",RooArgList(*sinhGConv,*cosGConv,*sinGConv),RooArgList(fsinh4,fcos4tag,fsin4tag));
    RooAddition *myAmpa5tag=new RooAddition("myAmpa5tag","myAmpa5tag",RooArgList(*sinhGConv,*coshGConv,*sinGConv),RooArgList(fsinh5,RooConst(1),fsin5tag));
    RooAddition *myAmpa6tag=new RooAddition("myAmpa6tag","myAmpa6tag",RooArgList(*sinhGConv,*cosGConv,*sinGConv),RooArgList(fsinh6,fcos6tag,fsin6tag));
    RooAddition *myAmp0untag=new RooAddition("myAmp0untag","myAmp0untag",RooArgList(*sinhGConv,*coshGConv),RooArgList(fsinh1,RooConst(1)));
    RooAddition *myAmpeuntag=new RooAddition("myAmpeunatag","myAmpeuntag",RooArgList(*sinhGConv,*coshGConv),RooArgList(fsinh2,RooConst(1)));
    RooAddition *myAmpa4untag=new RooAddition("myAmpa4untag","myAmpa4untag",RooArgList(*sinhGConv),RooArgList(fsinh4));
    RooAddition *myAmpa5untag=new RooAddition("myAmpa5untag","myAmpa5untag",RooArgList(*sinhGConv,*coshGConv),RooArgList(fsinh5,fcosh5));
    RooAddition *myAmpa6untag=new RooAddition("myAmpa6untag","myAmpa6untag",RooArgList(*sinhGConv),RooArgList(fsinh6));
//RooBDecay *myAmp0untag=new  RooBDecay("Amp0","Amp0",*BsCt2DMC,*tau,*dGam,RooConst(1),fsinh1,RooConst(0),fsin1,*dm,truth,RooBDecay::SingleSided);
//RooBDecay *myAmpeuntag=new  RooBDecay("Ampe","Ampe",*BsCt2DMC,*tau,*dGam,RooConst(1),fsinh2,RooConst(0),fsin2,*dm,truth,RooBDecay::SingleSided);
//RooBDecay *myAmpa4untag=new RooBDecay("Ampa4","Ampa4",*BsCt2DMC,*tau,*dGam,RooConst(0),fsinh4,fcos4,fsin4,*dm,truth,RooBDecay::SingleSided);   
//RooBDecay *myAmpa5untag=new RooBDecay("Ampa5","Ampa5",*BsCt2DMC,*tau,*dGam,fcosh5,fsinh5,RooConst(0),fsin5,*dm,truth,RooBDecay::SingleSided);   
//RooBDecay *myAmpa6untag=new RooBDecay("Ampa6","Ampa6",*BsCt2DMC,*tau,*dGam,RooConst(0),fsinh6,fcos6,fsin6,*dm,truth,RooBDecay::SingleSided);    
//RooBDecay *myAmp0untag=new  RooBDecay("Amp0","Amp0",*BsCt2DMC,*tau,*dGam,RooConst(1),fsinh1,RooConst(0),RooConst(0),*dm,truth,RooBDecay::SingleSided);
//RooBDecay *myAmpeuntag=new  RooBDecay("Ampe","Ampe",*BsCt2DMC,*tau,*dGam,RooConst(1),fsinh2,RooConst(0),RooConst(0),*dm,truth,RooBDecay::SingleSided);
//RooBDecay *myAmpa4untag=new RooBDecay("Ampa4","Ampa4",*BsCt2DMC,*tau,*dGam,RooConst(0),fsinh4,RooConst(0),RooConst(0),*dm,truth,RooBDecay::SingleSided);   
//RooBDecay *myAmpa5untag=new RooBDecay("Ampa5","Ampa5",*BsCt2DMC,*tau,*dGam,fcosh5,fsinh5,RooConst(0),RooConst(0),*dm,truth,RooBDecay::SingleSided);   
//RooBDecay *myAmpa5untag=new RooBDecay("Ampa5","Ampa5",*BsCt2DMC,*tau,*dGam,RooConst(1),fsinh5bis,RooConst(0),RooConst(0),*dm,truth,RooBDecay::SingleSided);   
//RooBDecay *myAmpa6untag=new RooBDecay("Ampa6","Ampa6",*BsCt2DMC,*tau,*dGam,RooConst(0),fsinh6,RooConst(0),RooConst(0),*dm,truth,RooBDecay::SingleSided);    
    RooFormulaVar *cosdpa=new RooFormulaVar("cosdpa","cosdpa","cos(@0)",RooArgList(*deltaPa));
    RooFormulaVar *sindpe=new RooFormulaVar("sindpe","sindpe","sin(@0)",RooArgList(*deltaPe));
    RooFormulaVar *sindpadpe=new RooFormulaVar("sindpadpe","sindpadpe","sin(@0-@1)",RooArgList(*deltaPe,*deltaPa));
    //Angular dependent functions funiT
    
/*    RooFun1TRPdf* fun1T=new RooFun1TRPdf("fun1T","fun1T",*BscospsiMC,*BscosthetaMC,*BsphiMC);
    RooFun2TRPdf* fun2T=new RooFun2TRPdf("fun2T","fun2T",*BscospsiMC,*BscosthetaMC,*BsphiMC);
    RooFun3TRPdf* fun3T=new RooFun3TRPdf("fun3T","fun3T",*BscospsiMC,*BscosthetaMC);
    RooFun4TRPdf* fun4T=new RooFun4TRPdf("fun4T","fun4T",*BscospsiMC,*BscosthetaMC,*BsphiMC);
    RooFun5TRPdf* fun5T=new RooFun5TRPdf("fun5T","fun5T",*BscospsiMC,*BscosthetaMC,*BsphiMC);
    RooFun6TRPdf* fun6T=new RooFun6TRPdf("fun6T","fun6T",*BscospsiMC,*BscosthetaMC,*BsphiMC);*/

    RooFormulaVar phinew("phinew","BsphiMC/TMath::Pi()",*BsphiMC);
    RooLegendre phizero("phizero","phizero",phinew,0,0);
    RooFormulaVar* fun1T = new RooFormulaVar("fun1T","2*BscospsiMC*BscospsiMC*(1-(1-BscosthetaMC*BscosthetaMC)*cos(BsphiMC)*cos(BsphiMC))",RooArgSet(*BscosthetaMC,*BsphiMC,*BscospsiMC));
    RooFormulaVar* fun2T = new RooFormulaVar("fun2T","(1-BscospsiMC*BscospsiMC)*(1-(1-BscosthetaMC*BscosthetaMC)*sin(BsphiMC)*sin(BsphiMC))",RooArgSet(*BscosthetaMC,*BsphiMC,*BscospsiMC));
    RooFormulaVar* fun3T = new RooFormulaVar("fun3T","(1-BscospsiMC*BscospsiMC)*(1-BscosthetaMC*BscosthetaMC)*phizero",RooArgSet(*BscosthetaMC,*BscospsiMC,phizero));
    RooFormulaVar* fun4T = new RooFormulaVar("fun4T","-(1-BscospsiMC*BscospsiMC)*2*BscosthetaMC*sqrt(1-BscosthetaMC*BscosthetaMC)*sin(BsphiMC)",RooArgSet(*BscosthetaMC,*BsphiMC,*BscospsiMC));
    RooFormulaVar* fun5T = new RooFormulaVar("fun5T","2/sqrt(2.)*BscospsiMC*sqrt(1-BscospsiMC*BscospsiMC)*(1-BscosthetaMC*BscosthetaMC)*sin(2*BsphiMC)",RooArgSet(*BscosthetaMC,*BsphiMC,*BscospsiMC));
    RooFormulaVar* fun6T = new RooFormulaVar("fun6T","2/sqrt(2.)*BscospsiMC*sqrt(1-BscospsiMC*BscospsiMC)*2*BscosthetaMC*sqrt(1-BscosthetaMC*BscosthetaMC)*cos(BsphiMC)",RooArgSet(*BscosthetaMC,*BsphiMC,*BscospsiMC));
    RooFormulaVar* fun7T = new RooFormulaVar("fun7T","2/3*(1-(1-@0*@0)*cos(@1)*cos(@1))",RooArgSet(*BscosthetaMC,*BsphiMC));
    RooFormulaVar* fun8T = new RooFormulaVar("fun8T","sqrt(6.)/3*sqrt(1-BscospsiMC*BscospsiMC)*(1-BscosthetaMC*BscosthetaMC)*sin(2*BsphiMC)",RooArgSet(*BscosthetaMC,*BsphiMC,*BscospsiMC));
    RooFormulaVar* fun9T = new RooFormulaVar("fun9T","sqrt(6.)/3*sqrt(1-BscospsiMC*BscospsiMC)*2*BscosthetaMC*sqrt(1-BscosthetaMC*BscosthetaMC)*cos(BsphiMC)",RooArgSet(*BscosthetaMC,*BsphiMC,*BscospsiMC));
    RooFormulaVar* fun10T = new RooFormulaVar("fun10T","sqrt(3.)*4/3*BscospsiMC*(1-(1-BscosthetaMC*BscosthetaMC)*cos(BsphiMC)*cos(BsphiMC))",RooArgSet(*BscosthetaMC,*BsphiMC,*BscospsiMC));

    //RooFormulaVar *coef11=new RooFormulaVar("coef11","coef11","@0",RooArgList(*A_0));
    RooFormulaVar *coef12=new RooFormulaVar("coef12","coef12","-cos(@0)",RooArgList(*phi_s));
    RooFormulaVar *coef14=new RooFormulaVar("coef14","coef14","-sin(@0)*@1*(1-2*@2)",RooArgList(*phi_s,*tag,*mistag));
    //RooFormulaVar *coef21=new RooFormulaVar("coef21","coef21","@0",RooArgList(*A_pa));
    RooFormulaVar *coef22=new RooFormulaVar("coef22","coef22","-cos(@0)",RooArgList(*phi_s));
    RooFormulaVar *coef24=new RooFormulaVar("coef24","coef24","-sin(@0)*@1*(1-2*@2)",RooArgList(*phi_s,*tag,*mistag));
    //RooFormulaVar *coef31=new RooFormulaVar("coef31","coef31","@0",RooArgList(*A_pe));
    RooFormulaVar *coef32=new RooFormulaVar("coef32","coef32","cos(@0)",RooArgList(*phi_s));
    RooFormulaVar *coef34=new RooFormulaVar("coef34","coef34","sin(@0)*@1*(1-2*@2)",RooArgList(*phi_s,*tag,*mistag));
    RooFormulaVar *coef42=new RooFormulaVar("coef42","coef42","-sin(@0)*cos(@1-@2)",RooArgList(*phi_s,*deltaPe,*deltaPa));
    RooFormulaVar *coef43=new RooFormulaVar("coef43","coef43","-sin(@0-@1)*@2*(1-2*@3)",RooArgList(*deltaPe,*deltaPa,*tag,*mistag));
    RooFormulaVar *coef44=new RooFormulaVar("coef44","coef44","cos(@0)*cos(@1-@2)*@3*(1-2*@4)",RooArgList(*phi_s,*deltaPe,*deltaPa,*tag,*mistag));
    RooFormulaVar *coef51=new RooFormulaVar("coef51","coef51","cos(@0)",RooArgList(*deltaPa));
    RooFormulaVar *coef52=new RooFormulaVar("coef52","coef52","-cos(@0)*cos(@1)",RooArgList(*deltaPa,*phi_s));
    RooFormulaVar *coef54=new RooFormulaVar("coef54","coef54","-cos(@0)*sin(@1)*@2*(1-2*@3)",RooArgList(*deltaPa,*phi_s,*tag,*mistag));
    RooFormulaVar *coef62=new RooFormulaVar("coef62","coef62","-sin(@0)*cos(@1)",RooArgList(*phi_s,*deltaPe));
    RooFormulaVar *coef63=new RooFormulaVar("coef63","coef63","-sin(@0)*@1*(1-2*@2)",RooArgList(*deltaPe,*tag,*mistag));
    RooFormulaVar *coef64=new RooFormulaVar("coef64","coef64","cos(@0)*cos(@1)*@2*(1-2*@3)",RooArgList(*phi_s,*deltaPe,*tag,*mistag));
    RooFormulaVar *coef72=new RooFormulaVar("coef72","coef72","cos(@0)",RooArgList(*phi_s));
    RooFormulaVar *coef74=new RooFormulaVar("coef74","coef74","sin(@0)*@1*(1-2*@2)",RooArgList(*phi_s,*tag,*mistag));
    RooFormulaVar *coef82=new RooFormulaVar("coef82","coef82","-sin(@0)*sin(@1-@2-@3)",RooArgList(*phi_s,*deltaPa,*deltaSPe,*deltaPe));
    RooFormulaVar *coef83=new RooFormulaVar("coef83","coef83","-cos(@0-@1-@2)*@3*(1-2*@4)",RooArgList(*deltaPa,*deltaSPe,*deltaPe,*tag,*mistag));
    RooFormulaVar *coef84=new RooFormulaVar("coef84","coef84","sin(@0-@1-@2)*@3*(1-2*@4)*cos(@5)",RooArgList(*deltaPa,*deltaSPe,*deltaPe,*tag,*mistag,*phi_s));
    RooFormulaVar *coef91=new RooFormulaVar("coef91","coef91","sin(-@0)",RooArgList(*deltaSPe));
    RooFormulaVar *coef92=new RooFormulaVar("coef92","coef92","sin(-@0)*cos(@1)",RooArgList(*deltaSPe,*phi_s));
    RooFormulaVar *coef94=new RooFormulaVar("coef94","coef94","sin(-@0)*sin(@1)*@2*(1-2*@3)",RooArgList(*deltaSPe,*phi_s,*tag,*mistag));
    RooFormulaVar *coef102=new RooFormulaVar("coef102","coef102","-sin(@0)*sin(-@1-@2)",RooArgList(*phi_s,*deltaSPe,*deltaPe));
    RooFormulaVar *coef103=new RooFormulaVar("coef103","coef103","-cos(-@0-@1)*@2*(1-2*@3)",RooArgList(*deltaSPe,*deltaPe,*tag,*mistag));
    RooFormulaVar *coef104=new RooFormulaVar("coef104","coef104","sin(-@0-@1)*@2*(1-2*@3)*cos(@4)",RooArgList(*deltaSPe,*deltaPe,*tag,*mistag,*phi_s));

    RooProduct *eq11=new RooProduct("eq11","amp0*g1",RooArgSet(*coshGConv,*fun1T));
    RooProduct *eq12=new RooProduct("eq12","amp0*g1",RooArgSet(*coef12,*sinhGConv,*fun1T));
    RooProduct *eq14=new RooProduct("eq14","amp0*g1",RooArgSet(*coef14,*sinGConv,*fun1T));
    RooProduct *eq21=new RooProduct("eq21","amp0*g1",RooArgSet(*coshGConv,*fun2T));
    RooProduct *eq22=new RooProduct("eq22","amp0*g1",RooArgSet(*coef22,*sinhGConv,*fun2T));
    RooProduct *eq24=new RooProduct("eq24","amp0*g1",RooArgSet(*coef24,*sinGConv,*fun2T));
    RooProduct *eq31=new RooProduct("eq31","amp0*g1",RooArgSet(*coshGConv,*fun3T));
    RooProduct *eq32=new RooProduct("eq32","amp0*g1",RooArgSet(*coef32,*sinhGConv,*fun3T));
    RooProduct *eq34=new RooProduct("eq34","amp0*g1",RooArgSet(*coef34,*sinGConv,*fun3T));
    RooProduct *eq42=new RooProduct("eq42","amp0*g1",RooArgSet(*coef42,*sinhGConv,*fun4T));
    RooProduct *eq43=new RooProduct("eq43","amp0*g1",RooArgSet(*coef43,*cosGConv,*fun4T));
    RooProduct *eq44=new RooProduct("eq44","amp0*g1",RooArgSet(*coef44,*sinGConv,*fun4T));
    RooProduct *eq51=new RooProduct("eq51","amp0*g1",RooArgSet(*coef51,*coshGConv,*fun5T));
    RooProduct *eq52=new RooProduct("eq52","amp0*g1",RooArgSet(*coef52,*sinhGConv,*fun5T));
    RooProduct *eq54=new RooProduct("eq54","amp0*g1",RooArgSet(*coef54,*sinGConv,*fun5T));
    RooProduct *eq62=new RooProduct("eq62","amp0*g1",RooArgSet(*coef62,*sinhGConv,*fun6T));
    RooProduct *eq63=new RooProduct("eq63","amp0*g1",RooArgSet(*coef63,*cosGConv,*fun6T));
    RooProduct *eq64=new RooProduct("eq64","amp0*g1",RooArgSet(*coef64,*sinGConv,*fun6T));
    RooProduct *eq71=new RooProduct("eq71","amp0*g1",RooArgSet(*coshGConv,*fun7T));
    RooProduct *eq72=new RooProduct("eq72","amp0*g1",RooArgSet(*coef72,*sinhGConv,*fun7T));
    RooProduct *eq74=new RooProduct("eq74","amp0*g1",RooArgSet(*coef74,*sinGConv,*fun7T));
    RooProduct *eq82=new RooProduct("eq82","amp0*g1",RooArgSet(*coef82,*sinhGConv,*fun8T));
    RooProduct *eq83=new RooProduct("eq83","amp0*g1",RooArgSet(*coef83,*cosGConv,*fun8T));
    RooProduct *eq84=new RooProduct("eq84","amp0*g1",RooArgSet(*coef84,*sinGConv,*fun8T));
    RooProduct *eq91=new RooProduct("eq91","amp0*g1",RooArgSet(*coef91,*coshGConv,*fun9T));
    RooProduct *eq92=new RooProduct("eq92","amp0*g1",RooArgSet(*coef92,*sinhGConv,*fun9T));
    RooProduct *eq94=new RooProduct("eq94","amp0*g1",RooArgSet(*coef94,*sinGConv,*fun9T));
    RooProduct *eq102=new RooProduct("eq102","amp0*g1",RooArgSet(*coef102,*sinhGConv,*fun10T));
    RooProduct *eq103=new RooProduct("eq103","amp0*g1",RooArgSet(*coef103,*cosGConv,*fun10T));
    RooProduct *eq104=new RooProduct("eq104","amp0*g1",RooArgSet(*coef104,*sinGConv,*fun10T));

    RooArgList funtot(*eq11,*eq12,*eq14,*eq21,*eq22,*eq24,*eq31,*eq32,*eq34);
    funtot.add(*eq42);
    funtot.add(*eq43);
    funtot.add(*eq44);
    funtot.add(*eq51);
    funtot.add(*eq52);
    funtot.add(*eq54);
    funtot.add(*eq62);
    funtot.add(*eq63);
    funtot.add(*eq64);
    RooArgList coeftot(*A_0,*A_0,*A_0,*A_pa,*A_pa,*A_pa,*A_pe,*A_pe,*A_pe);
    coeftot.add(*A_4);
    coeftot.add(*A_4);
    coeftot.add(*A_4);
    coeftot.add(*A_5);
    coeftot.add(*A_5);
    coeftot.add(*A_5);
    coeftot.add(*A_6);
    coeftot.add(*A_6);
    coeftot.add(*A_6);

    RooArgList funtotP(*eq11,*eq12,*eq14,*eq21,*eq22,*eq24,*eq31,*eq32,*eq34);
    funtotP.add(*eq42);
    funtotP.add(*eq43);
    funtotP.add(*eq44);
    funtotP.add(*eq51);
    funtotP.add(*eq52);
    funtotP.add(*eq54);
    funtotP.add(*eq62);
    funtotP.add(*eq63);
    funtotP.add(*eq64);
    funtotP.add(*eq71);
    funtotP.add(*eq72);
    funtotP.add(*eq74);
    funtotP.add(*eq82);
    funtotP.add(*eq83);
    funtotP.add(*eq84);
    funtotP.add(*eq91);
    funtotP.add(*eq92);
    funtotP.add(*eq94);
    funtotP.add(*eq102);
    funtotP.add(*eq103);
    funtotP.add(*eq104);
    RooArgList coeftotP(*A_0,*A_0,*A_0,*A_pa,*A_pa,*A_pa,*A_pe,*A_pe,*A_pe);
    coeftotP.add(*A_4);
    coeftotP.add(*A_4);
    coeftotP.add(*A_4);
    coeftotP.add(*A_5);
    coeftotP.add(*A_5);
    coeftotP.add(*A_5);
    coeftotP.add(*A_6);
    coeftotP.add(*A_6);
    coeftotP.add(*A_6);
    coeftotP.add(*A_S);
    coeftotP.add(*A_S);
    coeftotP.add(*A_S);
    coeftotP.add(*A_8);
    coeftotP.add(*A_8);
    coeftotP.add(*A_8);
    coeftotP.add(*A_9);
    coeftotP.add(*A_9);
    coeftotP.add(*A_9);
    coeftotP.add(*A_10);
    coeftotP.add(*A_10);
    coeftotP.add(*A_10);

    RooRealSumPdf *PDFdefP=new RooRealSumPdf("PDFdefP","Signal PDF Amp_i*funiT",funtotP,coeftotP);
    RooRealSumPdf *PDFdef=new RooRealSumPdf("PDFdef","Signal PDF Amp_i*funiT",funtot,coeftot);
    //Constructions of the Signal PDF*efficiencies=Time_dependent_functions* Angles_dependent_functions= Amp_i*funiT*eff(t)

   
    RooFormulaVar tagfun("tagfun","tagfum","(1-@0*(1-2*@1))/2",RooArgList(*tag,*mistag));
    RooFormulaVar tagfunbar("tagfunbar","tagfumbar","(1+@0*(1-2*@1))/2",RooArgList(*tag,*mistag));
    RooProduct *pdf1bar=new RooProduct("pdf1bar","amp0*g1",RooArgSet(*myAmp0bar,*fun1T,tagfunbar));// eff_ct + *poly3D
    RooProduct *pdf2bar=new RooProduct("pdf2bar","amp0*g2",RooArgSet(*myAmp0bar,*fun2T,tagfunbar));//
    RooProduct *pdf3bar=new RooProduct("pdf3bar","ampe*g3",RooArgSet(*myAmpebar,*fun3T,tagfunbar));//
    RooProduct *pdf4bar=new RooProduct("pdf4bar","Ampa4*fun4T",RooArgSet(*myAmpa4bar,*fun4T,tagfunbar));   
    RooProduct *pdf5bar=new RooProduct("pdf5bar","ampa0*g5",RooArgSet(*myAmp0bar,*fun5T,tagfunbar));//
    RooProduct *pdf6bar=new RooProduct("pdf6bar","Ampa6*fun6T",RooArgSet(*myAmpa6bar,*fun6T,tagfunbar)); 
    RooProduct *pdf1=new RooProduct("pdf1","amp0*g1",RooArgSet(*myAmp0,*fun1T,tagfun));// eff_ct + *poly3D
    RooProduct *pdf2=new RooProduct("pdf2","amp0*g2",RooArgSet(*myAmp0,*fun2T,tagfun));//
    RooProduct *pdf3=new RooProduct("pdf3","ampe*g3",RooArgSet(*myAmpe,*fun3T,tagfun));//
    RooProduct *pdf4=new RooProduct("pdf4","Ampa4*fun4T",RooArgSet(*myAmpa4,*fun4T,tagfun));   
    RooProduct *pdf5=new RooProduct("pdf5","ampa0*g5",RooArgSet(*myAmp0,*fun5T,tagfun));//
    RooProduct *pdf6=new RooProduct("pdf6","Ampa6*fun6T",RooArgSet(*myAmpa6,*fun6T,tagfun));
    RooProduct *pdf1tag=new RooProduct("pdf1tag","amp0*g1",RooArgSet(*myAmp0tag,*fun1T));//,*eff_ct));// eff_ct + *poly3D
    RooProduct *pdf2tag=new RooProduct("pdf2tag","amp0*g2",RooArgSet(*myAmp0tag,*fun2T));//,*eff_ct));//
    RooProduct *pdf3tag=new RooProduct("pdf3tag","ampe*g3",RooArgSet(*myAmpetag,*fun3T));//,*eff_ct));//
    RooProduct *pdf4tag=new RooProduct("pdf4tag","Ampa4*fun4T",RooArgSet(*myAmpa4tag,*fun4T,*sindpadpe));//,*eff_ct));   
    RooProduct *pdf5tag=new RooProduct("pdf5tag","ampa0*g5",RooArgSet(*myAmpa5tag,*fun5T,*cosdpa));//,*eff_ct));//
    RooProduct *pdf6tag=new RooProduct("pdf6tag","Ampa6*fun6T",RooArgSet(*myAmpa6tag,*fun6T,*sindpe));//,*eff_ct)); 
    RooProduct *pdf1untag=new RooProduct("pdf1untag","pdf1untag",RooArgSet(*myAmp0untag,*fun1T));//,*eff_ct));//,tagfun));// eff_ct + *poly3D
    RooProduct *pdf2untag=new RooProduct("pdf2untag","pdf2untag",RooArgSet(*myAmp0untag,*fun2T));//,*eff_ct));//,tagfun));//
    RooProduct *pdf3untag=new RooProduct("pdf3untag","pdf3untag",RooArgSet(*myAmpeuntag,*fun3T));//,*eff_ct));//,tagfun));//
    RooProduct *pdf4untag=new RooProduct("pdf4untag","pdf4untag",RooArgSet(*myAmpa4untag,*fun4T));//,*eff_ct));//,tagfun));   
    RooProduct *pdf5untag=new RooProduct("pdf5untag","pdf5untag",RooArgSet(*myAmpa5untag,*fun5T));//,*eff_ct));//,tagfun));//
    RooProduct *pdf6untag=new RooProduct("pdf6untag","pdf6untag",RooArgSet(*myAmpa6untag,*fun6T));//,*eff_ct));//,tagfun)); 

/*
    RooProduct *pdf1=new RooProduct("pdf1","amp0*g1",RooArgSet(*myAmp0,*fun1T));// eff_ct + *poly3D
    RooProduct *pdf2=new RooProduct("pdf2","amp0*g2",RooArgSet(*myAmp0,*fun2T));//
    RooProduct *pdf3=new RooProduct("pdf3","ampe*g3",RooArgSet(*myAmpe,*fun3T));//
    RooProduct *pdf4=new RooProduct("pdf4","Ampa4*fun4T",RooArgSet(*myAmpa4,*fun4T));
    RooProduct *pdf5=new RooProduct("pdf5","ampa0*g5",RooArgSet(*myAmp0,*fun5T,cosdpa));//
    RooProduct *pdf6=new RooProduct("pdf6","Ampa6*fun6T",RooArgSet(*myAmpa6,*fun6T));
*/
    //RooRealSumPdf *PDF=new RooRealSumPdf("PDF","Signal PDF Amp_i*funiT",RooArgList(*pdf1,*pdf2,*pdf3,*pdf4,*pdf5,*pdf6),RooArgList(*A_0,*A_pa,*A_pe,*A_4,*A_5,*A_6));//,Conditional(gm,*BdCt2DBS)); //FULL MODEL
//    RooRealSumPdf *PDF=new RooRealSumPdf("PDF","Signal PDF Amp_i*funiT",RooArgList(*pdf1,*pdf2,*pdf3,*pdf4,*pdf5,*pdf6),RooArgList(*A_0,*A_3,*A_2,*A_4,*A_5,*A_6));//,Conditional(gm,*BdCt2DBS));
    RooRealSumPdf *PDFBsbar=new RooRealSumPdf("PDFBsbar","Signal PDF Amp_i*funiT",RooArgList(*pdf1bar,*pdf2bar,*pdf3bar,*pdf4bar,*pdf5bar,*pdf6bar),RooArgList(*A_0,*A_pa,*A_pe,*A_4,*A_5,*A_6));//,Conditional(gm,*BdCt2DBS));  
    RooRealSumPdf *PDFBs=new RooRealSumPdf("PDFBs","Signal PDF Amp_i*funiT",RooArgList(*pdf1,*pdf2,*pdf3,*pdf4,*pdf5,*pdf6),RooArgList(*A_0,*A_pa,*A_pe,*A_4,*A_5,*A_6));//,Conditional(gm,*BdCt2DBS));  
    RooRealSumPdf *PDFtag=new RooRealSumPdf("PDFtag","Signal PDF Amp_i*funiT",RooArgList(*pdf1tag,*pdf2tag,*pdf3tag,*pdf4tag,*pdf5tag,*pdf6tag),RooArgList(*A_0,*A_pa,*A_pe,*A_4,*A_5,*A_6));//,Conditional(gm,*BdCt2DBS));  
    RooRealSumPdf *PDFuntag=new RooRealSumPdf("PDFuntag","Signal PDF Amp_i*funiT",RooArgList(*pdf1untag,*pdf2untag,*pdf3untag,*pdf4untag,*pdf5untag,*pdf6untag),RooArgList(*A_0,*A_pa,*A_pe,*A_4,*A_5,*A_6));//,Conditional(gm,*BdCt2DBS));  
    RooRealVar *nBs=new RooRealVar("nBs","nBs",0,1000000);
    RooRealVar *nBsbar=new RooRealVar("nBsbar","nBsbar",0,1000000);
    RooAddPdf *PDF=new RooAddPdf("PDF","Signal model",RooArgSet(*PDFBsbar,*PDFBs),RooArgSet(*nBsbar,*nBs));
//    RooRealSumPdf *PDF=new RooRealSumPdf("PDF","Signal model",RooArgSet(*PDFBsbar,*PDFBs),RooArgList(RooConst(1),RooConst(1)));


    RooWorkspace *phismodelMC = new RooWorkspace("phismodelMC","phismodelMC") ;
    phismodelMC->import(*PDFtag);
    phismodelMC->import(*PDFuntag,RecycleConflictNodes());
    phismodelMC->import(*PDFdef,RecycleConflictNodes());
    phismodelMC->import(*PDFdefP,RecycleConflictNodes());
    phismodelMC->autoImportClassCode("true");
 //   phismodelMC->import(Alldataset);
    phismodelMC->Print("v");
    cout<<"why?"<<endl;
    phismodelMC->writeToFile((target + "phismodelMC.root").c_str());

}
