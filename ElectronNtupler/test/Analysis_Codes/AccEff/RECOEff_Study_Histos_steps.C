#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <math.h>
#include <vector>
#include <iostream>

Double_t deltaPhi(Double_t phi1, Double_t phi2)
{
  Double_t pi = 3.1415927;
  Double_t dphi = fabs(phi1 - phi2);
  if(dphi >= pi) dphi = 2. * pi - dphi;
  return dphi;
}

Double_t deltaEta(Double_t eta1, Double_t eta2)
{
  Double_t deta = fabs(eta1 - eta2);
  return deta;
}

Double_t deltaR(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2)
{
  Double_t dEta = deltaEta(eta1, eta2);
  Double_t dPhi = deltaPhi(phi1, phi2);
  Double_t dr = sqrt(dEta*dEta + dPhi*dPhi);
  return dr;
}

Bool_t isPassAccCondition_GenLepton_ECALGAP(double leadpt, double leadeta, double subleadpt, double subleadeta)
{
	Bool_t isPassAcc = kFALSE;


	if( leadpt > 30. && fabs(leadeta) < 2.5 && !( fabs(leadeta) > 1.4442 && fabs(leadeta) < 1.566 ) &&
		subleadpt  > 10.  && fabs(subleadeta) < 2.5 && !( fabs(subleadeta) > 1.4442 && fabs(subleadeta) < 1.566 ) )
		isPassAcc = 1;

	return isPassAcc;
}

Bool_t isPassAccCondition_Electron(double leadpt, double leadeta, double subleadpt, double subleadeta)
{

	Bool_t isPassAcc = kFALSE;
	if( leadpt > 30 && fabs(leadeta) < 2.5 && !( fabs(leadeta) > 1.4442 && fabs(leadeta) < 1.566 ) &&
		subleadpt  > 10  && fabs(subleadeta)  < 2.5 && !( fabs(subleadeta) > 1.4442 && fabs(subleadeta) < 1.566 ) )
		isPassAcc = kTRUE;

	return isPassAcc;
}

void RECOEff_Study_Histos_steps(int debug=0) {


     TString workdir;
     std::vector<TFile*> InputFiles_signal_DY;

double xsec[11] = {18609.9/3,5789./3,226./3,7.68/3,0.423/3,0.24/3,0.036/3,0.03/3,0.0159/3,0.00201/3,0.00054/3};
double sumweights[11] = {758771667841.549683,150453917158.292908,219884886.214768,7008726.024242,122990.185607,69998.861766,11080.956088,9415.627747,4893.463037,664.033593,159.615701};
int mass[11] = {10,50,100,200,400,500,700,800,1000,1500,2000};

 workdir = "/tmp/arun/eos/cms/store/group/phys_higgs/cmshww/arun/DYAnalysis_76X_Calibrated/DY_Signal/";
 //workdir = "/home/andriusj/DY13TeV/DYAnalysis_76X_Calibrated/";
 workdir = "/media/sf_Share2/DYee_76X_Calibrated/";
 InputFiles_signal_DY.clear();

 InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_10to50.root"));
 InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_50toInf.root"));
 InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_100to200.root"));
 InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_200to400.root"));
 InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_400to500.root"));
 InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_500to700.root"));
 InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_700to800.root"));
 InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_800to1000.root"));
 InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_1000to1500.root"));
 InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_1500to2000.root"));
 InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_2000to3000.root"));

 const Int_t nMassBin = 43;

 Double_t MassBinEdges[44] = {15,20,25,30,35,40,45,50,55,60,64,68,72,76,81,86,91,96,101,106,110,115,120,126,133,141,150,160,171,185,200,220,243,273,320,380,440,510,600,700,830,1000,1500,3000};

 TH1D* h1preFSR_1GeV= new TH1D("h1preFSR_1GeV","preFSR (no cuts);M_{ee} [GeV];weighted count",3100,0,3100);
 TH1D* h1preFSR= new TH1D("h1preFSR","preFSR;M_{ee} [GeV];weighted count",nMassBin,MassBinEdges);

 h1preFSR_1GeV->Sumw2();
 h1preFSR_1GeV->SetDirectory(0);
 h1preFSR->Sumw2();
 h1preFSR->SetDirectory(0);

 TH1D* h1_eff_sumPass= new TH1D("h1_eff_sumPass","selection passing events;M_{ee} [GeV];weighted count",nMassBin,MassBinEdges);
 TH1D* h1_eff_sumTot= new TH1D("h1_eff_sumTot","total events in acceptance;M_{ee} [GeV];weighted count",nMassBin,MassBinEdges);
 TH1D* h1_acc_sumPass= new TH1D("h1_acc_sumPass","sum events passing acceptance cuts;M_{ee} [GeV];weighted count",nMassBin,MassBinEdges);
 TH1D* h1_acc_sumTot= new TH1D("h1_acc_sumTot","total sum postFSR events;M_{ee} [GeV];weighted count",nMassBin,MassBinEdges);
 h1_eff_sumPass->Sumw2(); h1_eff_sumPass->SetDirectory(0);
 h1_eff_sumTot->Sumw2(); h1_eff_sumTot->SetDirectory(0);
 h1_acc_sumPass->Sumw2(); h1_acc_sumPass->SetDirectory(0);
 h1_acc_sumTot->Sumw2(); h1_acc_sumTot->SetDirectory(0);

 TH1D* h1_eff_sumPass_chk= new TH1D("h1_eff_sumPass_chk","selection passing events chk;M_{ee} [GeV];weighted count",nMassBin,MassBinEdges);
 TH1D* h1_eff_sumPass_noPU= new TH1D("h1_eff_sumPass_noPU","selection passing events (noPU);M_{ee} [GeV];weighted (noPU) count",nMassBin,MassBinEdges);
 TH1D *h1_eff_sumPass_noTrigObjMatching= new TH1D("h1_eff_sumPass_noTrigObjMatching","selection passing events (noTrigObjMatching);M_{ee} [GeV];weighted (noPU) count",nMassBin,MassBinEdges);
 TH1D* h1_eff_sumPass_matchDR3= new TH1D("h1_eff_sumPass_matchDR3","selection passing events (matchDR3);M_{ee} [GeV];weighted (matchDR3) count",nMassBin,MassBinEdges);
 TH1D* h1_eff_sumPass_noTrig= new TH1D("h1_eff_sumPass_noTrig","selection passing events (noTrig);M_{ee} [GeV];weighted count",nMassBin,MassBinEdges);
 TH1D *h1_eff_sumPass_noTrigNoPU= new TH1D("h1_eff_sumPass_noTrigNoPU","selection passing events (noTrigNoPU);M_{ee} [GeV];weighted (noPU) count",nMassBin,MassBinEdges);
 h1_eff_sumPass_chk->Sumw2(); h1_eff_sumPass_chk->SetDirectory(0);
 h1_eff_sumPass_noPU->Sumw2(); h1_eff_sumPass_noPU->SetDirectory(0);
 h1_eff_sumPass_noTrigObjMatching->Sumw2(); h1_eff_sumPass_noTrigObjMatching->SetDirectory(0);
 h1_eff_sumPass_matchDR3->Sumw2(); h1_eff_sumPass_matchDR3->SetDirectory(0);
 h1_eff_sumPass_noTrig->Sumw2(); h1_eff_sumPass_noTrig->SetDirectory(0);
 h1_eff_sumPass_noTrigNoPU->Sumw2(); h1_eff_sumPass_noTrigNoPU->SetDirectory(0);

 //unsigned int nsample = InputFiles_signal_DY.size();
  TFile* file[11];

  for(unsigned int j = 0; j < InputFiles_signal_DY.size(); ++j ) {

    TString codePath="/afs/cern.ch/work/a/arun/DYAnalysis2015_76X_v1/CMSSW_7_6_3/src/EgammaWork/";
    codePath="/home/andriusj/DY13TeV/DYanalysis-20160817/";
    TFile *f1 = TFile::Open(codePath + "ElectronNtupler/test/Analysis_Codes/Data-MC/dataPUDist.root");
    TFile *f2 = TFile::Open(codePath + "ElectronNtupler/test/Analysis_Codes/Data-MC/PileUp_MC.root");

//data PU histo
    TH1F *DATA_puDist = (TH1F*)f1->Get("pileup");
    DATA_puDist->Scale(1/DATA_puDist->Integral());

//mc PU histo
    TH1F *MC_puDist = (TH1F*)f2->Get("pileup_MC");
    TH1F *weights = (TH1F*)DATA_puDist->Clone("weights");
    weights->Divide(MC_puDist);


  TTree * tmpTree = (TTree*)InputFiles_signal_DY.at(j)->Get("ntupler/ElectronTree");

  vector<float>   *genPreFSR_Pt;
  vector<float>   *genPreFSR_Eta;
  vector<float>   *genPreFSR_Rap;
  vector<float>   *genPreFSR_Phi;
  vector<float>   *genPreFSR_En;
  vector<float>   *genPostFSR_Pt;
  vector<float>   *genPostFSR_Eta;
  vector<float>   *genPostFSR_Rap;
  vector<float>   *genPostFSR_Phi;
  vector<float>   *genPostFSR_En;
  vector<int>     *passMediumId;
  vector<int>     *passTightId;
  vector<int>     *passHEEPId;
  vector<float>   *ptElec;
  vector<float>   *etaElec;
  vector<float>   *rapElec;
  vector<float>   *phiElec;
  vector<float>   *energyElec;
  vector<float>   *etaSC;
  Int_t           tauFlag; 
  Double_t        theWeight;
  Int_t           nPUTrue;
  Bool_t          Ele23_WPLoose;
  vector<double>  *pt_Ele23;
  vector<double>  *eta_Ele23;
  vector<double>  *phi_Ele23;

  genPreFSR_Pt = 0;
  genPreFSR_Eta = 0;
  genPreFSR_Rap = 0;
  genPreFSR_Phi = 0;
  genPreFSR_En = 0;
  genPostFSR_Pt = 0;
  genPostFSR_Eta = 0;
  genPostFSR_Rap = 0;
  genPostFSR_Phi = 0;
  genPostFSR_En = 0;
  ptElec = 0;
  etaElec = 0;
  rapElec = 0;
  phiElec = 0;
  energyElec = 0;
  etaSC = 0;
  passMediumId = 0;
  passTightId = 0;
  passHEEPId = 0;
  pt_Ele23 = 0;
  eta_Ele23 = 0;
  phi_Ele23 = 0;

    tmpTree->SetBranchStatus("*", 0);
    tmpTree->SetBranchStatus("genPreFSR_Pt", 1);
    tmpTree->SetBranchStatus("genPreFSR_Eta", 1);
    tmpTree->SetBranchStatus("genPreFSR_Rap", 1);
    tmpTree->SetBranchStatus("genPreFSR_Phi", 1);
    tmpTree->SetBranchStatus("genPreFSR_En", 1);
    tmpTree->SetBranchStatus("genPostFSR_Pt", 1);
    tmpTree->SetBranchStatus("genPostFSR_Eta", 1);
    tmpTree->SetBranchStatus("genPostFSR_Rap", 1);
    tmpTree->SetBranchStatus("genPostFSR_Phi", 1);
    tmpTree->SetBranchStatus("genPostFSR_En", 1);
    tmpTree->SetBranchStatus("ptElec", 1);
    tmpTree->SetBranchStatus("etaElec", 1);
    tmpTree->SetBranchStatus("phiElec", 1);
    tmpTree->SetBranchStatus("energyElec", 1);
    tmpTree->SetBranchStatus("etaSC", 1);
    tmpTree->SetBranchStatus("passMediumId", 1);
    tmpTree->SetBranchStatus("tauFlag", 1);
    tmpTree->SetBranchStatus("theWeight", 1);
    tmpTree->SetBranchStatus("nPUTrue", 1);
    tmpTree->SetBranchStatus("pt_Ele23", 1);
    tmpTree->SetBranchStatus("eta_Ele23", 1);
    tmpTree->SetBranchStatus("phi_Ele23", 1);
    tmpTree->SetBranchStatus("Ele23_WPLoose", 1);

  tmpTree->SetBranchAddress("genPreFSR_Pt", &genPreFSR_Pt);
  tmpTree->SetBranchAddress("genPreFSR_Eta", &genPreFSR_Eta);
  tmpTree->SetBranchAddress("genPreFSR_Rap", &genPreFSR_Rap);
  tmpTree->SetBranchAddress("genPreFSR_Phi", &genPreFSR_Phi);
  tmpTree->SetBranchAddress("genPreFSR_En", &genPreFSR_En);
  tmpTree->SetBranchAddress("genPostFSR_Pt", &genPostFSR_Pt);
  tmpTree->SetBranchAddress("genPostFSR_Eta", &genPostFSR_Eta);
  tmpTree->SetBranchAddress("genPostFSR_Rap", &genPostFSR_Rap);
  tmpTree->SetBranchAddress("genPostFSR_Phi", &genPostFSR_Phi);
  tmpTree->SetBranchAddress("genPostFSR_En", &genPostFSR_En);
  tmpTree->SetBranchAddress("ptElec", &ptElec);
  tmpTree->SetBranchAddress("etaElec", &etaElec);
  tmpTree->SetBranchAddress("rapElec", &rapElec);
  tmpTree->SetBranchAddress("phiElec", &phiElec);
  tmpTree->SetBranchAddress("energyElec", &energyElec);
  tmpTree->SetBranchAddress("etaSC", &etaSC);
  tmpTree->SetBranchAddress("passMediumId", &passMediumId);
  tmpTree->SetBranchAddress("passTightId", &passTightId);
  tmpTree->SetBranchAddress("passHEEPId", &passHEEPId);
  tmpTree->SetBranchAddress("tauFlag", &tauFlag);
  tmpTree->SetBranchAddress("theWeight", &theWeight);
  tmpTree->SetBranchAddress("nPUTrue", &nPUTrue);
  tmpTree->SetBranchAddress("Ele23_WPLoose", &Ele23_WPLoose);
  tmpTree->SetBranchAddress("pt_Ele23", &pt_Ele23);
  tmpTree->SetBranchAddress("eta_Ele23", &eta_Ele23);
  tmpTree->SetBranchAddress("phi_Ele23", &phi_Ele23);

file[j] = new TFile(Form("DY_%d_EE_76X_forAccEff_v1steps.root",mass[j]),"RECREATE");

 	TH1D *h_mass_AccTotal = new TH1D("h_mass_AccTotal", "", nMassBin, MassBinEdges);
	TH1D *h_mass_AccPass = new TH1D("h_mass_AccPass", "", nMassBin, MassBinEdges);
	TH1D *h_mass_EffTotal = new TH1D("h_mass_EffTotal", "", nMassBin, MassBinEdges);
	TH1D *h_mass_EffPass = new TH1D("h_mass_EffPass", "", nMassBin, MassBinEdges);

  double dR, massGen;
  TLorentzVector gen1,gen2,diGen;
  TLorentzVector numEle1, numEle2, numDiele, numEle1HEEP, numEle2HEEP, numDieleHEEP, denEle1, denEle2, denDiele;
  TLorentzVector genEle1, genEle2, genDiele;
  TLorentzVector ele1,ele2,dielectron;


// for id corr
int mediumId,passId;
int heepId;

//for sorted electrons
  vector <double> newelePt; 
  vector <double> neweleEta; 
  vector <double> neweleEnr; 
  vector <double> newelePhi;
  vector <double> neweleMediumId;
  vector <double> newSCeta;

//for gen electrons, postFSR
  vector <double> newgenPt; 
  vector <double> newgenEta; 
  vector <double> newgenEnr; 
  vector <double> newgenPhi;

//for gen matched sorted electrons
  vector <double> elePt; 
  vector <double> eleEta; 
  vector <double> eleEnr; 
  vector <double> elePhi;
  vector <double> eleMediumId;

double lumiWeight = (xsec[j]/sumweights[j])*2316.969;
//double lumiWeight = 6015/451498239235.;
//double lumiWeight = 0.036/11571.8555;//700-800
//double lumiWeight = 0.03/9268.827694;// 800-1000
//double lumiWeight = 0.0159/5101.239783;// 1000-1500
//double lumiWeight = 0.00201/733.451968; //1500-2000
//double lumiWeight = 0.00054/177.265536; //2000-3000

cout << "Lumiweight = " << lumiWeight << endl;

 unsigned int nentries = tmpTree->GetEntries();
 cout<<"entries: "<<nentries<<endl;

 for (unsigned int k=0; k < nentries; k++) {
    tmpTree->GetEntry(k);

    if (debug && (k>=10000)) {
      std::cout << "debug mode: stopping loop at k=" << k << "\n";
      break;
    }
    if(k%1000000 == 0){
      cout << "Events Processed :  " << k
	   << Form(" (%4.1lf%c)",k*100./double(nentries),'%') << endl;
    }

//PUWeight
      int bin = 0;
      double puWeights = 1.0;
      bin = weights->GetXaxis()->FindBin(nPUTrue);
      puWeights = weights->GetBinContent(bin);


// sorting of reco electrons
    int index1[ptElec->size()];
    float pt1[ptElec->size()];

    for(unsigned int el=0; el<ptElec->size(); el++)
    {
      pt1[el]=ptElec->at(el);
    }

    int sizer = sizeof(pt1)/sizeof(pt1[0]);
    TMath::Sort(sizer,pt1,index1,true);

// sorting of gen electrons : postFSR
    int index2[genPostFSR_Pt->size()];
    float pt2[genPostFSR_Pt->size()];

    for(unsigned int b=0; b<genPostFSR_Pt->size(); b++)
    {
      pt2[b]=genPostFSR_Pt->at(b);
    }
    int sizen = sizeof(pt2)/sizeof(pt2[0]);
    TMath::Sort(sizen,pt2,index2,true);

// sorting of gen electrons : preFSR

    int index3[genPreFSR_Pt->size()];
    float pt3[genPreFSR_Pt->size()];
    for(unsigned int c=0; c<genPreFSR_Pt->size(); c++)
    { pt3[c]=genPreFSR_Pt->at(c);}
    int size1 = sizeof(pt3)/sizeof(pt3[0]);
    TMath::Sort(size1,pt3,index3,true);

//clearing of vectors
    mediumId = 0,  passId = 0, heepId = 0;
    dR = 0.; massGen = 0.0;
    newelePt.clear(); 
    neweleEta.clear(); 
    neweleEnr.clear(); 
    newelePhi.clear(); 
    neweleMediumId.clear(); 
    newSCeta.clear();

    newgenPt.clear(); 
    newgenEta.clear(); 
    newgenEnr.clear(); 
    newgenPhi.clear();

    elePt.clear(); 
    eleEta.clear(); 
    eleEnr.clear(); 
    elePhi.clear(); 
    eleMediumId.clear();

    //unsigned int matchGen = 0;
    //unsigned int matchReco = 0;
    double postFSRmass = 0.;
    double recomass = 0.;

    if(genPreFSR_Pt->size() == 2){
      gen1.SetPtEtaPhiE(genPreFSR_Pt->at(0),genPreFSR_Eta->at(0),genPreFSR_Phi->at(0),genPreFSR_En->at(0));
      gen2.SetPtEtaPhiE(genPreFSR_Pt->at(1),genPreFSR_Eta->at(1),genPreFSR_Phi->at(1),genPreFSR_En->at(1));
      diGen=gen1+gen2;
      massGen=diGen.M();
    }

if(j==1 && massGen >= 100) continue;
double dR1(0.),dR2(0.);
    if(!tauFlag && genPreFSR_Pt->size() == 2){

      h1preFSR_1GeV->Fill( massGen, theWeight*lumiWeight );
      h1preFSR->Fill( massGen, theWeight*lumiWeight );


      if(genPostFSR_Pt->size() >= 2.){
	for(unsigned int j=0;j<genPostFSR_Eta->size();j++){
	    newgenPt.push_back(genPostFSR_Pt->at(index2[j]));
	    newgenEta.push_back(genPostFSR_Eta->at(index2[j]));
	    newgenEnr.push_back(genPostFSR_En->at(index2[j]));
	    newgenPhi.push_back(genPostFSR_Phi->at(index2[j]));
	  }
 	}

       Bool_t Flag_PassAcc = kFALSE;
       genEle1.SetPtEtaPhiE(newgenPt.at(0),newgenEta.at(0),newgenPhi.at(0),newgenEnr.at(0));
       genEle2.SetPtEtaPhiE(newgenPt.at(1),newgenEta.at(1),newgenPhi.at(1),newgenEnr.at(1));

          genDiele = genEle1+genEle2;
          postFSRmass = genDiele.M();

          Flag_PassAcc = isPassAccCondition_GenLepton_ECALGAP(newgenPt.at(0),newgenEta.at(0),newgenPt.at(1),newgenEta.at(1));
      if(Flag_PassAcc) {
         h_mass_AccTotal->Fill(postFSRmass,theWeight*lumiWeight);
         h_mass_AccPass->Fill(postFSRmass,theWeight*lumiWeight);
	 h1_acc_sumPass->Fill(postFSRmass,theWeight*lumiWeight);
	 h1_acc_sumTot->Fill(postFSRmass, theWeight*lumiWeight);
         }
      else
      {
           h_mass_AccTotal->Fill(postFSRmass,theWeight*lumiWeight);
	   h1_acc_sumTot->Fill(postFSRmass, theWeight*lumiWeight);
      }

if(Flag_PassAcc == kTRUE) {
Bool_t Flag_PassEff = kFALSE;

      if(ptElec->size()>=2.){
        for(unsigned int i=0;i<ptElec->size();i++){
 //         heepId = passHEEPId->at(index1[i]);
          mediumId = passMediumId->at(index1[i]);
          passId = mediumId; 
     //     passId = heepId;            
          if(passId) {
	    if(fabs(etaSC->at(index1[i])) < 2.5 && !(fabs(etaSC->at(index1[i])) > 1.4442 && fabs(etaSC->at(index1[i])) < 1.566)){
            newelePt.push_back(ptElec->at(index1[i]));
            neweleEta.push_back(etaElec->at(index1[i]));
            neweleEnr.push_back(energyElec->at(index1[i]));
            newelePhi.push_back(phiElec->at(index1[i]));
            neweleMediumId.push_back(passMediumId->at(index1[i]));
            newSCeta.push_back(etaSC->at(index1[i]));
          } //eta    
          } // id
          } //nEle
          } //size >=2

      if(newelePt.size()==2){
      int countTrig = 0;
      int countTrigDR3 = 0;
      for(unsigned int j = 0; j < pt_Ele23->size(); j++){
          dR1 = deltaR(neweleEta.at(0), newelePhi.at(0), eta_Ele23->at(j), phi_Ele23->at(j));
          dR2 = deltaR(neweleEta.at(1), newelePhi.at(1), eta_Ele23->at(j), phi_Ele23->at(j));
          if(dR1 < 0.1 || dR2 < 0.1) countTrig++;
	  if(dR1 < 0.3 || dR2 < 0.3) countTrigDR3++;
        } // trigger_object loop
          Bool_t isPassAcc = kFALSE;
          isPassAcc = isPassAccCondition_Electron(newelePt.at(0),newSCeta.at(0),newelePt.at(1),newSCeta.at(1));
	  ele1.SetPtEtaPhiE(newelePt.at(0),neweleEta.at(0),newelePhi.at(0),neweleEnr.at(0));
          ele2.SetPtEtaPhiE(newelePt.at(1),neweleEta.at(1),newelePhi.at(1),neweleEnr.at(1));
	  dielectron=ele1+ele2;
          recomass = dielectron.M();

//cout << "count trig = " << countTrig << "  " << "Ele23_WPLoose" << Ele23_WPLoose << endl;
//cout << "pu weight = " << puWeights << endl;
if(isPassAcc == kTRUE && recomass > 10. && countTrig != 0 && Ele23_WPLoose) {
Flag_PassEff = kTRUE;
} // acceptance on RECO

      if (isPassAcc == kTRUE && recomass > 10.) {
	double w=theWeight*lumiWeight;
	double wPU= w*puWeights;
	if (countTrig!=0 && Ele23_WPLoose) {
	  h1_eff_sumPass_chk->Fill(postFSRmass, wPU);
	  h1_eff_sumPass_noPU->Fill(postFSRmass, w);
	}
	if (Ele23_WPLoose) {
	  h1_eff_sumPass_noTrigObjMatching->Fill(postFSRmass, wPU);
	  if (countTrigDR3!=0)
	    h1_eff_sumPass_matchDR3->Fill(postFSRmass, wPU);
	}
	h1_eff_sumPass_noTrig->Fill(postFSRmass, wPU);
	h1_eff_sumPass_noTrigNoPU->Fill(postFSRmass, w);
      }
} //size == 2

if(Flag_PassEff == kTRUE){
h_mass_EffPass->Fill(postFSRmass,theWeight*lumiWeight*puWeights);
h_mass_EffTotal->Fill(postFSRmass,theWeight*lumiWeight);
 h1_eff_sumPass->Fill(postFSRmass,theWeight*lumiWeight*puWeights);
 h1_eff_sumTot->Fill(postFSRmass,theWeight*lumiWeight);
}
else {
h_mass_EffTotal->Fill(postFSRmass,theWeight*lumiWeight);
 h1_eff_sumTot->Fill(postFSRmass,theWeight*lumiWeight);
}

} //Acceptance Condition
}  // tauFlag
} //event loop
  file[j]->Write();
}//file loop


  TH1D* h1eff=(TH1D*)h1_eff_sumPass->Clone("h1eff");
  h1eff->SetTitle("Efficiency;M_{ee} [GeV];efficiency");
  h1eff->Divide(h1_eff_sumPass,h1_eff_sumTot,1,1,"B");
  TH1D* h1acc=(TH1D*)h1_acc_sumPass->Clone("h1acc");
  h1acc->SetTitle("Acceptance;M_{ee} [GeV];acceptance");
  h1acc->Divide(h1_acc_sumPass,h1_acc_sumTot,1,1,"B");
  TH1D* h1effAcc=(TH1D*)h1_eff_sumPass->Clone("h1effAcc");
  h1effAcc->SetTitle("Efficiency times Acceptance;M_{ee} [GeV];eff #times A");
  h1effAcc->Divide(h1_eff_sumPass,h1_acc_sumTot,1,1,"B");

  TFile fout("dyee_preFSR_forAccEff_v1steps.root","RECREATE");
  h1preFSR_1GeV->Write();
  h1preFSR->Write();
  h1_eff_sumPass->Write();
  h1_eff_sumTot->Write();
  h1_acc_sumPass->Write();
  h1_acc_sumTot->Write();
  h1eff->Write();
  h1acc->Write();
  h1effAcc->Write();

  h1_eff_sumPass_chk->Write();
  h1_eff_sumPass_noPU->Write();
  h1_eff_sumPass_noTrigObjMatching->Write();
  h1_eff_sumPass_matchDR3->Write();
  h1_eff_sumPass_noTrig->Write();
  h1_eff_sumPass_noTrigNoPU->Write();
  fout.Close();

}
