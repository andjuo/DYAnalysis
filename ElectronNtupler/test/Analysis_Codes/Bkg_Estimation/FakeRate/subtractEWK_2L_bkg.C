#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TLorentzVector.h>
#include <math.h>

void subtractEWK_2L_bkg() {

  TString workdir;
  std::vector<TFile*> InputFiles_bkg;
  const char *bkg[6] = {"TTbar", "diBoson_WW", "diBoson_WZ", "diBoson_ZZ", "Single_antiTop", "SingleTop"};
  double xsec[6] = {831.76,118.7,66.1,15.4,35.6,35.6};
  double noEvts[6] = {85849690.,988418.,1000000.,985600.,999400.,1000000.};
  
  workdir = "/tmp/rchawla/eos/cms/store/group/phys_higgs/cmshww/arun/DYAnalysis_76X_Calibrated/Backgrounds/";
  //workdir = "/tmp/rchawla/";
  InputFiles_bkg.clear();

  InputFiles_bkg.push_back(TFile::Open(workdir+"TTbar.root"));
  InputFiles_bkg.push_back(TFile::Open(workdir+"diBoson_WW.root"));
  InputFiles_bkg.push_back(TFile::Open(workdir+"diBoson_WZ.root"));
  InputFiles_bkg.push_back(TFile::Open(workdir+"diBoson_ZZ.root"));
  InputFiles_bkg.push_back(TFile::Open(workdir+"Single_antiTop.root"));
  InputFiles_bkg.push_back(TFile::Open(workdir+"SingleTop.root"));

  int nsample = InputFiles_bkg.size();
  TFile* file[6];

  for(unsigned int jentry = 0; jentry < InputFiles_bkg.size(); ++jentry) {
    TTree * T1 = (TTree*)InputFiles_bkg.at(jentry)->Get("ntupler/ElectronTree");

    TFile *f1 = TFile::Open("../../../dataPUDist.root");
    TFile *f2 = TFile::Open("../../../PileUp_MC.root");

    // data histogram 
    TH1F *DATA_puDist = (TH1F*)f1->Get("pileup");
    DATA_puDist->Scale(1/DATA_puDist->Integral());

    // mc histogram 
    TH1F *MC_puDist = (TH1F*)f2->Get("pileup_MC");
    TH1F *weights = (TH1F*)DATA_puDist->Clone("weights");
    weights->Divide(MC_puDist);

    vector<float>   *ptElec;
    vector<float>   *etaElec;
    vector<float>   *phiElec;
    vector<float>   *energyElec;
    vector<float>   *etaSC;
    vector<int>     *passMediumId;
    Double_t        theWeight;
    Bool_t          Ele23_WPLoose;
    Int_t           nPUTrue;

    ptElec = 0;
    etaElec = 0;
    phiElec = 0;
    energyElec = 0;
    etaSC = 0;
    passMediumId = 0;

    T1->SetBranchStatus("*", 0);
    T1->SetBranchStatus("ptElec", 1);
    T1->SetBranchStatus("etaElec", 1);
    T1->SetBranchStatus("phiElec", 1);
    T1->SetBranchStatus("energyElec", 1);
    T1->SetBranchStatus("etaSC", 1);
    T1->SetBranchStatus("passMediumId", 1);
    T1->SetBranchStatus("theWeight", 1);
    T1->SetBranchStatus("Ele23_WPLoose", 1);
    T1->SetBranchStatus("nPUTrue", 1);

    T1->SetBranchAddress("ptElec", &ptElec);
    T1->SetBranchAddress("etaElec", &etaElec);
    T1->SetBranchAddress("phiElec", &phiElec);
    T1->SetBranchAddress("energyElec", &energyElec);
    T1->SetBranchAddress("etaSC", &etaSC);
    T1->SetBranchAddress("passMediumId", &passMediumId);
    T1->SetBranchAddress("theWeight", &theWeight);
    T1->SetBranchAddress("Ele23_WPLoose", &Ele23_WPLoose);
    T1->SetBranchAddress("nPUTrue", &nPUTrue);

    file[jentry] = new TFile(Form("2Loose/%s.root",bkg[jentry]),"RECREATE");
    TTree *tree = new TTree("tree"," after preselections tree");

    int count;
    int mediumId, dzMaskId, passId;
    bool passKin;
    TLorentzVector ele1,ele2,diReco;

    // Branch variable declaration
    double Ele1PT; double Ele2PT;
    double Ele1Eta; double Ele2Eta;
    double qcdEstMass;
    double lumiWeight, genWeight, PUWeight;
    bool BB, BE, EE;

    // Branch declaration
    tree->Branch("Ele1PT", &Ele1PT, "Ele1PT/D");
    tree->Branch("Ele2PT", &Ele2PT, "Ele2PT/D");
    tree->Branch("Ele1Eta", &Ele1Eta, "Ele1Eta/D");
    tree->Branch("Ele2Eta", &Ele2Eta, "Ele2Eta/D");
    tree->Branch("qcdEstMass", &qcdEstMass, "qcdEstMass/D");
    tree->Branch("lumiWeight", &lumiWeight, "lumiWeight/D");
    tree->Branch("genWeight", &genWeight, "genWeight/D");
    tree->Branch("PUWeight", &PUWeight, "PUWeight/D");
    tree->Branch("BB", &BB, "BB/B");
    tree->Branch("BE", &BE, "BE/B");
    tree->Branch("EE", &EE, "EE/B");

    double lumi_Weight = xsec[jentry]/noEvts[jentry];
    cout<<"Background Sample: "<<bkg[jentry]<<endl;

    vector<int> idx; 

    int nentries = T1->GetEntries();
    //int nentries = 50000;
    cout<<"entries: "<<nentries<<endl;
    for (unsigned int i=0; i < nentries; i++) {
      T1->GetEntry(i);

      if(i%1000000 == 0){
	cout << "Events Processed :  " << i << endl;
      }

      int index[ptElec->size()];
      float pt1[ptElec->size()];

      for(unsigned int el=0; el<ptElec->size(); el++) {
	pt1[el]=ptElec->at(el); }

      int size = sizeof(pt1)/sizeof(pt1[0]);
      TMath::Sort(size,pt1,index,true);

      count = 0;
      mediumId = 0; dzMaskId = 0; passId =0;
      passKin = false;
      bool etacut = false;
      idx.clear();
      BB=false; BE=false; EE=false;

      // PU Weight
      int bin = 0;
      bin = weights->GetXaxis()->FindBin(nPUTrue);
      double puweight = weights->GetBinContent(bin);

      if(!Ele23_WPLoose) continue;

      for(int j=0;j<ptElec->size();j++){

	mediumId = passMediumId->at(index[j]);
	//passKin = (fabs(etaSC->at(index[j])) < 2.5 && !(fabs(etaSC->at(index[j])) > 1.4442 && fabs(etaSC->at(index[j])) < 1.566));

	if(!mediumId){
	  count++;
	  idx.push_back(index[j]);
	}

      } //pt size > 2.

      if(count == 2){

	if((fabs(etaSC->at(idx[0])) < 2.5 && !(fabs(etaSC->at(idx[0])) > 1.4442 && fabs(etaSC->at(idx[0])) < 1.566)) && (fabs(etaSC->at(idx[1])) < 2.5 && !(fabs(etaSC->at(idx[1])) > 1.4442 && fabs(etaSC->at(idx[1])) < 1.566))) etacut = true;
	
	if(ptElec->at(idx[0]) > 30 && ptElec->at(idx[1]) > 10 && etacut){
	  Ele1PT  = ptElec->at(idx[0]);
	  Ele1Eta = etaElec->at(idx[0]);
	  Ele2PT  = ptElec->at(idx[1]);
	  Ele2Eta = etaElec->at(idx[1]);

	  ele1.SetPtEtaPhiE(ptElec->at(idx[0]),etaElec->at(idx[0]),phiElec->at(idx[0]),energyElec->at(idx[0]));
	  ele2.SetPtEtaPhiE(ptElec->at(idx[1]),etaElec->at(idx[1]),phiElec->at(idx[1]),energyElec->at(idx[1]));

	  diReco=ele1+ele2;
	  qcdEstMass = diReco.M();

	  lumiWeight = lumi_Weight;
	  genWeight  = theWeight;
	  PUWeight = puweight;

	  if(fabs(etaSC->at(idx[0])) < 1.4442 && fabs(etaSC->at(idx[1])) < 1.4442) BB = true;
	  if((fabs(etaSC->at(idx[0])) < 1.4442 && fabs(etaSC->at(idx[1])) > 1.566 && fabs(etaSC->at(idx[1])) < 2.5) || (fabs(etaSC->at(idx[0])) > 1.566 && fabs(etaSC->at(idx[0])) < 2.5 && fabs(etaSC->at(idx[1])) < 1.4442)) BE = true;
	  if(fabs(etaSC->at(idx[0])) > 1.566 && fabs(etaSC->at(idx[0])) < 2.5 && fabs(etaSC->at(idx[1])) > 1.566 && fabs(etaSC->at(idx[1])) < 2.5) EE = true;

	  tree->Fill();

	}
      }// count==1
    } // event

    file[jentry]->Write();
    file[jentry]->Close();
  } // file Loop
}
