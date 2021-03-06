#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TLorentzVector.h>
#include <math.h>

void subtractEWK_dyee() {

  TString workdir;
  std::vector<TFile*> InputFiles_signal_DY;
  int mass[12] = {10,50,100,200,400,500,700,800,1000,1500,2000,3000};

  double xsec[11] = {18609.9/3,5789./3,226./3,7.68/3,0.423/3,0.24/3,0.036/3,0.03/3,0.0159/3,0.00201/3,0.00054/3};
  double sumofWts[11] = {758771667841.549683,150453917158.292908,219884886.214768,7008726.024242,122990.185607,69998.861766,11080.956088,9415.627747,4893.463037,664.033593,159.615701};

  workdir = "/tmp/rchawla/eos/cms/store/group/phys_higgs/cmshww/arun/DYAnalysis_76X_Calibrated/DY_Signal/";
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

  int nsample = InputFiles_signal_DY.size();
  TFile* file[11];

  for(unsigned int jentry = 1; jentry < 3; ++jentry) {
    TTree * T1 = (TTree*)InputFiles_signal_DY.at(jentry)->Get("ntupler/ElectronTree");

    TFile *f1 = TFile::Open("../../../dataPUDist.root");
    TFile *f2 = TFile::Open("../../../PileUp_MC.root");

    // data histogram 
    TH1F *DATA_puDist = (TH1F*)f1->Get("pileup");
    DATA_puDist->Scale(1/DATA_puDist->Integral());

    // mc histogram 
    TH1F *MC_puDist = (TH1F*)f2->Get("pileup_MC");
    TH1F *weights = (TH1F*)DATA_puDist->Clone("weights");
    weights->Divide(MC_puDist);

    vector<float>   *genPreFSR_Pt;
    vector<float>   *genPreFSR_Eta;
    vector<float>   *genPreFSR_Phi;
    vector<float>   *genPreFSR_En;
    vector<float>   *ptElec;
    vector<float>   *etaElec;
    vector<float>   *phiElec;
    vector<float>   *energyElec;
    vector<float>   *etaSC;
    vector<int>     *passMediumId;
    Int_t           tauFlag;
    Double_t        theWeight;
    Bool_t          Ele23_WPLoose;
    Int_t           nPUTrue;

    genPreFSR_Pt = 0;
    genPreFSR_Eta = 0;
    genPreFSR_Phi = 0;
    genPreFSR_En = 0;
    ptElec = 0;
    etaElec = 0;
    phiElec = 0;
    energyElec = 0;
    etaSC = 0;
    passMediumId = 0;

    T1->SetBranchStatus("*",0);
    T1->SetBranchStatus("genPreFSR_Pt", 1);
    T1->SetBranchStatus("genPreFSR_Eta", 1);
    T1->SetBranchStatus("genPreFSR_Phi", 1);
    T1->SetBranchStatus("genPreFSR_En", 1);
    T1->SetBranchStatus("tauFlag", 1);
    T1->SetBranchStatus("theWeight", 1);
    T1->SetBranchStatus("ptElec", 1);
    T1->SetBranchStatus("etaElec", 1);
    T1->SetBranchStatus("phiElec", 1);
    T1->SetBranchStatus("energyElec", 1);
    T1->SetBranchStatus("etaSC", 1);
    T1->SetBranchStatus("passMediumId", 1);
    T1->SetBranchStatus("Ele23_WPLoose", 1);
    T1->SetBranchStatus("nPUTrue", 1);

    T1->SetBranchAddress("genPreFSR_Pt", &genPreFSR_Pt);
    T1->SetBranchAddress("genPreFSR_Eta", &genPreFSR_Eta);
    T1->SetBranchAddress("genPreFSR_Phi", &genPreFSR_Phi);
    T1->SetBranchAddress("genPreFSR_En", &genPreFSR_En);
    T1->SetBranchAddress("ptElec", &ptElec);
    T1->SetBranchAddress("etaElec", &etaElec);
    T1->SetBranchAddress("phiElec", &phiElec);
    T1->SetBranchAddress("energyElec", &energyElec);
    T1->SetBranchAddress("etaSC", &etaSC);
    T1->SetBranchAddress("passMediumId", &passMediumId);
    T1->SetBranchAddress("tauFlag", &tauFlag);
    T1->SetBranchAddress("theWeight", &theWeight);
    T1->SetBranchAddress("Ele23_WPLoose", &Ele23_WPLoose);
    T1->SetBranchAddress("nPUTrue", &nPUTrue);

    file[jentry] = new TFile(Form("1Tight_1Loose/DYEE_M%dto%d.root",mass[jentry],mass[jentry+1]),"RECREATE");
    TTree *tree = new TTree("tree"," after preselections tree");

    int count, count1;
    int mediumId, dzMaskId, passId;
    bool passKin;
    bool ptcut;
    double massGen;
    TLorentzVector ele1,ele2,diReco;
    TLorentzVector gen1,gen2,diGen;

    // Branch variable declaration
    double ElePT_Fail, ElePT_Pass;
    double EleEta_Fail, EleEta_Pass;
    double wjetEstMass;
    double lumiWeight, genWeight, PUWeight;
    bool BB, BE, EE;

    // Branch declaration
    tree->Branch("ElePT_Fail", &ElePT_Fail, "ElePT_Fail/D");
    tree->Branch("EleEta_Fail", &EleEta_Fail, "EleEta_Fail/D");
    tree->Branch("ElePT_Pass", &ElePT_Pass, "ElePT_Pass/D");
    tree->Branch("EleEta_Pass", &EleEta_Pass, "EleEta_Pass/D");
    tree->Branch("wjetEstMass", &wjetEstMass, "wjetEstMass/D");
    tree->Branch("lumiWeight", &lumiWeight, "lumiWeight/D");
    tree->Branch("genWeight", &genWeight, "genWeight/D");
    tree->Branch("PUWeight", &PUWeight, "PUWeight/D");
    tree->Branch("BB", &BB, "BB/B");
    tree->Branch("BE", &BE, "BE/B");
    tree->Branch("EE", &EE, "EE/B");

    double lumi_Weight = xsec[jentry]/sumofWts[jentry];
    cout<<"DY Sample: "<<mass[jentry]<<"to"<<mass[jentry+1]<<endl;

    vector<int> idx; vector<int> idx1;

    int nentries = T1->GetEntries();
    //int nentries = 50000;
    cout<<"entries: "<<nentries<<endl;
    for (unsigned int i=0; i < nentries; i++) {
      T1->GetEntry(i);

      if(i%1000000 == 0){
	cout << "Events Processed :  " << i << endl;
      }

      int index[ptElec->size()];
      float pt[ptElec->size()];

      for(unsigned int el=0; el<ptElec->size(); el++) {
	pt[el]=ptElec->at(el); }

      int size = sizeof(pt)/sizeof(pt[0]);
      TMath::Sort(size,pt,index,true);

      count = 0; count1 = 0;
      mediumId = 0; dzMaskId = 0; passId =0;
      passKin = false;
      massGen = 0.0;
      ptcut = false;
      bool etacut = false;
      idx.clear();
      idx1.clear();
      BB=false; BE=false; EE=false;

      // PU Weight
      int bin = 0;
      bin = weights->GetXaxis()->FindBin(nPUTrue);
      double puweight = weights->GetBinContent(bin);

      if(genPreFSR_Pt->size() == 2){
	gen1.SetPtEtaPhiE(genPreFSR_Pt->at(0),genPreFSR_Eta->at(0),genPreFSR_Phi->at(0),genPreFSR_En->at(0));
	gen2.SetPtEtaPhiE(genPreFSR_Pt->at(1),genPreFSR_Eta->at(1),genPreFSR_Phi->at(1),genPreFSR_En->at(1));

	diGen=gen1+gen2;
	massGen=diGen.M();
      }

      if(!Ele23_WPLoose) continue;
      if(jentry==1 && massGen > 100.) continue;        // Gen Mass cut ----- for 50 to inf sample
      if(tauFlag) continue;               // taus

      for(int j=0;j<ptElec->size();j++){

	mediumId = passMediumId->at(index[j]);
	//passKin = (fabs(etaSC->at(index[j])) < 2.5 && !(fabs(etaSC->at(index[j])) > 1.4442 && fabs(etaSC->at(index[j])) < 1.566));

	if(mediumId){
	  count++;
	  idx.push_back(index[j]);
	}

	if(!mediumId){
	  count1++;
	  idx1.push_back(index[j]);
	}

      } //pt size > 2.

      if(count == 1 && count1 == 1){

	//cout<<"idx[0]: "<<idx[0]<<"   "<<ptElec->at(idx[0])<<"   "<<"idx1[0]: "<<idx1[0]<<"   "<<ptElec->at(idx1[0])<<endl;

	if((idx[0] < idx1[0]) && ptElec->at(idx[0]) > 30 && ptElec->at(idx1[0]) > 10.) ptcut = true;
	if((idx[0] > idx1[0]) && ptElec->at(idx1[0]) > 30 && ptElec->at(idx[0]) > 10.) ptcut = true;

	if((fabs(etaSC->at(idx[0])) < 2.5 && !(fabs(etaSC->at(idx[0])) > 1.4442 && fabs(etaSC->at(idx[0])) < 1.566)) && (fabs(etaSC->at(idx1[0])) < 2.5 && !(fabs(etaSC->at(idx1[0])) > 1.4442 && fabs(etaSC->at(idx1[0])) < 1.566))) etacut = true;

	if(ptcut && etacut) {

	  ElePT_Fail  = ptElec->at(idx1[0]);
	  EleEta_Fail = etaElec->at(idx1[0]);
	  ElePT_Pass  = ptElec->at(idx[0]);
	  EleEta_Pass = etaElec->at(idx[0]);

	  ele1.SetPtEtaPhiE(ptElec->at(idx[0]),etaElec->at(idx[0]),phiElec->at(idx[0]),energyElec->at(idx[0]));
	  ele2.SetPtEtaPhiE(ptElec->at(idx1[0]),etaElec->at(idx1[0]),phiElec->at(idx1[0]),energyElec->at(idx1[0]));

	  diReco=ele1+ele2;
	  wjetEstMass = diReco.M();

	  lumiWeight = lumi_Weight;
	  genWeight  = theWeight;
	  PUWeight = puweight;

	  if(fabs(etaSC->at(idx[0])) < 1.4442 && fabs(etaSC->at(idx1[0])) < 1.4442) BB = true;
	  if((fabs(etaSC->at(idx[0])) < 1.4442 && fabs(etaSC->at(idx1[0])) > 1.566 && fabs(etaSC->at(idx1[0])) < 2.5) || (fabs(etaSC->at(idx[0])) > 1.566 && fabs(etaSC->at(idx[0])) < 2.5 && fabs(etaSC->at(idx1[0])) < 1.4442)) BE = true;
	  if(fabs(etaSC->at(idx[0])) > 1.566 && fabs(etaSC->at(idx[0])) < 2.5 && fabs(etaSC->at(idx1[0])) > 1.566 && fabs(etaSC->at(idx1[0])) < 2.5) EE = true;

	  tree->Fill();
	}
      }// count==1
    } // event

    file[jentry]->Write();
    file[jentry]->Close();
  }
}
