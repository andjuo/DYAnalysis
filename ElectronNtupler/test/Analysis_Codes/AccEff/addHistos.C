#include <TROOT.h>
#include <TH1D.h>
#include <TFile.h>
#include <iostream>

void addHistos() {

  int nMass=11;
  int mass[11] = {10,50,100,200,400,500,700,800,1000,1500,2000};

  TH1D *h1_mass_EffPass=NULL, *h1_mass_EffTotal=NULL;
  TH1D *h1_mass_AccPass=NULL, *h1_mass_AccTotal=NULL;

  for (int im=0; im<nMass; im++) {
    TString fname=Form("DY_%d_EE_76X_forAccEff.root",mass[im]);
    TFile fin(fname);
    if (!fin.IsOpen()) {
      std::cout << "failed to open the file <" << fin.GetName() << "\n";
      return;
    }

    TH1D *h1effPass=(TH1D*)fin.Get("h_mass_EffPass");
    TH1D *h1effTot =(TH1D*)fin.Get("h_mass_EffTotal");
    TH1D *h1accPass=(TH1D*)fin.Get("h_mass_AccPass");
    TH1D *h1accTot =(TH1D*)fin.Get("h_mass_AccTotal");
    if (!h1effPass || !h1effTot || !h1accPass || !h1accTot) {
      std::cout << "failed to get histo\n";
      fin.Close();
      return;
    }
    if (h1_mass_EffPass==NULL) {
      h1_mass_EffPass=(TH1D*)h1effPass->Clone("h1_mass_EffPass");
      h1_mass_EffPass->SetDirectory(0);
    }
    else {
      h1_mass_EffPass->Add(h1effPass);
    }
    if (h1_mass_EffTotal==NULL) {
      h1_mass_EffTotal=(TH1D*)h1effTot->Clone("h1_mass_EffTotal");
      h1_mass_EffTotal->SetDirectory(0);
    }
    else {
      h1_mass_EffTotal->Add(h1effTot);
    }
    if (h1_mass_AccPass==NULL) {
      h1_mass_AccPass=(TH1D*)h1accPass->Clone("h1_mass_AccPass");
      h1_mass_AccPass->SetDirectory(0);
    }
    else {
      h1_mass_AccPass->Add(h1accPass);
    }
    if (h1_mass_AccTotal==NULL) {
      h1_mass_AccTotal=(TH1D*)h1accTot->Clone("h1_mass_AccTotal");
      h1_mass_AccTotal->SetDirectory(0);
    }
    else {
      h1_mass_AccTotal->Add(h1accTot);
    }
    delete h1effPass;
    delete h1effTot;
    delete h1accPass;
    delete h1accTot;
    fin.Close();
  }

  TH1D* h1eff=(TH1D*)h1_mass_EffPass->Clone("h1eff");
  h1eff->Divide(h1_mass_EffPass,h1_mass_EffTotal,1,1,"B");
  TH1D *h1acc=(TH1D*)h1_mass_AccPass->Clone("h1acc");
  h1acc->Divide(h1_mass_AccPass,h1_mass_AccTotal,1,1,"B");

  TFile fout("dyee_addHistos_accEff.root","RECREATE");
  h1eff->Write();
  h1acc->Write();
  h1_mass_EffPass->Write();
  h1_mass_EffTotal->Write();
  h1_mass_AccPass->Write();
  h1_mass_AccTotal->Write();
  fout.Close();
}
