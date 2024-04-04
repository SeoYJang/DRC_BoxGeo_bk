#include "RootInterface.h"
#include "RecoInterface.h"
#include "DRsimInterface.h"
#include "functions.h"

#include "TROOT.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TPaveStats.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TGraph.h"
#include "TGraphErrors.h"

#include <iostream>
#include <string>
#include "Riostream.h"

int main(int argc, char* argv[]) {
  TString filename = argv[1];
  float low = std::stof(argv[2]);
  float high = std::stof(argv[3]);
  float en_val = std::stof(argv[4]);
  TString material = argv[5];
  TString particle = argv[6];

  int mNum = 24;

  gStyle->SetOptFit(1);

  RootInterface<DRsimInterface::DRsimEventData>* drInterface = new RootInterface<DRsimInterface::DRsimEventData>(std::string(filename)+".root", true);
  drInterface->set("DRsim","DRsimEventData");
  //drInterface->GetChain("DRsim");

  TH1F* tEdep = new TH1F("totEdep",";MeV;Evt",100,low*1000.,high*1000.);
  tEdep->Sumw2(); tEdep->SetLineColor(kRed); tEdep->SetLineWidth(2);
  TH1F* tHit_C = new TH1F("Hit_C",";# of p.e.;Evt",200,0,3000);
  tHit_C->Sumw2(); tHit_C->SetLineColor(kBlue); tHit_C->SetLineWidth(2);
  TH1F* tHit_S = new TH1F("Hit_S",";# of p.e.;Evt",200,0,40000);
  tHit_S->Sumw2(); tHit_S->SetLineColor(kRed); tHit_S->SetLineWidth(2);
  TH1F* tE_C = new TH1F("E_C",";GeV;Evt",100,low,high);
  tE_C->Sumw2(); tE_C->SetLineColor(kBlue); tE_C->SetLineWidth(2);
  TH1F* tE_S = new TH1F("E_S",";GeV;Evt",100,low,high);
  tE_S->Sumw2(); tE_S->SetLineColor(kRed); tE_S->SetLineWidth(2);
  TH1F* tE_Scorr = new TH1F("E_Scorr",";GeV;Evt",100,low,high);
  tE_Scorr->Sumw2(); tE_Scorr->SetLineColor(kRed); tE_Scorr->SetLineWidth(2);
  TH1F* tE_SC = new TH1F("E_SC",";GeV;Evt",100,2.*low,2.*high);
  tE_SC->Sumw2(); tE_SC->SetLineColor(kBlack); tE_SC->SetLineWidth(2);
  TH1F* tE_DR = new TH1F("E_DR",";Evt",100,low,high);
  tE_DR->Sumw2(); tE_DR->SetLineColor(kBlack); tE_DR->SetLineWidth(2);
  TH1F* tE_DRcorr = new TH1F("E_DRcorr",";GeV;Evt",100,low,high);
  tE_DRcorr->Sumw2(); tE_DRcorr->SetLineColor(kBlack); tE_DRcorr->SetLineWidth(2);
  TH1F* tP_leak = new TH1F("Pleak",";MeV;Evt",100,0.,1000.*high);
  tP_leak->Sumw2(); tP_leak->SetLineWidth(2);
  TH1F* tP_leak_nu = new TH1F("Pleak_nu",";MeV;Evt",100,0.,1000.*high);
  tP_leak_nu->Sumw2(); tP_leak_nu->SetLineWidth(2);
  TH1F* tDepth = new TH1F("depth",";m;Evt",70,-1.0,3.0);
  tDepth->Sumw2(); tDepth->SetLineWidth(2);
  TH1F* tCDepth = new TH1F("Cdepth",";m;Evt",70,-1.0,3.0);
  tCDepth->Sumw2(); tCDepth->SetLineWidth(2);

  TH1F* tP_leak_theta = new TH1F("Pleak_theta",";degree;Evt",180,-180.,180.);
  tP_leak_theta->Sumw2(); tP_leak_theta->SetLineWidth(2); tP_leak_theta->SetStats(0);
  TH1F* tP_leak_phi = new TH1F("Pleak_phi",";degree;Evt",180,-180.,180.);
  tP_leak_phi->Sumw2(); tP_leak_phi->SetLineWidth(2); tP_leak_phi->SetStats(0);

  TH2F* tScvsC = new TH2F("SvsC", ";E_{Scorr};E_{C}", 100, low, high, 100, low, high);
  tScvsC->Sumw2(); tScvsC->SetStats(0);

  TH1F* tT_C = new TH1F("time_C",";ns;p.e.",700,0.,70.);
  tT_C->Sumw2(); tT_C->SetLineColor(kBlue); tT_C->SetLineWidth(2);
  TH1F* tT_S = new TH1F("time_S",";ns;p.e.",700,0.,70.);
  tT_S->Sumw2(); tT_S->SetLineColor(kRed); tT_S->SetLineWidth(2);
  TH1F* tWav_S = new TH1F("wavlen_S",";nm;p.e.",120,300.,900.);
  tWav_S->Sumw2(); tWav_S->SetLineColor(kRed); tWav_S->SetLineWidth(2);
  TH1F* tWav_C = new TH1F("wavlen_C",";nm;p.e.",120,300.,900.);
  tWav_C->Sumw2(); tWav_C->SetLineColor(kBlue); tWav_C->SetLineWidth(2);
  TH1F* tNhit_S = new TH1F("nHits_S",";n",200,0.,200.);
  tNhit_S->Sumw2(); tNhit_S->SetLineColor(kRed); tNhit_S->SetLineWidth(2);
  TH1F* tNhit_C = new TH1F("nHits_C",";n",50,0.,50.);
  tNhit_C->Sumw2(); tNhit_C->SetLineColor(kBlue); tNhit_C->SetLineWidth(2);
  TH1F* oEdep = new TH1F("oEdep",";MeV;Evt",100,low*1000.,high*1000.);
  oEdep->Sumw2(); oEdep->SetLineColor(kBlack); oEdep->SetLineWidth(2);
  TH1I* Chit = new TH1I("C_Hit","",100,0.,3000.);
  Chit->Sumw2(); Chit->SetLineColor(kBlue); Chit->SetLineWidth(2);
  TH1I* Shit = new TH1I("S_Hit","",100,0.,40000.);
  Shit->Sumw2(); Shit->SetLineColor(kRed); Shit->SetLineWidth(2);

  TH2F* t2DDES = new TH2F("2D Depth Scint", "", 70, -1.0, 3.0, 50, low, high); t2DDES->Sumw2(); t2DDES->SetStats(0);
  TH2F* t2DDEC = new TH2F("2D Depth Ceren", "", 70, -1.0, 3.0, 50, low, high); t2DDEC->Sumw2(); t2DDEC->SetStats(0);

  TH2D* t2DhitC = new TH2D("2D Hit C", "", 420, -0.5, 419.5, 420, -0.5, 419.5); t2DhitC->Sumw2(); t2DhitC->SetStats(0);
  TH2D* t2DhitS = new TH2D("2D Hit S", "", 420, -0.5, 419.5, 420, -0.5, 419.5); t2DhitS->Sumw2(); t2DhitS->SetStats(0);

  std::string filename_var = "/u/user/syjang/scratch/DRC_generic/results/ele/Variables.csv";
  std::ifstream in;
  std::string mat;
  float fCalibC, fCalibS, Cs, Ss, es, al, mtmax, rl, ch, Cpi;
  float Cthres, Cth, chi, Cpithres;

  float effspeed, attlen, matTmax, radilen;

  in.open(filename_var, std::ios::in);
  while (true) {
    in >> mat >> Cth >> Cs >> Ss >> es >> al >> mtmax >> rl >> ch >> Cpi;
    if (mat == material) {
      fCalibC = Cs;
      fCalibS = Ss;
      Cthres = Cth;
      effspeed = es;
      attlen = al;
      matTmax = mtmax;
      radilen = rl;
      chi = ch;
      Cpithres = Cpi;
    }
    if (!in.good()) break;
  }

  if (particle == "pi") Cthres = Cpithres;

  std::vector<float> E_Ss,E_Cs,E_Sscorr;

  unsigned int entries = drInterface->entries();
  if (entries > 3000) entries = 3000;
  while (drInterface->numEvt() < entries) {
    if (drInterface->numEvt() % 100 == 0) printf("Analyzing %dth event ...\n", drInterface->numEvt());

    DRsimInterface::DRsimEventData drEvt;
    drInterface->read(drEvt);

    float Edep = 0.; float fEdep = 0.;
    int moduleNum = 0;

    float rE_C, rE_S;
    float sE_C, sE_S;
    sE_C = 0.; sE_S = 0.;

    for (auto edepItr = drEvt.Edeps.begin(); edepItr != drEvt.Edeps.end(); ++edepItr) {
      auto edep = *edepItr;
      Edep += edep.Edep;
      moduleNum = edep.ModuleNum;
      if (moduleNum==mNum) fEdep += edep.Edep;
    }
    tEdep->Fill(Edep);
    oEdep->Fill(fEdep);

    TH1F* tT_max = new TH1F("tmax","",600,10.,70.);
    TH1F* tC_max = new TH1F("cmax","",600,10.,70.);

    float Pleak = 0.;
    float Eleak_nu = 0.;
    float Pleak_theta = 0., Pleak_phi = 0.;
    for (auto leak : drEvt.leaks) {
      TLorentzVector leak4vec;
      leak4vec.SetPxPyPzE(leak.px,leak.py,leak.pz,leak.E);
      if ( std::abs(leak.pdgId)==12 || std::abs(leak.pdgId)==14 || std::abs(leak.pdgId)==16 ) {
        Eleak_nu += leak4vec.P();
      } else {
        Pleak += leak4vec.P();
        if (leak4vec.Pz() >= 0) {
          Pleak_theta = TMath::ATan(leak4vec.Py()/leak4vec.Pz())*180/TMath::Pi();
          Pleak_phi   = TMath::ATan(leak4vec.Px()/leak4vec.Pz())*180/TMath::Pi();
        } else {
          Pleak_theta = (leak4vec.Py() >= 0) ? (TMath::ATan(leak4vec.Py()/leak4vec.Pz())*180/TMath::Pi()) + 180 : (TMath::ATan(leak4vec.Py()/leak4vec.Pz())*180/TMath::Pi() - 180);
          Pleak_phi   = (leak4vec.Px() >= 0) ? (TMath::ATan(leak4vec.Px()/leak4vec.Pz())*180/TMath::Pi()) + 180 : (TMath::ATan(leak4vec.Px()/leak4vec.Pz())*180/TMath::Pi() - 180);
        }
        tP_leak_theta->Fill(Pleak_theta);
        tP_leak_phi->Fill(Pleak_phi);
      }
    }
    tP_leak->Fill(Pleak);
    tP_leak_nu->Fill(Eleak_nu);

    int nHitC = 0; int nHitS = 0;
    int fC_hits = 0; int fS_hits = 0;
    int sum;

    for (auto tower = drEvt.towers.begin(); tower != drEvt.towers.end(); ++tower) {

      moduleNum = tower->ModuleNum;
      rE_C = 0.;
      rE_S = 0.;

      for (auto sipm = tower->SiPMs.begin(); sipm != tower->SiPMs.end(); ++sipm) {
        int plateNum = sipm->x; int fiberNum = sipm->y; 
        if ( RecoInterface::IsCerenkov(sipm->x,sipm->y) ) {
          tNhit_C->Fill(sipm->count);
          sum = 0;
          for (const auto timepair : sipm->timeStruct) {
            tT_C->Fill(timepair.first.first+0.05,timepair.second);
            tC_max->Fill(timepair.first.first+0.05,timepair.second);
            if (timepair.first.first < Cthres) {
              nHitC += timepair.second;
              sum += timepair.second;
              if (moduleNum==mNum) {
                fC_hits += timepair.second;
              }
              t2DhitC->Fill(60*(moduleNum%7)+fiberNum, 60*(moduleNum/7)+plateNum, timepair.second);
            }
          }

          for (const auto wavpair : sipm->wavlenSpectrum) {
            tWav_C->Fill(wavpair.first.first,wavpair.second);
          }

          rE_C += (float)sum / fCalibC;

        } else {
          tNhit_S->Fill(sipm->count);
          nHitS += sipm->count;
          rE_S += (float)sipm->count / fCalibS;
          if (moduleNum==mNum) {
            fS_hits += sipm->count;
          }
          t2DhitS->Fill(60*(moduleNum%7)+fiberNum, 60*(moduleNum/7)+plateNum, sipm->count);
          for (const auto timepair : sipm->timeStruct) {
            tT_S->Fill(timepair.first.first+0.05,timepair.second);
            tT_max->Fill(timepair.first.first+0.05,timepair.second);
          }
          for (const auto wavpair : sipm->wavlenSpectrum) {
            tWav_S->Fill(wavpair.first.first,wavpair.second);
          }
        }
      }
      sE_C += rE_C;
      sE_S += rE_S;
    }

    float T_max = tT_max->GetBinCenter( tT_max->GetMaximumBin() );
    float depth = (matTmax + radilen*(1./effspeed - 1./0.3) - T_max)/(1./effspeed - 1./0.3);

    float C_max = tC_max->GetBinCenter( tC_max->GetMaximumBin() );
    float Cdepth = (17.95 + radilen*(1.76) - C_max)/(1.76);

    float sE_Scorr = sE_S*std::exp((radilen-depth)/attlen);

    delete tT_max;
    delete tC_max;

    //int sum;

    /*if (C_max >= 17.5) {
      for (auto tower = drEvt.towers.begin(); tower != drEvt.towers.end(); ++tower) {
        rE_C = 0;

        for (auto sipm = tower->SiPMs.begin(); sipm != tower->SiPMs.end(); ++sipm) {
          sum = 0;
          int plateNum = sipm->x; int fiberNum = sipm->y; 
          if ( RecoInterface::IsCerenkov(sipm->x,sipm->y) ) {
            for (const auto timepair : sipm->timeStruct) {
              tT_C->Fill(timepair.first.first+0.05,timepair.second);
              if (timepair.first.first < 32.5) {
                sum += timepair.second;
              }
            }
          }
          rE_C += (float)sum / fCalibC;
        }
        sE_C += rE_C;
      }
    } else {
      for (auto tower = drEvt.towers.begin(); tower != drEvt.towers.end(); ++tower) {
        rE_C = 0;

        for (auto sipm = tower->SiPMs.begin(); sipm != tower->SiPMs.end(); ++sipm) {
          sum = 0;
          int plateNum = sipm->x; int fiberNum = sipm->y; 
          if ( RecoInterface::IsCerenkov(sipm->x,sipm->y) ) {
            for (const auto timepair : sipm->timeStruct) {
              //tT_C->Fill(timepair.first.first+0.05,timepair.second);
              if (timepair.first.first < 29.5) {
                sum += timepair.second;
              }
            }
          }
          rE_C += (float)sum / fCalibC;
        }
        sE_C += rE_C;
      }
    }*/
    
    E_Cs.push_back(sE_C);
    E_Ss.push_back(sE_S);
    E_Sscorr.push_back(sE_Scorr);

    tHit_C->Fill(nHitC);
    tHit_S->Fill(nHitS);
    Chit->Fill(fC_hits);
    Shit->Fill(fS_hits);

    tE_C->Fill(sE_C);
    tE_S->Fill(sE_S);
    tE_Scorr->Fill(sE_Scorr);
    tE_SC->Fill(sE_C + sE_S);

    tScvsC->Fill(sE_Scorr, sE_C);

    tE_DR->Fill(functions::E_DR(sE_C,sE_S,chi));
    tE_DRcorr->Fill(functions::E_DR(sE_C,sE_Scorr,chi));
    tDepth->Fill(depth);
    tCDepth->Fill(Cdepth);

    t2DDES->Fill(depth, sE_S);
    t2DDEC->Fill(Cdepth, sE_C);
  } // event loop
  //drInterface->close();
  //std::cout << "changed!" << std::endl;

  /*float maxtim = 9999.;
  for (int i=0; i<tT_C->GetNbinsX(); i++) {
    if (maxtim < tT_C->GetBinContent(i+1) && tT_C->GetBinContent(i+1) < tT_C->GetBinContent(i+2)) {
      std::cout << "Minimum C Timing : " << tT_C->GetBinCenter(i+1) << std::endl;
      break;
    }
    maxtim = tT_C->GetBinContent(i+1);
  }*/

  TCanvas* c = new TCanvas("c","");

  c->cd();

  //tT_C->Draw("Hist"); c->SaveAs(filename+"_Ctime_fortest.png");
  tScvsC->Draw("COLZ"); c->SaveAs(filename+"_2DScorrC.pdf");

  float TmaxC = tT_C->GetBinCenter(tT_C->GetMaximumBin());
  std::cout<< "T_max C  " << TmaxC << std::endl;

  float timi = 0.;
  float deltimi = tT_S->GetMaximum();
  int tii = 0.;
  while (true) {
    timi = tT_S->GetBinContent(tT_S->GetMaximumBin()-tii);
    if (deltimi <= std::abs(timi-(tT_S->GetMaximum()/10.))) break;
    deltimi = timi-(tT_S->GetMaximum()/10.);
    tii += 1;
  }
  std::cout<< "T_lead10%  " << tT_S->GetBinCenter(tT_S->GetMaximumBin()-tii) << std::endl;

  c->cd();
  t2DDES->Draw("COL"); t2DDES->SetMarkerStyle(kFullDotMedium);
  t2DDES->GetXaxis()->SetRangeUser(-0.5, 2.5);
  c->SaveAs(filename+"_2DdepthScint.pdf");

  c->cd();
  t2DDEC->Draw("COL"); t2DDEC->SetMarkerStyle(kFullDotMedium);
  t2DDEC->GetXaxis()->SetRangeUser(-0.5, 2.5);
  c->SaveAs(filename+"_2DdepthCeren.pdf");

  int t2Dxbins = t2DDES->GetNbinsX();
  float tdepval[t2Dxbins], tdeperr[t2Dxbins];
  float tDEx[t2Dxbins], tDExerr[t2Dxbins];
  int tcount = 0;
  int SgraphN = 0;

  TGraphErrors* tDEgraph = new TGraphErrors();
  
  for (int i=0; i<t2DDES->GetNbinsX(); i++) {
    TH1F* tdepth_bin = new TH1F("tdepth","",t2DDES->GetNbinsY(),low,high);
    for (int j=0; j<t2DDES->GetNbinsY(); j++) {
      if (t2DDES->GetBinContent(i+1, j+1) != 0) {
        tdepth_bin->Fill(t2DDES->GetYaxis()->GetBinCenter(j+1),t2DDES->GetBinContent(i+1, j+1));
        tcount += t2DDES->GetBinContent(i+1, j+1);
      }
    }
    if (tcount != 0) {
      tdepval[SgraphN] = tdepth_bin->GetMean();
      tdeperr[SgraphN] = tdepth_bin->GetRMS()/sqrt(tcount);
    //  std::cout << "Bin Num " << i << "  Mean Error : "<< tdepth_bin->GetMean() << " and " << tdepth_bin->GetRMS()/sqrt(tcount) << std::endl;
    } else {
      delete tdepth_bin;
      continue;
    }
    tDEx[SgraphN] = t2DDES->GetXaxis()->GetBinCenter(i+1);
    tDExerr[SgraphN] = 0.0285;
    tDEgraph->SetPoint(SgraphN, tDEx[SgraphN], tdepval[SgraphN]);
    tDEgraph->SetPointError(SgraphN, tDExerr[SgraphN], tdeperr[SgraphN]);
    SgraphN++;
    tcount = 0;
    delete tdepth_bin;
  }

  TF1 *f1 = new TF1("f1", "exp([0] + [1]*x)");
  tDEgraph->Fit(f1, "", "", 0., 1.3);
  c->cd();
  tDEgraph->Draw("AP");
  tDEgraph->SetTitle(""); 
  tDEgraph->GetXaxis()->SetRangeUser(-0.5, 2.5);
  tDEgraph->GetYaxis()->SetRangeUser(low, high); 

  c->SaveAs(filename+"_AttenLengthplotS.pdf");

  int t2Cxbins = t2DDEC->GetNbinsX();
  float tCdepval[t2Cxbins], tCdeperr[t2Cxbins];
  float tCDEx[t2Cxbins], tCDExerr[t2Cxbins];
  int tCount = 0;
  int CgraphN = 0;

  TGraphErrors* tCEgraph = new TGraphErrors();
  
  for (int i=0; i<t2DDEC->GetNbinsX(); i++) {
    TH1F* tdepth_bin = new TH1F("tdepth","",t2DDEC->GetNbinsY(),low,high);
    for (int j=0; j<t2DDEC->GetNbinsY(); j++) {
      if (t2DDEC->GetBinContent(i+1, j+1) != 0) {
        tdepth_bin->Fill(t2DDEC->GetYaxis()->GetBinCenter(j+1),t2DDEC->GetBinContent(i+1, j+1));
        tCount += t2DDEC->GetBinContent(i+1, j+1);
      }
    }
    if (tCount != 0) {
      tCdepval[CgraphN] = tdepth_bin->GetMean();
      tCdeperr[CgraphN] = tdepth_bin->GetRMS()/sqrt(tCount);
    //  std::cout << "Bin Num " << i << "  Mean Error : "<< tdepth_bin->GetMean() << " and " << tdepth_bin->GetRMS()/sqrt(tcount) << std::endl;
    } else {
      delete tdepth_bin;
      continue;
    }
    tCDEx[CgraphN] = t2DDEC->GetXaxis()->GetBinCenter(i+1);
    tCDExerr[CgraphN] = 0.0285;
    tCEgraph->SetPoint(CgraphN, tCDEx[CgraphN], tCdepval[CgraphN]);
    tCEgraph->SetPointError(CgraphN, tCDExerr[CgraphN], tCdeperr[CgraphN]);
    CgraphN += 1;
    tCount = 0;
    delete tdepth_bin;
  }

  TF1 *f2 = new TF1("f2", "exp([0] + [1]*x)");
  f2->SetParLimits(1, 0.00, 0.30);
  f2->SetLineColor(kRed);
  tCEgraph->Fit(f2, "MR", "", 0.0, 1.3);
  c->cd();
  tCEgraph->Draw("AP");
  tCEgraph->SetTitle(""); 
  tCEgraph->GetXaxis()->SetRangeUser(-0.5, 2.5);
  tCEgraph->GetYaxis()->SetRangeUser(low, high); 

  c->SaveAs(filename+"_AttenLengthplotC.pdf");

  int t2Dscbins = tScvsC->GetNbinsX();
  float tCval[t2Dscbins], tCerr[t2Dscbins];
  float tSval[t2Dscbins], tSerr[t2Dscbins];
  int tSvsCcount = 0;
  int tSvsCgraphN = 0;

  TGraphErrors* tSCgraph = new TGraphErrors();
  
  for (int i=0; i<t2Dscbins; i++) {
    TH1F* tSC_bin = new TH1F("tSCbin","",tScvsC->GetNbinsY(),low,high);
    for (int j=0; j<tScvsC->GetNbinsY(); j++) {
      if (tScvsC->GetBinContent(i+1, j+1) != 0) {
        tSC_bin->Fill(tScvsC->GetYaxis()->GetBinCenter(j+1),tScvsC->GetBinContent(i+1, j+1));
        tSvsCcount += tScvsC->GetBinContent(i+1, j+1);
      }
    }
    if (tSvsCcount != 0) {
      tCval[tSvsCgraphN] = tSC_bin->GetMean();
      tCerr[tSvsCgraphN] = tSC_bin->GetRMS()/sqrt(tSvsCcount);
    //  std::cout << "Bin Num " << i << "  Mean Error : "<< tdepth_bin->GetMean() << " and " << tdepth_bin->GetRMS()/sqrt(tcount) << std::endl;
    } else {
      delete tSC_bin;
      continue;
    }
    tSval[tSvsCgraphN] = tScvsC->GetXaxis()->GetBinCenter(i+1);
    tSerr[tSvsCgraphN] = (high-low)/(100*2);
    tSCgraph->SetPoint(tSvsCgraphN, tSval[tSvsCgraphN], tCval[tSvsCgraphN]);
    tSCgraph->SetPointError(tSvsCgraphN, tSerr[tSvsCgraphN], tCerr[tSvsCgraphN]);
    tSvsCgraphN++;
    tSvsCcount = 0;
    delete tSC_bin;
  }

  TF1 *f0 = new TF1("f0", "x*[0] + [1]");
  tSCgraph->Fit(f0, "MR", "", 16, 20);
  c->cd();
  tSCgraph->Draw("AP");
  tSCgraph->GetXaxis()->SetLimits(low, high);
  tSCgraph->SetMinimum(low);
  tSCgraph->SetMaximum(high);
  tSCgraph->SetTitle(""); 

  c->SaveAs(filename+"_SvsCfit.pdf");

  tEdep->Draw("Hist"); c->SaveAs(filename+"_Edep.png");

  c->cd();
  tE_S->SetTitle("");
  tE_S->Draw("Hist"); c->Update();
  TPaveStats* statsE_S = (TPaveStats*)c->GetPrimitive("stats");
  statsE_S->SetName("Scint");
  statsE_S->SetTextColor(kRed);
  statsE_S->SetY1NDC(.6); statsE_S->SetY2NDC(.8);

  tE_C->Draw("Hist&sames"); c->Update();
  TPaveStats* statsE_C = (TPaveStats*)c->GetPrimitive("stats");
  statsE_C->SetName("Cerenkov");
  statsE_C->SetTextColor(kBlue);
  statsE_C->SetY1NDC(.8); statsE_C->SetY2NDC(1.);
  c->SaveAs(filename+"_EcsHist.pdf");

  c->cd();
  tE_Scorr->SetTitle("");
  tE_Scorr->Draw("Hist"); c->Update();
  TPaveStats* statsE_Scorr = (TPaveStats*)c->GetPrimitive("stats");
  statsE_Scorr->SetName("Scint_corr");
  statsE_Scorr->SetTextColor(kRed);
  statsE_Scorr->SetY1NDC(.6); statsE_Scorr->SetY2NDC(.8);

  tE_C->Draw("Hist&sames"); c->Update();
  c->SaveAs(filename+"_EcsCorrHist.pdf");

  TF1* grE_C = new TF1("Cfit","gaus",low,high); grE_C->SetLineColor(kBlue);
  TF1* grE_S = new TF1("Sfit","gaus",low,high); grE_S->SetLineColor(kRed);
  tE_C->SetOption("p"); tE_C->Fit(grE_C,"R+&same");
  tE_S->SetOption("p"); tE_S->Fit(grE_S,"R+&same");

  c->cd();
  tE_S->SetTitle("");
  tE_S->Draw(""); c->Update();
  statsE_S->SetName("Scint");
  statsE_S->SetTextColor(kRed);
  statsE_S->SetX1NDC(.7);
  statsE_S->SetY1NDC(.4); statsE_S->SetY2NDC(.7);

  tE_C->Draw("sames"); c->Update();
  statsE_C->SetName("Cerenkov");
  statsE_C->SetTextColor(kBlue);
  statsE_C->SetX1NDC(.7);
  statsE_C->SetY1NDC(.7); statsE_C->SetY2NDC(1.);

  c->SaveAs(filename+"_Ecs.pdf");

  c->cd();
  tE_SC->Draw("Hist"); c->Update();
  c->SaveAs(filename+"_EsumHist.pdf");
  tE_DR->Draw("Hist"); c->Update();
  c->SaveAs(filename+"_EcorrHist.pdf");
  
  tE_DRcorr->Draw("Hist"); c->Update();
  c->SaveAs(filename+"_EcorrCorrHist.pdf");

  TF1* grE_SC = new TF1("S+Cfit","gaus",2.*low,2.*high); grE_SC->SetLineColor(kBlack);
  TF1* grE_DR = new TF1("Efit","gaus",low,high); grE_DR->SetLineColor(kBlack);
  tE_SC->SetOption("p"); tE_SC->Fit(grE_SC,"R+&same");
  tE_DRcorr->SetOption("p"); tE_DRcorr->Fit(grE_DR,"R+&same");

  tE_SC->Draw(""); c->SaveAs(filename+"_Esum.pdf");
  tE_DRcorr->Draw(""); c->SaveAs(filename+"_EDRcorr.pdf");

  tDepth->Draw("Hist"); tDepth->GetXaxis()->SetRangeUser(-0.5, 2.5); 
  c->SaveAs(filename+"_depth.png");

  tCDepth->Draw("Hist"); tCDepth->GetXaxis()->SetRangeUser(-0.5, 2.5); 
  c->SaveAs(filename+"_Cdepth.png");

  c->SetLogy(1);
  tP_leak->Draw("Hist"); c->SaveAs(filename+"_Pleak.png");
  tP_leak_nu->Draw("Hist"); c->SaveAs(filename+"_Pleak_nu.png");
  c->SetLogy(0);

  tP_leak_theta->Draw("Hist"); c->SaveAs(filename+"_Pleak_theta.png");
  tP_leak_phi->Draw("Hist"); c->SaveAs(filename+"_Pleak_phi.png");

  tHit_C->Draw("Hist"); c->SaveAs(filename+"_nHitpEventC.png");
  tHit_S->Draw("Hist"); c->SaveAs(filename+"_nHitpEventS.png");
  Chit->Draw(); c->SaveAs(filename+"_Chit.png");
  Shit->Draw(); c->SaveAs(filename+"_Shit.png");

  t2DhitS->Draw("COLZ"); c->SaveAs(filename+"_n2DHitS.png");
  t2DhitC->Draw("COLZ"); c->SaveAs(filename+"_n2DHitC.png");

  TGraph* grSvsC = new TGraph(entries,&(E_Ss[0]),&(E_Cs[0]));
  grSvsC->SetTitle("");
  grSvsC->SetMarkerSize(0.5); grSvsC->SetMarkerStyle(20);
  grSvsC->GetXaxis()->SetLimits(0.,high);
  grSvsC->GetYaxis()->SetRangeUser(0.,high);
  grSvsC->SetMaximum(high);
  grSvsC->SetMinimum(0.);
  grSvsC->Draw("ap");
  c->SaveAs(filename+"_SvsC.png");

  TGraph* grScorrvsC = new TGraph(entries,&(E_Sscorr[0]),&(E_Cs[0]));
  grScorrvsC->SetTitle("");
  grScorrvsC->SetMarkerSize(0.5); grScorrvsC->SetMarkerStyle(20);
  grScorrvsC->GetXaxis()->SetLimits(0.,high);
  grScorrvsC->GetXaxis()->SetTitle("E_{S}");
  grScorrvsC->GetYaxis()->SetRangeUser(0.,high);
  grScorrvsC->GetYaxis()->SetTitle("E_{C}");
  grScorrvsC->SetMaximum(high);
  grScorrvsC->SetMinimum(0.);
  grScorrvsC->Draw("ap");
  c->SaveAs(filename+"_ScorrvsC.png");

  tT_C->Draw("Hist"); c->SaveAs(filename+"_tC.png");
  tT_S->Draw("Hist"); c->SaveAs(filename+"_tS.png");
  tWav_C->Draw("Hist"); c->SaveAs(filename+"_wavC.png");
  tWav_S->Draw("Hist"); c->SaveAs(filename+"_wavS.png");
  tNhit_C->Draw("Hist"); c->SaveAs(filename+"_nhitC.png");
  tNhit_S->Draw("Hist"); c->SaveAs(filename+"_nhitS.png");

  /*std::string fileout_CSS = "/u/user/syjang/scratch/DRC_generic/results/ele/length";
  std::ofstream ofstream1;

  ofstream1.open(fileout_CSS+"/"+material+"_CSS_Error.csv", std::ios::out | std::ios::app);
  ofstream1 << en_val << "   C    mean    "<<tE_C->GetMean() << "    " << tE_C->GetMeanError() << std::endl;
  ofstream1 << en_val << "   C    sig     "<<tE_C->GetRMS() << "    " << tE_C->GetRMSError() << std::endl;
  ofstream1 << en_val << "   S        mean    "<<tE_S->GetMean() << "    " << tE_S->GetMeanError() << std::endl;
  ofstream1 << en_val << "   S        sig     "<<tE_S->GetRMS() << "    "  << tE_S->GetRMSError() << std::endl;
  //ofstream1 << en_val << "   S    mean    "<<tE_Scorr->GetMean() << "    " << tE_Scorr->GetMeanError() << std::endl;
  //ofstream1 << en_val << "   S    sig     "<<tE_Scorr->GetRMS() << "    "  << tE_Scorr->GetRMSError() << std::endl;
  ofstream1 << en_val << "   Sum      mean    "<<grE_SC->GetParameter(1) << "    " << grE_SC->GetParError(1) << std::endl;
  ofstream1 << en_val << "   Sum      sig     "<<grE_SC->GetParameter(2) << "    " << grE_SC->GetParError(2) << std::endl;
  //ofstream1 << en_val << "   Sum  mean    "<<grE_DR->GetParameter(1) << "    " << grE_DR->GetParError(1) << std::endl;
  //ofstream1 << en_val << "   Sum  sig     "<<grE_DR->GetParameter(2) << "    " << grE_DR->GetParError(2) << std::endl;*/

}
