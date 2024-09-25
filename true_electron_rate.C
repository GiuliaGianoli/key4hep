
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <utility>
#include <map>
#include <cmath>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TVector3.h>
#include <TEfficiency.h>
#include <marlinutil/HelixClass_double.h>
#include "podio/ROOTReader.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/MCRecoTrackParticleAssociation.h"
#include "edm4hep/TrackCollection.h"
#include "edm4hep/TrackState.h"
#include "edm4hep/MCRecoParticleAssociationCollection.h"
#include "edm4hep/ReconstructedParticleCollection.h"
#include "edm4hep/utils/kinematics.h"
#include "podio/Frame.h"
#include "podio/ROOTFrameReader.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include <string>
#include <vector>
#include "TStyle.h"
#include "TGraph.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TLegend.h"
#include "getParticle.h"

// Main function
int true_electron_rate(const char* input, const char* output, int pdg_number) {
    gStyle->SetOptStat(0000);
    double Bz = 2;
    std::cout << Bz << std::endl;

    TFile* outputfile = new TFile(output, "RECREATE");

    // To read the generated file
    podio::ROOTReader reader;
    reader.openFile(input);

    // Map for the particles
    std::map<int, std::string> Particle_Map = {
    {11, "e^{#pm}"},
    {13, "#mu^{#pm}"},
    {211, "#pi^{#pm}"}
    };
    auto name = Particle_Map[pdg_number];

    // Maps for the histograms
    std::map<std::string,TH1D*> Histo_Map_theta; 
    Histo_Map_theta["h_pass_theta"] = new TH1D("h_pass_theta", ";#theta [mrad]", 30, 0, TMath::Pi()*1000);
    Histo_Map_theta["h_photon_theta"] = new TH1D("h_photon_theta", ";#theta [mrad]", 30, 0, TMath::Pi()*1000);
    Histo_Map_theta["h_neutron_theta"] = new TH1D("h_neutron_theta", ";#theta [mrad]", 30, 0, TMath::Pi()*1000);
    Histo_Map_theta["h_electron_theta"] = new TH1D("h_electron_theta", ";#theta [mrad]", 30, 0, TMath::Pi()*1000); 
    Histo_Map_theta["h_pion_theta"] = new TH1D("h_pion_theta", ";#theta [mrad]", 30, 0, TMath::Pi()*1000);
    Histo_Map_theta["h_other_theta"] = new TH1D("h_other_theta", ";#theta [mrad]", 30, 0, TMath::Pi()*1000);
    Histo_Map_theta["h_notmatched_theta"] = new TH1D("h_notmatched_theta", ";#theta [mrad]", 30, 0, TMath::Pi()*1000);

    std::map<std::string,TH1D*> Histo_Map_pt;
    Histo_Map_pt["h_pass_pt"] = new TH1D("h_pass_pt", ";p_{T} [GeV]", 30, 0, 100);  
    Histo_Map_pt["h_photon_pt"] = new TH1D("h_photon_pt", ";p_{T} [GeV]", 30, 0, 100);
    Histo_Map_pt["h_neutron_pt"] = new TH1D("h_neutron_pt", ";p_{T} [GeV]", 30, 0, 100);
    Histo_Map_pt["h_electron_pt"] = new TH1D("h_electron_pt", ";p_{T} [GeV]", 30, 0, 100);
    Histo_Map_pt["h_pion_pt"] = new TH1D("h_pion_pt", ";p_{T} [GeV]", 30, 0, 100);
    Histo_Map_pt["h_other_pt"] = new TH1D("h_other_pt", ";p_{T} [GeV]", 30, 0, 100);
    Histo_Map_pt["h_notmatched_pt"] = new TH1D("h_notmatched_pt", ";p_{T} [GeV]", 30, 0, 100);

    std::map<std::string,TH1D*> Histo_Map_p;
    Histo_Map_p["h_pass_p"] = new TH1D("h_pass_p", ";p [GeV]", 30, 0, 100);   
    Histo_Map_p["h_photon_p"] = new TH1D("h_photon_p", ";p [GeV]", 30, 0, 100);
    Histo_Map_p["h_neutron_p"] = new TH1D("h_neutron_p", ";p [GeV]", 30, 0, 100);
    Histo_Map_p["h_electron_p"] = new TH1D("h_electron_p", ";p [GeV]", 30, 0, 100);
    Histo_Map_p["h_pion_p"] = new TH1D("h_pion_p", ";p [GeV]", 30, 0, 100);
    Histo_Map_p["h_other_p"] = new TH1D("h_other_p", ";p [GeV]", 30, 0, 100);
    Histo_Map_p["h_notmatched_p"] = new TH1D("h_notmatched_p", ";p_{T} [GeV]", 30, 0, 100);

    std::map<std::string,TH1D*> Histo_Map_void;
    Histo_Map_void["void_theta"] = new TH1D("void_theta", ";#theta [mrad]", 30, 0, TMath::Pi()*1000);
    Histo_Map_void["void_pt"] = new TH1D("void_pt", ";p_{T} [GeV]", 30, 0, TMath::Pi()*1000);
    Histo_Map_void["void_p"] = new TH1D("void_p", ";p [GeV]", 30, 0, TMath::Pi()*1000);

    TH1D* h_total_p = new TH1D("h_total_p", ";p [GeV]", 30, 0, 100);
    TH1D* h_total_theta = new TH1D("h_total_theta", ";#theta [mrad]", 30, 0, TMath::Pi()*1000);
    TH1D* h_total_pt = new TH1D("h_total_pt", ";p_{T} [GeV]", 30, 0, 100);

    int all = 0;
  
    for (size_t i = 0; i < reader.getEntries("events"); ++i) {
        if (all >= 10000) break;
        all++;

        const auto event = podio::Frame(reader.readNextEntry("events"));
        const auto& mcparticles = event.get<edm4hep::MCParticleCollection>("MCParticles");
        const auto& mclinks = event.get<edm4hep::MCRecoParticleAssociationCollection>("MCTruthRecoLink");
    
        // Search for particles
        for (const auto& mcp: mcparticles) {
            MC_Particle particle;
            particle = getStableMCParticle(mcp, particle, pdg_number);
            if (particle.flag == true) {
                h_total_theta->Fill(particle.theta*1000);
                h_total_pt->Fill(particle.pt);
                h_total_p->Fill(particle.p);
                Reconstructed_Particle reco_p;
                reco_p = getReco(mcp, mclinks, reco_p); 
                if(reco_p.reco.first.getPDG() == mcp.getPDG() && reco_p.trwgt.first > 0.75) {
                    Histo_Map_theta["h_pass_theta"]->Fill(particle.theta *1000);
                    Histo_Map_pt["h_pass_pt"]->Fill(particle.pt);
                    Histo_Map_p["h_pass_p"]->Fill(particle.p);
                }
                else {
                    if (std::abs(reco_p.reco.first.getPDG()) == 11) Histo_Map_theta["h_electron_theta"]->Fill(particle.theta*1000);
                    else if (std::abs(reco_p.reco.first.getPDG()) == 211) Histo_Map_theta["h_pion_theta"]->Fill(particle.theta*1000);
                    else if (std::abs(reco_p.reco.first.getPDG()) == 22) Histo_Map_theta["h_photon_theta"]->Fill(particle.theta*1000);
                    else if (std::abs(reco_p.reco.first.getPDG()) == 2112) Histo_Map_theta["h_neutron_theta"]->Fill(particle.theta*1000);
                    else if (std::abs(reco_p.reco.first.getPDG()) == 0) Histo_Map_theta["h_notmatched_theta"]->Fill(particle.theta*1000);
                    else Histo_Map_theta["h_other_theta"]->Fill(particle.theta*1000);
                        
                                       
                    if (std::abs(reco_p.reco.first.getPDG()) == 11) Histo_Map_pt["h_electron_pt"]->Fill(particle.pt);
                    else if (std::abs(reco_p.reco.first.getPDG()) == 211) Histo_Map_pt["h_pion_pt"]->Fill(particle.pt);
                    else if (std::abs(reco_p.reco.first.getPDG()) == 22) Histo_Map_pt["h_photon_pt"]->Fill(particle.pt);
                    else if (std::abs(reco_p.reco.first.getPDG()) == 2112) Histo_Map_pt["h_neutron_pt"]->Fill(particle.pt);
                    else if (std::abs(reco_p.reco.first.getPDG()) == 0) Histo_Map_pt["h_notmatched_pt"]->Fill(particle.pt);
                    else Histo_Map_pt["h_other_pt"]->Fill(particle.pt);

                    if (std::abs(reco_p.reco.first.getPDG()) == 11) Histo_Map_p["h_electron_p"]->Fill(particle.p);
                    else if (std::abs(reco_p.reco.first.getPDG()) == 211) Histo_Map_p["h_pion_p"]->Fill(particle.p);
                    else if (std::abs(reco_p.reco.first.getPDG()) == 22) Histo_Map_p["h_photon_p"]->Fill(particle.p);
                    else if (std::abs(reco_p.reco.first.getPDG()) == 2112) Histo_Map_p["h_neutron_p"]->Fill(particle.p);
                    else if (std::abs(reco_p.reco.first.getPDG()) == 0) Histo_Map_p["h_notmatched_p"]->Fill(particle.p);
                    else Histo_Map_p["h_other_p"]->Fill(particle.p);
                }               
            }
        }
    }

    // Map for the  TEfficency
    std::map<std::string,TEfficiency*> Histo_Map_eff; 
    Histo_Map_eff["eff_theta"] = new TEfficiency(*Histo_Map_theta["h_pass_theta"], *h_total_theta);
    Histo_Map_eff["eff_theta"]->CreateGraph();
    Histo_Map_eff["eff_pt"] = new TEfficiency(*Histo_Map_pt["h_pass_pt"], *h_total_pt);
    Histo_Map_eff["eff_pt"]->CreateGraph();
    Histo_Map_eff["eff_p"] = new TEfficiency(*Histo_Map_p["h_pass_p"], *h_total_p);
    Histo_Map_eff["eff_p"]->CreateGraph(); 

    // Percentage of reconstructed particles
    for (auto& histo: Histo_Map_theta) {   
        histo.second->Divide(h_total_theta);
    }
    for (auto& histo: Histo_Map_pt) {   
        histo.second->Divide(h_total_pt);
    }
    for (auto& histo: Histo_Map_p) {   
        histo.second->Divide(h_total_p);
    }

    // Create graphs 
    for (auto& histo: Histo_Map_eff) {   
        histo.second->SetMarkerSize(2);
        histo.second->SetLineWidth(2);
        histo.second->SetLineColor(kRed);
        histo.second->SetMarkerColor(kRed);
    }
    for (auto& histo: Histo_Map_void) {
        histo.second->GetYaxis()->SetRangeUser(0.0001, 1.); 
        histo.second->GetYaxis()->SetTitle(Form("True %s rate", name.c_str()));
        histo.second->GetYaxis()->SetTitleSize(0.05);
        histo.second->GetXaxis()->SetTitleSize(0.05);
        histo.second->GetXaxis()->SetTitleOffset(0.8);
        histo.second->GetYaxis()->SetTitleOffset(0.8);
    }
    Histo_Map_void["void_theta"]->GetXaxis()->SetRangeUser(0, TMath::Pi()*1000);
    Histo_Map_void["void_pt"]->GetXaxis()->SetRangeUser(0, 100); 
    Histo_Map_void["void_p"]->GetXaxis()->SetRangeUser(0, 100); 
  
    for (auto& histo: Histo_Map_theta) {   
        histo.second->SetMarkerSize(2);
        histo.second->SetLineWidth(2);
    }
    Histo_Map_theta["h_photon_theta"]->GetYaxis()->SetRangeUser(0.0001, 1.05);
    Histo_Map_theta["h_photon_theta"]->GetXaxis()->SetRangeUser(0, TMath::Pi()*1000); 
    Histo_Map_theta["h_photon_theta"]->GetYaxis()->SetTitle("Reconstructed particles (%)");
    Histo_Map_theta["h_photon_theta"]->GetYaxis()->SetTitleSize(0.05);
    Histo_Map_theta["h_photon_theta"]->GetXaxis()->SetTitleSize(0.05);
    Histo_Map_theta["h_photon_theta"]->GetXaxis()->SetTitleOffset(0.8);
    Histo_Map_theta["h_photon_theta"]->GetYaxis()->SetTitleOffset(0.8);
    Histo_Map_theta["h_photon_theta"]->SetLineColor(kYellow+1);
    Histo_Map_theta["h_photon_theta"]->SetMarkerColor(kYellow+1);
    Histo_Map_theta["h_neutron_theta"]->SetLineColor(kCyan+1);
    Histo_Map_theta["h_neutron_theta"]->SetMarkerColor(kCyan+1);
    Histo_Map_theta["h_electron_theta"]->SetLineColor(kRed-7);
    Histo_Map_theta["h_electron_theta"]->SetMarkerColor(kRed-7);
    Histo_Map_theta["h_pion_theta"]->SetLineColor(kBlue+1);
    Histo_Map_theta["h_pion_theta"]->SetMarkerColor(kBlue+1);
    Histo_Map_theta["h_other_theta"]->SetLineColor(kRed+2);
    Histo_Map_theta["h_other_theta"]->SetMarkerColor(kRed+2);
    Histo_Map_theta["h_pass_theta"]->SetLineColor(kBlack);
    Histo_Map_theta["h_pass_theta"]->SetMarkerColor(kBlack);
    Histo_Map_theta["h_notmatched_theta"]->SetLineColor(kGreen+2);
    Histo_Map_theta["h_notmatched_theta"]->SetMarkerColor(kGreen+2);

    for (auto& histo: Histo_Map_pt) {   
        histo.second->SetMarkerSize(2);
        histo.second->SetLineWidth(2);
    }
    Histo_Map_pt["h_photon_pt"]->GetYaxis()->SetRangeUser(0.0001, 1.05);
    Histo_Map_pt["h_photon_pt"]->GetXaxis()->SetRangeUser(0, 100); 
    Histo_Map_pt["h_photon_pt"]->GetYaxis()->SetTitle("Reconstructed particles (%)");
    Histo_Map_pt["h_photon_pt"]->GetYaxis()->SetTitleSize(0.05);
    Histo_Map_pt["h_photon_pt"]->GetXaxis()->SetTitleSize(0.05);
    Histo_Map_pt["h_photon_pt"]->GetXaxis()->SetTitleOffset(0.8);
    Histo_Map_pt["h_photon_pt"]->GetYaxis()->SetTitleOffset(0.8);
    Histo_Map_pt["h_photon_pt"]->SetLineColor(kYellow+1);
    Histo_Map_pt["h_photon_pt"]->SetMarkerColor(kYellow+1);
    Histo_Map_pt["h_neutron_pt"]->SetLineColor(kCyan+1);
    Histo_Map_pt["h_neutron_pt"]->SetMarkerColor(kCyan+1);
    Histo_Map_pt["h_electron_pt"]->SetLineColor(kRed-7);
    Histo_Map_pt["h_electron_pt"]->SetMarkerColor(kRed-7);
    Histo_Map_pt["h_pion_pt"]->SetLineColor(kBlue+1);
    Histo_Map_pt["h_pion_pt"]->SetMarkerColor(kBlue+1);
    Histo_Map_pt["h_other_pt"]->SetLineColor(kRed+2);
    Histo_Map_pt["h_other_pt"]->SetMarkerColor(kRed+2);
    Histo_Map_pt["h_pass_pt"]->SetLineColor(kBlack);
    Histo_Map_pt["h_pass_pt"]->SetMarkerColor(kBlack);
    Histo_Map_pt["h_notmatched_pt"]->SetLineColor(kGreen+2);
    Histo_Map_pt["h_notmatched_pt"]->SetMarkerColor(kGreen+2);

    for (auto& histo: Histo_Map_p) {   
        histo.second->SetMarkerSize(2);
        histo.second->SetLineWidth(2);
    }
    Histo_Map_p["h_photon_p"]->GetYaxis()->SetRangeUser(0.0001, 1.05);
    Histo_Map_p["h_photon_p"]->GetXaxis()->SetRangeUser(0, 100); 
    Histo_Map_p["h_photon_p"]->GetYaxis()->SetTitle("Reconstructed particles (%)");
    Histo_Map_p["h_photon_p"]->GetYaxis()->SetTitleSize(0.05);
    Histo_Map_p["h_photon_p"]->GetXaxis()->SetTitleSize(0.05);
    Histo_Map_p["h_photon_p"]->GetXaxis()->SetTitleOffset(0.8);
    Histo_Map_p["h_photon_p"]->GetYaxis()->SetTitleOffset(0.8);
    Histo_Map_p["h_photon_p"]->SetLineColor(kYellow+1);
    Histo_Map_p["h_photon_p"]->SetMarkerColor(kYellow+1);
    Histo_Map_p["h_neutron_p"]->SetLineColor(kCyan+1);
    Histo_Map_p["h_neutron_p"]->SetMarkerColor(kCyan+1);
    Histo_Map_p["h_electron_p"]->SetLineColor(kRed-7);
    Histo_Map_p["h_electron_p"]->SetMarkerColor(kRed-7);
    Histo_Map_p["h_pion_p"]->SetLineColor(kBlue+1);
    Histo_Map_p["h_pion_p"]->SetMarkerColor(kBlue+1);
    Histo_Map_p["h_other_p"]->SetLineColor(kRed+2);
    Histo_Map_p["h_other_p"]->SetMarkerColor(kRed+2);
    Histo_Map_p["h_pass_p"]->SetLineColor(kBlack);
    Histo_Map_p["h_pass_p"]->SetMarkerColor(kBlack);
    Histo_Map_p["h_notmatched_p"]->SetLineColor(kGreen+2);
    Histo_Map_p["h_notmatched_p"]->SetMarkerColor(kGreen+2);

    // Create canvas and draw histograms
    TLatex t;
    t.SetTextSize(30);
    t.SetTextFont(63);
    gStyle->Reset("Modern");
    
    TCanvas* c_track_theta = new TCanvas();    
    Histo_Map_void["void_theta"]->Draw("p");
    Histo_Map_eff["eff_theta"]->Draw("psame");

    t.DrawLatexNDC(0.15, 0.93935, "CLD #font[52]{work in progress}");
    t.DrawLatexNDC(0.3, 0.7, "#scale[0.8]{Highest track weight}");
    //gPad->SetLogy();
    c_track_theta->Draw();
    c_track_theta->SaveAs("plots/True_electron_rate/ElectronRate_theta_trwgt.pdf"); 

    TCanvas* c_track_pt = new TCanvas();   
    Histo_Map_void["void_pt"]->Draw("p");
    Histo_Map_eff["eff_pt"]->Draw("psame");

    t.DrawLatexNDC(0.15, 0.93935, "CLD #font[52]{work in progress}");
    t.DrawLatexNDC(0.3, 0.7, "#scale[0.8]{Highest track weight}");
    //gPad->SetLogy();
    c_track_pt->Draw();
    c_track_pt->SaveAs("plots/True_electron_rate/ElectronRate_pt_trwgt.pdf");

    TCanvas* c_track_p = new TCanvas();   
    Histo_Map_void["void_p"]->Draw("p");
    Histo_Map_eff["eff_p"]->Draw("psame");


    t.DrawLatexNDC(0.15, 0.93935, "CLD #font[52]{work in progress}");
    t.DrawLatexNDC(0.3, 0.7, "#scale[0.8]{Highest track weight}");
    //gPad->SetLogy();
    c_track_p->Draw();
    c_track_p->SaveAs("plots/True_electron_rate/ElectronRate_p_trwgt.pdf");

    TCanvas* c_other_theta = new TCanvas();
    Histo_Map_theta["h_photon_theta"]->Draw("hist");
    for (auto& histo: Histo_Map_theta) {   
        histo.second->Draw("histsame");
    }

    t.DrawLatexNDC(0.15, 0.93935, "CLD #font[52]{work in progress}");
    t.DrawLatexNDC(0.3, 0.7, "#scale[0.8]{Highest track weight}");
    TLegend* l = new TLegend(0.65, 0.65, 1, 1);
    l->AddEntry(Histo_Map_theta["h_photon_theta"], "#gamma");
    l->AddEntry(Histo_Map_theta["h_neutron_theta"], "n");
    l->AddEntry(Histo_Map_theta["h_electron_theta"], "e^{#pm} (wgt < 0.75)");
    l->AddEntry(Histo_Map_theta["h_pion_theta"], "#pi^{#pm}");
    l->AddEntry(Histo_Map_theta["h_other_theta"], "Other particles");
    l->AddEntry(Histo_Map_theta["h_pass_theta"], Form("%s correctly reconstructed (wgt > 0.75)", name.c_str()));
    l->AddEntry(Histo_Map_theta["h_notmatched_theta"], Form("%s with no reconstructed particle matched", name.c_str()));
    l->Draw();
    //gPad->SetLogy();
    c_other_theta->Draw();
    c_other_theta->SaveAs("plots/True_electron_rate/NotReconstructed_theta_trwgt.pdf");

    TCanvas* c_other_pt = new TCanvas();
    Histo_Map_pt["h_photon_pt"]->Draw("hist");
    for (auto& histo: Histo_Map_pt) {   
        histo.second->Draw("histsame");
    }

    t.DrawLatexNDC(0.15, 0.93935, "CLD #font[52]{work in progress}");
    t.DrawLatexNDC(0.3, 0.7, "#scale[0.8]{Highest track weight}");
    TLegend* l2 = new TLegend(0.65, 0.65, 1, 1);
    l2->AddEntry(Histo_Map_pt["h_photon_pt"], "#gamma");
    l2->AddEntry(Histo_Map_pt["h_neutron_pt"], "n");
    l2->AddEntry(Histo_Map_pt["h_electron_pt"], "e^{#pm} (wgt < 0.75)");
    l2->AddEntry(Histo_Map_pt["h_pion_pt"], "#pi^{#pm}");
    l2->AddEntry(Histo_Map_pt["h_other_pt"], "Other particles");
    l2->AddEntry(Histo_Map_pt["h_pass_theta"], Form("%s correctly reconstructed (wgt > 0.75)", name.c_str()));
    l2->AddEntry(Histo_Map_pt["h_notmatched_pt"], Form("%s with no reconstructed particle matched", name.c_str()));
    l2->Draw(); 
    //gPad->SetLogy();
    c_other_pt->Draw();
    c_other_pt->SaveAs("plots/True_electron_rate/NotReconstructed_pt_trwgt.pdf");

    TCanvas* c_other_p = new TCanvas();
    Histo_Map_p["h_photon_p"]->Draw("hist");
       for (auto& histo: Histo_Map_p) {   
        histo.second->Draw("histsame");
    }

    t.DrawLatexNDC(0.15, 0.93935, "CLD #font[52]{work in progress}");
    t.DrawLatexNDC(0.3, 0.7, "#scale[0.8]{Highest track weight}");
    TLegend* l3 = new TLegend(0.65, 0.65, 1, 1);
    l3->AddEntry(Histo_Map_p["h_photon_p"], "#gamma");
    l3->AddEntry(Histo_Map_p["h_neutron_p"], "n");
    l3->AddEntry(Histo_Map_p["h_electron_p"], "e^{#pm} (wgt < 0.75)");
    l3->AddEntry(Histo_Map_p["h_pion_p"], "#pi^{#pm}");
    l3->AddEntry(Histo_Map_p["h_other_p"], "Other particles");
    l3->AddEntry(Histo_Map_p["h_pass_theta"], Form("%s correctly reconstructed (wgt > 0.75)", name.c_str()));
    l3->AddEntry(Histo_Map_p["h_notmatched_p"], Form("%s with no reconstructed particle matched", name.c_str()));
    l3->Draw(); 
    //gPad->SetLogy();
    c_other_p->Draw();
    c_other_p->SaveAs("plots/True_electron_rate/NotReconstructed_p_trwgt.pdf");

    outputfile->Write();
    outputfile->Close(); 
    return 0;
}

