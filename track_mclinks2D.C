
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
int track_mclinks2D(const char* input, const char* output, int pdg_number) {
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

    // Maps for the misreconstructed particles
    std::map<std::string, std::string> MisReco_Map_pt = {
    {"h_pass_pt", Form("%s correctly reconstructed (wgt > 0.75)",name.c_str())},
    {"h_photon_pt", "#gamma"},
    {"h_neutron_pt", "n"},
    {"h_electron_pt", "e^{#pm} (wgt < 0.75)"},
    {"h_pion_pt", "#pi^{#pm}"},
    {"h_other_pt", "Other particles"},
    {"h_notmatched_pt", Form("%s with no reconstructed particle matched",name.c_str())}
    };

    std::map<std::string, std::string> MisReco_Map_p = {
    {"h_pass_p", Form("%s correctly reconstructed (wgt > 0.75)",name.c_str())},
    {"h_photon_p", "#gamma"},
    {"h_neutron_p", "n"},
    {"h_electron_p", "e^{#pm} (wgt < 0.75)"},
    {"h_pion_p", "#pi^{#pm}"},
    {"h_other_p", "Other particles"},
    {"h_notmatched_p", Form("%s with no reconstructed particle matched",name.c_str())}
    };
    
    // Maps for the histograms
    std::map<std::string,TH2D*> Histo_Map_pt;
    Histo_Map_pt["h_pass_pt"] = new TH2D("h_pass_pt", ";#theta [mrad]", 30, 0, TMath::Pi()*1000, 30, 0, 100);  
    Histo_Map_pt["h_photon_pt"] = new TH2D("h_photon_pt", ";#theta [mrad]", 30, 0, TMath::Pi()*1000, 30, 0, 100);
    Histo_Map_pt["h_neutron_pt"] = new TH2D("h_neutron_pt", ";#theta [mrad]", 30, 0, TMath::Pi()*1000, 30, 0, 100);
    Histo_Map_pt["h_electron_pt"] = new TH2D("h_electron_pt", ";#theta [mrad]", 30, 0, TMath::Pi()*1000, 30, 0, 100);
    Histo_Map_pt["h_pion_pt"] = new TH2D("h_pion_pt", ";#theta [mrad]", 30, 0, TMath::Pi()*1000, 30, 0, 100);
    Histo_Map_pt["h_other_pt"] = new TH2D("h_other_pt",";#theta [mrad]", 30, 0, TMath::Pi()*1000, 30, 0, 100);
    Histo_Map_pt["h_notmatched_pt"] = new TH2D("h_notmatched_pt",";#theta [mrad]", 30, 0, TMath::Pi()*1000, 30, 0, 100);

    std::map<std::string,TH2D*> Histo_Map_p;
    Histo_Map_p["h_pass_p"] = new TH2D("h_pass_p", ";#theta [mrad]", 30, 0, TMath::Pi()*1000, 30, 0, 100);   
    Histo_Map_p["h_photon_p"] = new TH2D("h_photon_p",";#theta [mrad]", 30, 0, TMath::Pi()*1000, 30, 0, 100);
    Histo_Map_p["h_neutron_p"] = new TH2D("h_neutron_p", ";#theta [mrad]",30, 0, TMath::Pi()*1000, 30, 0, 100);
    Histo_Map_p["h_electron_p"] = new TH2D("h_electron_p",";#theta [mrad]", 30, 0, TMath::Pi()*1000, 30, 0, 100);
    Histo_Map_p["h_pion_p"] = new TH2D("h_pion_p", ";#theta [mrad]", 30, 0, TMath::Pi()*1000, 30, 0, 100);
    Histo_Map_p["h_other_p"] = new TH2D("h_other_p", ";#theta [mrad]", 30, 0, TMath::Pi()*1000, 30, 0, 100);
    Histo_Map_p["h_notmatched_p"] = new TH2D("h_notmatched_p", ";#theta [mrad]", 30, 0, TMath::Pi()*1000, 30, 0, 100);

    TH2D* h_total_p = new TH2D("h_total_p", ";#theta [mrad]",30, 0, TMath::Pi()*1000, 30, 0, 100);
    TH2D* h_total_pt = new TH2D("h_total_pt", ";#theta [mrad]",30, 0, TMath::Pi()*1000, 30, 0, 100);

    int all = 0;
  
    for (size_t i = 0; i < reader.getEntries("events"); ++i) {
        if (all >= 100000) break;
        all++;

        const auto event = podio::Frame(reader.readNextEntry("events"));
        const auto& mcparticles = event.get<edm4hep::MCParticleCollection>("MCParticles");
        const auto& mclinks = event.get<edm4hep::MCRecoParticleAssociationCollection>("MCTruthRecoLink");
    
        // Search for particles
        for (const auto& mcp: mcparticles) {
            MC_Particle particle;
            particle = getStableMCParticle(mcp, particle, pdg_number);
            if (particle.flag == true) {
                h_total_pt->Fill(particle.theta*1000,particle.pt);
                h_total_p->Fill(particle.theta*1000,particle.p);
                Reconstructed_Particle reco_p;
                reco_p = getReco(mcp, mclinks, reco_p); 
                if(reco_p.reco.first.getPDG() == mcp.getPDG() && reco_p.trwgt.first > 0.75) {
                    Histo_Map_pt["h_pass_pt"]->Fill(particle.theta*1000,particle.pt);
                    Histo_Map_p["h_pass_p"]->Fill(particle.theta*1000,particle.p);
                }
                else {                                     
                    if (std::abs(reco_p.reco.first.getPDG()) == 11) Histo_Map_pt["h_electron_pt"]->Fill(particle.theta*1000, particle.pt);
                    else if (std::abs(reco_p.reco.first.getPDG()) == 211) Histo_Map_pt["h_pion_pt"]->Fill(particle.theta*1000,particle.pt);
                    else if (std::abs(reco_p.reco.first.getPDG()) == 22) Histo_Map_pt["h_photon_pt"]->Fill(particle.theta*1000,particle.pt);
                    else if (std::abs(reco_p.reco.first.getPDG()) == 2112) Histo_Map_pt["h_neutron_pt"]->Fill(particle.theta*1000,particle.pt);
                    else if (std::abs(reco_p.reco.first.getPDG()) == 0) Histo_Map_pt["h_notmatched_pt"]->Fill(particle.theta*1000,particle.pt);
                    else Histo_Map_pt["h_other_pt"]->Fill(particle.theta*1000,particle.pt);

                    if (std::abs(reco_p.reco.first.getPDG()) == 11) Histo_Map_p["h_electron_p"]->Fill(particle.theta*1000,particle.p);
                    else if (std::abs(reco_p.reco.first.getPDG()) == 211) Histo_Map_p["h_pion_p"]->Fill(particle.theta*1000,particle.p);
                    else if (std::abs(reco_p.reco.first.getPDG()) == 22) Histo_Map_p["h_photon_p"]->Fill(particle.theta*1000,particle.p);
                    else if (std::abs(reco_p.reco.first.getPDG()) == 2112) Histo_Map_p["h_neutron_p"]->Fill(particle.theta*1000,particle.p);
                    else if (std::abs(reco_p.reco.first.getPDG()) == 0) Histo_Map_p["h_notmatched_p"]->Fill(particle.theta*1000,particle.p);
                    else Histo_Map_p["h_other_p"]->Fill(particle.theta*1000,particle.p);
                }               
            }
        }
    }

    // Percentage of reconstructed particles
    for (auto& histo: Histo_Map_pt) {   
        histo.second->Divide(h_total_pt);
    }
    for (auto& histo: Histo_Map_p) {   
        histo.second->Divide(h_total_p);
    }

    // Create graphs 
    for (auto& histo: Histo_Map_pt) {   
        histo.second->GetYaxis()->SetRangeUser(0, 100);
        histo.second->GetXaxis()->SetRangeUser(0, TMath::Pi()*1000); 
        histo.second->GetYaxis()->SetTitle("p_{T} [GeV]");
        histo.second->GetYaxis()->SetTitleSize(0.05);
        histo.second->GetXaxis()->SetTitleSize(0.05);
        histo.second->GetXaxis()->SetTitleOffset(0.8);
        histo.second->GetYaxis()->SetTitleOffset(0.8);
    }
    for (auto& histo: Histo_Map_p) {     
        histo.second->GetYaxis()->SetRangeUser(0, 100);
        histo.second->GetXaxis()->SetRangeUser(0, TMath::Pi()*1000); 
        histo.second->GetYaxis()->SetTitle("p [GeV]");
        histo.second->GetYaxis()->SetTitleSize(0.05);
        histo.second->GetXaxis()->SetTitleSize(0.05);
        histo.second->GetXaxis()->SetTitleOffset(0.8);
        histo.second->GetYaxis()->SetTitleOffset(0.8);
    }

    // Create canvas and draw histograms
    TLatex t;
    t.SetTextSize(30);
    t.SetTextFont(63);
    gStyle->Reset("Modern");
   
    for (auto& histo: Histo_Map_pt) { 
        auto misname = MisReco_Map_pt[histo.first.c_str()];
        TCanvas* c_other_pt = new TCanvas();
        histo.second->Draw("colz0");
        t.DrawLatexNDC(0.15, 0.93935, "CLD #font[52]{work in progress}");
        t.DrawLatexNDC(0.58, 0.93935, "#scale[0.8]{Highest track weight}");
        t.DrawLatexNDC(0.15, 0.02, Form("#scale[0.8]{%s}", misname.c_str()));
        c_other_pt->Draw();
        c_other_pt->SaveAs(Form("plots/trk_electrons/%s_thetapt_trwgt.pdf", misname.c_str()));
    }
    for (auto& histo: Histo_Map_p) { 
        auto misname = MisReco_Map_p[histo.first.c_str()];
        TCanvas* c_other_p = new TCanvas();
        histo.second->Draw("colz0");
        t.DrawLatexNDC(0.15, 0.93935, "CLD #font[52]{work in progress}");
        t.DrawLatexNDC(0.58, 0.93935, "#scale[0.8]{Highest track weight}");
        t.DrawLatexNDC(0.15, 0.02, Form("#scale[0.8]{%s}", misname.c_str()));
        c_other_p->Draw();
        c_other_p->SaveAs(Form("plots/trk_electrons/%s_thetap_trwgt.pdf", misname.c_str()));
    } 

    outputfile->Write();
    outputfile->Close(); 
    return 0;
}

