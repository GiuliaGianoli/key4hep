
#include <iostream>
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

#include "edm4hep/ReconstructedParticleCollection.h"
#include "edm4hep/utils/kinematics.h"
#include "podio/Frame.h"
#include "podio/ROOTFrameReader.h"
#include "TFile.h"
#include "TH1D.h"
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

void fakerate_plots() {
    gStyle->SetOptStat(0000);
    double Bz = 2;
    std::cout << Bz << std::endl;

    TFile* outputfile = new TFile("Outputfile.root", "RECREATE");

    // To read the generated file
    podio::ROOTReader reader;
    reader.openFile("1Step1_wogun_output__edm4hep.root");
    
    // Create histogram

    TH1D* h_pass_theta = new TH1D("h_pass_theta", ";#theta [mrad]", 30, 0, 250);
    TH1D* h_total_theta = new TH1D("h_total_theta", ";#theta [mrad]", 30, 0, 250);
    TH1D* h_pass_pt = new TH1D("h_pass_pt", ";p_{T} [GeV]", 30, 0, 100);
    TH1D* h_total_pt = new TH1D("h_total_pt", ";p_{T} [GeV]", 30, 0, 100);

    int all = 0;

    for (size_t i = 0; i < reader.getEntries("events"); ++i) {
        if (all >= 100000) break;
        all++;

        const auto event = podio::Frame(reader.readNextEntry("events"));
        const auto &mcparticles = event.get<edm4hep::MCParticleCollection>("MCParticles");
        const auto &mctruthrecolinks = event.get<edm4hep::MCRecoTrackParticleAssociationCollection>("SiTracksMCTruthLink");
        const auto &tracks = event.get<edm4hep::TrackCollection>("SiTracks_Refitted");
        
        // Search for tracks 
        for (const auto& track: tracks) {
            const edm4hep::TrackState& track_state = track.getTrackStates(0);
            auto h = HelixClass_double(); 
            h.Initialize_Canonical(track_state.phi, track_state.D0, track_state.Z0, track_state.omega, track_state.tanLambda, 2);
            auto mom = h.getMomentum(); 
            auto pt = h.getPXY();
            auto theta = std::atan2(pt,mom[2]);
            h_total_theta->Fill(theta *1000); 
            h_total_pt->Fill(pt);
            // Fake electron           
            for (const auto& link : mctruthrecolinks) {
                if (link.getRec() == track && link.getWeight() < 0.75) {    
                    h_pass_theta->Fill(theta *1000);
                    h_pass_pt->Fill(pt);
                }                                                        
            } 
        }
    }

    //Fake electron rate
    TH1D* u_theta = (TH1D*)h_pass_theta->Clone(); 
    u_theta->Scale(1 / h_total_theta->Integral());
    u_theta->SetLineColor(kRed);
    u_theta->SetMarkerColor(kRed);
    u_theta->SetMarkerSize(2);
    u_theta->SetLineWidth(2);
    u_theta->GetYaxis()->SetRangeUser(0.0001, 1.);
    u_theta->GetXaxis()->SetRangeUser(0, 250);
    u_theta->GetYaxis()->SetTitle("Fake rate"); 
    u_theta->GetYaxis()->SetTitleSize(0.05);
    u_theta->GetXaxis()->SetTitleSize(0.05);
    u_theta->GetXaxis()->SetTitleOffset(0.8);
    u_theta->GetYaxis()->SetTitleOffset(0.8);
 

    TH1D* u_pt = (TH1D*)h_pass_pt->Clone(); 
    u_pt->Scale(1 / h_total_pt->Integral());
    u_pt->SetLineColor(kRed);
    u_pt->SetMarkerColor(kRed);
    u_pt->SetMarkerSize(2);
    u_pt->SetLineWidth(2);
    u_pt->GetYaxis()->SetRangeUser(0.0001, 1.);
    u_pt->GetXaxis()->SetRangeUser(0, 100); 
    u_pt->GetYaxis()->SetTitle("Fake rate");
    u_pt->GetYaxis()->SetTitleSize(0.05);
    u_pt->GetXaxis()->SetTitleSize(0.05);
    u_pt->GetXaxis()->SetTitleOffset(0.8);
    u_pt->GetYaxis()->SetTitleOffset(0.8);

    // Create canvas and draw histogram 
    TLatex t;
    t.SetTextSize(30);
    t.SetTextFont(63);
    gStyle->Reset("Modern");
    
    TCanvas* c_track_theta = new TCanvas();    
    u_theta->Draw("p");
    t.DrawLatexNDC(0.15, 0.93935, "CLD #font[52]{work in progress}");
    gPad->SetLogy();
    c_track_theta->Draw();
    c_track_theta->SaveAs("plots/trk_hits_wogun_fake/FakeTrackingRate_theta.pdf"); 

    TCanvas* c_track_pt = new TCanvas();   
    u_pt->Draw("p");
    t.DrawLatexNDC(0.15, 0.93935, "CLD #font[52]{work in progress}");
    gPad->SetLogy();
    c_track_pt->Draw();
    c_track_pt->SaveAs("plots/trk_hits_wogun_fake/FakeTrackingRate_pt.pdf");

    outputfile->Write();
    outputfile->Close(); 
}

