
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
using namespace std;

void track_hits() {
    gStyle->SetOptStat(0000);
    double Bz = 2;
    std::cout << Bz << std::endl;

    TFile* outputfile = new TFile("Outputfile.root", "RECREATE");

    // To read the generated file
    podio::ROOTReader reader;
    reader.openFile("1Step1_wogun_output__edm4hep.root");
    
    // Create histogram
    TH1D* h_events = new TH1D("h_events", ";#theta [mrad]", 30, 0, 250);
    TH1D* h_itec_hits = new TH1D("h_itec_hits", ";#theta [mrad]", 30, 0, 250);
    TH1D* h_itbc_hits = new TH1D("h_itbc_hits", ";#theta [mrad]", 30, 0, 250);
    TH1D* h_otbc_hits = new TH1D("h_otbc_hits", ";#theta [mrad]", 30, 0, 250);
    TH1D* h_otec_hits = new TH1D("h_otec_hits", ";#theta [mrad]", 30, 0, 250);
    TH1D* h_vbc_hits = new TH1D("h_vbc_hits", ";#theta [mrad]", 30, 0, 250);
    TH1D* h_vec_hits = new TH1D("h_vec_hits", ";#theta [mrad]", 30, 0, 250);
    TH1D* h_total_hits = new TH1D("h_total_hits", ";#theta [mrad]", 30, 0, 250);

    TH1D* h_pass_theta = new TH1D("h_pass_theta", ";#theta [mrad]", 30, 0, 250);
    TH1D* h_total_theta = new TH1D("h_total_theta", ";#theta [mrad]", 30, 0, 250);

    int all = 0;
    int found = 0;
    int all_acc = 0;
    int all_c = 0;
    int prova = 0;
    

    for (size_t i = 0; i < reader.getEntries("events"); ++i) {
        if (all >= 100000) break;
        all++;

        const auto event = podio::Frame(reader.readNextEntry("events"));
        const auto &mcparticles = event.get<edm4hep::MCParticleCollection>("MCParticles");
        const auto &mctruthrecolinks = event.get<edm4hep::MCRecoTrackParticleAssociationCollection>("SiTracksMCTruthLink");

        const auto &itec = event.get<edm4hep::SimTrackerHitCollection>("InnerTrackerEndcapCollection");
        const auto &itbc = event.get<edm4hep::SimTrackerHitCollection>("InnerTrackerBarrelCollection");
        const auto &otbc = event.get<edm4hep::SimTrackerHitCollection>("OuterTrackerBarrelCollection");
        const auto &otec = event.get<edm4hep::SimTrackerHitCollection>("OuterTrackerEndcapCollection");
        const auto &vbc = event.get<edm4hep::SimTrackerHitCollection>("VertexBarrelCollection");
        const auto &vec = event.get<edm4hep::SimTrackerHitCollection>("VertexEndcapCollection");
        
        
        // Search for first genStat 1 electron
        for (const auto& mcp : mcparticles) {
            if (abs(mcp.getPDG()) == 11 && mcp.getGeneratorStatus() == 1) {  
                auto mom = mcp.getMomentum(); 
                TVector3 mom_v(mom.x, mom.y, mom.z);  
                h_events->Fill(mom_v.Theta() *1000); 
                h_total_theta->Fill(mom_v.Theta() *1000);                    
                for (const auto& hit : itec) {
                    if (hit.getParticle() == mcp) {                        
                        h_itec_hits->Fill(mom_v.Theta() *1000);
                        h_total_hits->Fill(mom_v.Theta() *1000);                                        
                    }
                }
                for (const auto& hit : itbc) {
                    if (hit.getParticle() == mcp) {                        
                        h_itbc_hits->Fill(mom_v.Theta() *1000); 
                        h_total_hits->Fill(mom_v.Theta() *1000);                      
                    }
                }
                for (const auto& hit : otbc) {
                    if (hit.getParticle() == mcp) {                        
                        h_otbc_hits->Fill(mom_v.Theta() *1000); 
                        h_total_hits->Fill(mom_v.Theta() *1000);                    
                    }
                }
                for (const auto& hit : otec) {
                    if (hit.getParticle() == mcp) {                        
                        h_otec_hits->Fill(mom_v.Theta() *1000); 
                        h_total_hits->Fill(mom_v.Theta() *1000);                     
                    }
                }
                for (const auto& hit : vbc) {
                    if (hit.getParticle() == mcp) {                        
                        h_vbc_hits->Fill(mom_v.Theta() *1000); 
                        h_total_hits->Fill(mom_v.Theta() *1000);                      
                    }
                }
                for (const auto& hit : vec) {
                    if (hit.getParticle() == mcp) {                        
                        h_vec_hits->Fill(mom_v.Theta() *1000);
                        h_total_hits->Fill(mom_v.Theta() *1000);                      
                    }
                }
               for (const auto& link : mctruthrecolinks) {
                    if (link.getSim() == mcp && link.getWeight() > 0.99) {
                        h_pass_theta->Fill(mom_v.Theta() *1000);
                    }
                }
            }
        }            
    }

    // TEfficency
    TEfficiency* eff_theta = new TEfficiency(*h_pass_theta, *h_total_theta);
    eff_theta->CreateGraph();

    // Average amount of hits
    TH1D* u_itec_hits = (TH1D*)h_itec_hits->Clone(); 
    u_itec_hits->Divide(h_events);
    TH1D* u_itbc_hits = (TH1D*)h_itbc_hits->Clone(); 
    u_itbc_hits->Divide(h_events);
    TH1D* u_otbc_hits = (TH1D*)h_otbc_hits->Clone(); 
    u_otbc_hits->Divide(h_events);
    TH1D* u_otec_hits = (TH1D*)h_otec_hits->Clone(); 
    u_otec_hits->Divide(h_events);
    TH1D* u_vbc_hits = (TH1D*)h_vbc_hits->Clone(); 
    u_vbc_hits->Divide(h_events);
    TH1D* u_vec_hits = (TH1D*)h_vec_hits->Clone(); 
    u_vec_hits->Divide(h_events);
    TH1D* u_total_hits = (TH1D*)h_total_hits->Clone(); 
    u_total_hits->Divide(h_events); 

    // Create canvas and draw histograms
    TLatex t;
    t.SetTextSize(30);
    t.SetTextFont(63);
    gStyle->Reset("Modern");
    TCanvas* c_hits = new TCanvas("c", "comp", 200, 10, 700, 500);
    
    t.DrawLatex(50, 21, "Single e^{-}");
    t.DrawLatexNDC(0.15, 0.93935, "CLD #font[52]{work in progress}");

    TPad *p1 = new TPad("p1", "", 0, 0, 1, 1);
    TPad *p2 = new TPad("p2", "", 0, 0, 1, 1);
    p1->SetFillStyle(4000); // will be transparent
    p2->SetFillStyle(4000); 

    p1->Draw();
    p1->cd(); 
     
    // Create the first graph
    u_itec_hits->SetLineColor(kCyan+1);
    u_itec_hits->SetLineWidth(2);
    u_itec_hits->GetYaxis()->SetRangeUser(0,25);
    u_itec_hits->GetXaxis()->SetRangeUser(0, 250);
    u_itec_hits->GetYaxis()->SetTitle("Number of hits per electron");
    u_itec_hits->GetYaxis()->SetTitleSize(0.05);
    u_itec_hits->GetXaxis()->SetTitleSize(0.05);
    u_itec_hits->GetXaxis()->SetTitleOffset(0.8);
    u_itec_hits->GetYaxis()->SetTitleOffset(0.8);
    u_itec_hits->Draw("hist");
    u_itbc_hits->GetYaxis()->SetRangeUser(0,25);
    u_itbc_hits->Draw("histsame");
    u_itbc_hits->SetLineColor(kGreen+2);
    u_itbc_hits->SetLineWidth(2);
    u_otbc_hits->GetYaxis()->SetRangeUser(0,25);
    u_otbc_hits->Draw("histsame");
    u_otbc_hits->SetLineColor(kRed+2);
    u_otbc_hits->SetLineWidth(2);
    u_otec_hits->GetYaxis()->SetRangeUser(0,25);
    u_otec_hits->Draw("histsame");
    u_otec_hits->SetLineColor(kYellow+1);
    u_otec_hits->SetLineWidth(2);
    u_vbc_hits->GetYaxis()->SetRangeUser(0,25);
    u_vbc_hits->Draw("histsame");
    u_vbc_hits->SetLineColor(kRed-7);
    u_vbc_hits->SetLineWidth(2);
    u_vec_hits->GetYaxis()->SetRangeUser(0,25);
    u_vec_hits->Draw("histsame");
    u_vec_hits->SetLineColor(kBlue+1);
    u_vec_hits->SetLineWidth(2);
    u_total_hits->SetLineColor(kBlack);
    u_total_hits->SetLineWidth(2);
    u_total_hits->GetYaxis()->SetRangeUser(0,25);
    u_total_hits->Draw("histsame");

    // Create the second graph
    Double_t xmin = 0;
    Double_t xmax = 250;
    Double_t dx = (xmax - xmin) / 0.8; // 10 percent margins left and right
    Double_t ymin = 0.;
    Double_t ymax = 1.05;
    Double_t dy = (ymax - ymin) / 0.8; // 10 percent margins top and bottom
    p2->Range(xmin-0.1*dx, ymin-0.1*dy, xmax+0.1*dx, ymax+0.1*dy);
    p2->Draw();
    p2->cd();
    eff_theta->SetLineColor(kViolet-7);
    eff_theta->SetMarkerColor(kViolet-7);
    eff_theta->SetMarkerSize(2);
    eff_theta->SetLineWidth(2);
    eff_theta->Draw("psame");

    //Draw an axis on the right side
    TGaxis *axis = new TGaxis(xmax, ymin, xmax, ymax, ymin, ymax, 510, "+L");
    axis->SetLineColor(kViolet-7);
    axis->SetLabelColor(kViolet-7);
    axis->SetTitleColor(kViolet-7);
    axis->SetTitle("Efficiency");
    axis->Draw();

    TLegend* l = new TLegend(0.1, 0.3, 0.35, 0.55);
    l->AddEntry(eff_theta, "Efficiency");
    l->AddEntry(u_itec_hits, "Inner Tracker Endcap");
    l->AddEntry(u_itbc_hits, "Inner Tracker Barrel");
    l->AddEntry(u_otbc_hits, "Outer Tracker Barrel");
    l->AddEntry(u_otec_hits, "Outer Tracker Endcap");
    l->AddEntry(u_vbc_hits, "Vertex Barrel");
    l->AddEntry(u_vec_hits, "Vertex Endcap");
    l->AddEntry(u_total_hits, "Total");
    l->Draw(); 

    //c_hits->Draw();
    c_hits->SaveAs("plots/trk_hits_wogun/hits.pdf");

    outputfile->Write();
    outputfile->Close(); 
}

