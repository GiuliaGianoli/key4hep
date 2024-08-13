
#include <iostream>
#include <fstream>
#include <sstream>
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
#include "getParticle.h"
using namespace std;

void tracker_hits(const char* input, const char* output, int pdg_number) {
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
    std::map<std::string,TH1D*> Histo_Map_itec;
    Histo_Map_itec["h_itec_radius"] = new TH1D("h_itec_radius", ";r [m]", 30, 0.127, 2.1);
    Histo_Map_itec["h_itec_z"] = new TH1D("h_itec_z", ";z [m]", 30, 0., 2.2);
    Histo_Map_itec["h_itec_path"] = new TH1D("h_itec_path", ";Path lenght", 30, 0., 0.35);
    Histo_Map_itec["h_itec_time"] = new TH1D("h_itec_time", ";Time [ns]", 30, 0., 30);

    std::map<std::string,TH1D*> Histo_Map_itbc;
    Histo_Map_itbc["h_itbc_radius"] = new TH1D("h_itbc_hits", ";r [m]", 30, 0.127, 2.1);
    Histo_Map_itbc["h_itbc_z"] = new TH1D("h_itbc_z", ";z [m]", 30, 0., 2.2);
    Histo_Map_itbc["h_itbc_path"] = new TH1D("h_itbc_path", ";Path lenght", 30, 0., 0.35);
    Histo_Map_itbc["h_itbc_time"] = new TH1D("h_itbc_time", ";Time [ns]", 30, 0., 30);
    
    std::map<std::string,TH1D*> Histo_Map_otbc;
    Histo_Map_otbc["h_otbc_radius"] = new TH1D("h_otbc_hits", ";r [m]", 30, 0.127, 2.1);
    Histo_Map_otbc["h_otbc_z"] = new TH1D("h_otbc_z", ";z [m]", 30, 0., 2.2);
    Histo_Map_otbc["h_otbc_path"] = new TH1D("h_otbc_path", ";Path lenght", 30, 0., 0.35);
    Histo_Map_otbc["h_otbc_time"] = new TH1D("h_otbc_time", ";Time [ns]", 30, 0., 30);
    
    std::map<std::string,TH1D*> Histo_Map_otec;
    Histo_Map_otec["h_otec_radius"] = new TH1D("h_otec_hits", ";r [m]", 30, 0.127, 2.1);
    Histo_Map_otec["h_otec_z"] = new TH1D("h_otec_z", ";z [m]", 30, 0., 2.2);
    Histo_Map_otec["h_otec_path"] = new TH1D("h_otec_path", ";Path lenght", 30, 0., 0.35);
    Histo_Map_otec["h_otec_time"] = new TH1D("h_otec_time", ";Time [ns]", 30, 0., 30);
    
    std::map<std::string,TH1D*> Histo_Map_vbc;
    Histo_Map_vbc["h_vbc_radius"] = new TH1D("h_vbc_hits", ";r [m]", 30, 0.127, 2.1);
    Histo_Map_vbc["h_vbc_z"] = new TH1D("h_vbc_z", ";z [m]", 30, 0., 2.2);
    Histo_Map_vbc["h_vbc_path"] = new TH1D("h_vbc_path", ";Path lenght", 30, 0., 0.35);
    Histo_Map_vbc["h_vbc_time"] = new TH1D("h_vbc_time", ";Time [ns]", 30, 0., 30);
    
    std::map<std::string,TH1D*> Histo_Map_vec;
    Histo_Map_vec["h_vec_radius"] = new TH1D("h_vec_hits", ";r [m]", 30, 0.127, 2.1);
    Histo_Map_vec["h_vec_z"] = new TH1D("h_vec_z", ";z [m]", 30, 0., 2.2);
    Histo_Map_vec["h_vec_path"] = new TH1D("h_vec_path", ";Path lenght", 30, 0., 0.35);
    Histo_Map_vec["h_vec_time"] = new TH1D("h_vec_time", ";Time [ns]", 30, 0., 30);
    
    std::map<std::string,TH1D*> Histo_Map_total;
    Histo_Map_total["h_total_radius"] = new TH1D("h_total_hits", ";r [m]", 30, 0.127, 2.1);
    Histo_Map_total["h_total_z"] = new TH1D("h_total_z", ";z [m]", 30, 0., 2.2);
    Histo_Map_total["h_total_path"] = new TH1D("h_total_path", ";Path lenght", 30, 0., 0.35);
    Histo_Map_total["h_total_time"] = new TH1D("h_total_time", ";Time [ns]", 30, 0., 30);

    int all = 0;

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
            if (abs(mcp.getPDG()) == pdg_number && mcp.getGeneratorStatus() == 1) {                      
                for (const auto& hit : itec) {
                    Tracker_Hit tracker_hit;
                    tracker_hit = getHit(mcp, hit, tracker_hit);
                    if(tracker_hit.flag == true) {                 
                        Histo_Map_itec["h_itec_radius"]->Fill(tracker_hit.radius/1000);
                        Histo_Map_itec["h_itec_z"]->Fill(tracker_hit.z/1000);
                        Histo_Map_itec["h_itec_path"]->Fill(tracker_hit.path);
                        Histo_Map_itec["h_itec_time"]->Fill(tracker_hit.time);
                        Histo_Map_total["h_total_radius"]->Fill(tracker_hit.radius/1000);  
                        Histo_Map_total["h_total_z"]->Fill(tracker_hit.z/1000);
                        Histo_Map_total["h_total_path"]->Fill(tracker_hit.path);
                        Histo_Map_total["h_total_time"]->Fill(tracker_hit.time);     
                    }                                                    
                }
                for (const auto& hit : itbc) {
                    Tracker_Hit tracker_hit;
                    tracker_hit = getHit(mcp, hit, tracker_hit);
                    if(tracker_hit.flag == true) {                        
                        Histo_Map_itbc["h_itbc_radius"]->Fill(tracker_hit.radius/1000);
                        Histo_Map_itbc["h_itbc_z"]->Fill(tracker_hit.z/1000);
                        Histo_Map_itbc["h_itbc_path"]->Fill(tracker_hit.path);
                        Histo_Map_itbc["h_itbc_time"]->Fill(tracker_hit.time);                        
                        Histo_Map_total["h_total_radius"]->Fill(tracker_hit.radius/1000);   
                        Histo_Map_total["h_total_z"]->Fill(tracker_hit.z/1000);
                        Histo_Map_total["h_total_path"]->Fill(tracker_hit.path);
                        Histo_Map_total["h_total_time"]->Fill(tracker_hit.time);                    
                    }
                }
                for (const auto& hit : otbc) {
                    Tracker_Hit tracker_hit;
                    tracker_hit = getHit(mcp, hit, tracker_hit);
                    if(tracker_hit.flag == true) {                   
                        Histo_Map_otbc["h_otbc_radius"]->Fill(tracker_hit.radius/1000);
                        Histo_Map_otbc["h_otbc_z"]->Fill(tracker_hit.z/1000);
                        Histo_Map_otbc["h_otbc_path"]->Fill(tracker_hit.path);
                        Histo_Map_otbc["h_otbc_time"]->Fill(tracker_hit.time);                        
                        Histo_Map_total["h_total_radius"]->Fill(tracker_hit.radius/1000);   
                        Histo_Map_total["h_total_z"]->Fill(tracker_hit.z/1000);
                        Histo_Map_total["h_total_path"]->Fill(tracker_hit.path);
                        Histo_Map_total["h_total_time"]->Fill(tracker_hit.time);                    
                    }
                }
                for (const auto& hit : otec) {
                    Tracker_Hit tracker_hit;
                    tracker_hit = getHit(mcp, hit, tracker_hit);
                    if(tracker_hit.flag == true) {                       
                        Histo_Map_otec["h_otec_radius"]->Fill(tracker_hit.radius/1000);
                        Histo_Map_otec["h_otec_z"]->Fill(tracker_hit.z/1000);
                        Histo_Map_otec["h_otec_path"]->Fill(tracker_hit.path);
                        Histo_Map_otec["h_otec_time"]->Fill(tracker_hit.time);                        
                        Histo_Map_total["h_total_radius"]->Fill(tracker_hit.radius/1000);   
                        Histo_Map_total["h_total_z"]->Fill(tracker_hit.z/1000);
                        Histo_Map_total["h_total_path"]->Fill(tracker_hit.path);
                        Histo_Map_total["h_total_time"]->Fill(tracker_hit.time);                        
                    }
                }
                for (const auto& hit : vbc) {
                    Tracker_Hit tracker_hit;
                    tracker_hit = getHit(mcp, hit, tracker_hit);
                    if(tracker_hit.flag == true) {                       
                        Histo_Map_vbc["h_vbc_radius"]->Fill(tracker_hit.radius/1000);
                        Histo_Map_vbc["h_vbc_z"]->Fill(tracker_hit.z/1000);
                        Histo_Map_vbc["h_vbc_path"]->Fill(tracker_hit.path);
                        Histo_Map_vbc["h_vbc_time"]->Fill(tracker_hit.time);                        
                        Histo_Map_total["h_total_radius"]->Fill(tracker_hit.radius/1000);   
                        Histo_Map_total["h_total_z"]->Fill(tracker_hit.z/1000);
                        Histo_Map_total["h_total_path"]->Fill(tracker_hit.path);
                        Histo_Map_total["h_total_time"]->Fill(tracker_hit.time);                       
                    }
                }
                for (const auto& hit : vec) {
                    Tracker_Hit tracker_hit;
                    tracker_hit = getHit(mcp, hit, tracker_hit);
                    if(tracker_hit.flag == true) {                       
                        Histo_Map_vec["h_vec_radius"]->Fill(tracker_hit.radius/1000);
                        Histo_Map_vec["h_vec_z"]->Fill(tracker_hit.z/1000);
                        Histo_Map_vec["h_vec_path"]->Fill(tracker_hit.path);
                        Histo_Map_vec["h_vec_time"]->Fill(tracker_hit.time);                        
                        Histo_Map_total["h_total_radius"]->Fill(tracker_hit.radius/1000);   
                        Histo_Map_total["h_total_z"]->Fill(tracker_hit.z/1000);
                        Histo_Map_total["h_total_path"]->Fill(tracker_hit.path);
                        Histo_Map_total["h_total_time"]->Fill(tracker_hit.time);                       
                    }
                }
            }
        }            
    }
    
    // Create graphs
    for (auto& histo: Histo_Map_itec) {   
        histo.second->SetLineWidth(2);
        histo.second->SetLineColor(kCyan+1);
        histo.second->GetYaxis()->SetTitle("Total number of hits");
        histo.second->GetYaxis()->SetTitleSize(0.05);
        histo.second->GetXaxis()->SetTitleSize(0.05);
        histo.second->GetXaxis()->SetTitleOffset(0.8);
        histo.second->GetXaxis()->SetTitleOffset(0.8);
    }
    Histo_Map_itec["h_itec_radius"]->GetXaxis()->SetRangeUser(0, 2.1);
    Histo_Map_itec["h_itec_z"]->GetXaxis()->SetRangeUser(0, 2.2);
    Histo_Map_itec["h_itec_path"]->GetXaxis()->SetRangeUser(0, 0.35);
    Histo_Map_itec["h_itec_time"]->GetXaxis()->SetRangeUser(0, 30);
    Histo_Map_itec["h_itec_radius"]->GetYaxis()->SetRangeUser(0, 120000);
    Histo_Map_itec["h_itec_z"]->GetYaxis()->SetRangeUser(0, 300000);
    Histo_Map_itec["h_itec_path"]->GetYaxis()->SetRangeUser(0, 300000);
    Histo_Map_itec["h_itec_time"]->GetYaxis()->SetRangeUser(0, 800000);
    for (auto& histo: Histo_Map_itbc) {   
        histo.second->SetLineWidth(2);
        histo.second->SetLineColor(kGreen+2);
    }
    Histo_Map_itbc["h_itbc_radius"]->GetXaxis()->SetRangeUser(0, 2.1);
    Histo_Map_itbc["h_itbc_z"]->GetXaxis()->SetRangeUser(0, 2.2);
    Histo_Map_itbc["h_itbc_path"]->GetXaxis()->SetRangeUser(0, 0.35);
    Histo_Map_itbc["h_itbc_time"]->GetXaxis()->SetRangeUser(0, 30);
    for (auto& histo: Histo_Map_otbc) {   
        histo.second->SetLineWidth(2);
        histo.second->SetLineColor(kRed+2);
    }
    Histo_Map_otbc["h_otbc_radius"]->GetXaxis()->SetRangeUser(0, 2.1);
    Histo_Map_otbc["h_otbc_z"]->GetXaxis()->SetRangeUser(0, 2.2);
    Histo_Map_otbc["h_otbc_path"]->GetXaxis()->SetRangeUser(0, 0.35);
    Histo_Map_otbc["h_otbc_time"]->GetXaxis()->SetRangeUser(0, 30);
    for (auto& histo: Histo_Map_otec) {   
        histo.second->SetLineWidth(2);
        histo.second->SetLineColor(kYellow+1);
    }
    Histo_Map_otec["h_otec_radius"]->GetXaxis()->SetRangeUser(0, 2.1);
    Histo_Map_otec["h_otec_z"]->GetXaxis()->SetRangeUser(0, 2.2);
    Histo_Map_otec["h_otec_path"]->GetXaxis()->SetRangeUser(0, 0.35);
    Histo_Map_otec["h_otec_time"]->GetXaxis()->SetRangeUser(0, 30);
    for (auto& histo: Histo_Map_vbc) {   
        histo.second->SetLineWidth(2);
        histo.second->SetLineColor(kRed-7);
    }
    Histo_Map_vbc["h_vbc_radius"]->GetXaxis()->SetRangeUser(0, 2.1);
    Histo_Map_vbc["h_vbc_z"]->GetXaxis()->SetRangeUser(0, 2.2);
    Histo_Map_vbc["h_vbc_path"]->GetXaxis()->SetRangeUser(0, 0.35);
    Histo_Map_vbc["h_vbc_time"]->GetXaxis()->SetRangeUser(0, 30);
    for (auto& histo: Histo_Map_vec) {   
        histo.second->SetLineWidth(2);
        histo.second->SetLineColor(kBlue+1);
    }
    Histo_Map_vec["h_vec_radius"]->GetXaxis()->SetRangeUser(0, 2.1);
    Histo_Map_vec["h_vec_z"]->GetXaxis()->SetRangeUser(0, 2.2);
    Histo_Map_vec["h_vec_path"]->GetXaxis()->SetRangeUser(0, 0.35);
    Histo_Map_vec["h_vec_time"]->GetXaxis()->SetRangeUser(0, 30);
    for (auto& histo: Histo_Map_total) {   
        histo.second->SetLineWidth(2);
        histo.second->SetLineColor(kBlack);
    }
    Histo_Map_total["h_total_radius"]->GetXaxis()->SetRangeUser(0, 2.1);
    Histo_Map_total["h_total_z"]->GetXaxis()->SetRangeUser(0, 2.2);
    Histo_Map_total["h_total_path"]->GetXaxis()->SetRangeUser(0, 0.35);
    Histo_Map_total["h_total_time"]->GetXaxis()->SetRangeUser(0, 30);

    // Create canvas and draw histograms
    TLatex t;
    t.SetTextSize(30);
    t.SetTextFont(63);
    gStyle->Reset("Modern");

    TCanvas* c_hits_radius = new TCanvas("c", "comp", 200, 10, 700, 500);

    Histo_Map_itec["h_itec_radius"]->Draw("hist");
    Histo_Map_itbc["h_itbc_radius"]->Draw("histsame"); 
    Histo_Map_otbc["h_otbc_radius"]->Draw("histsame");  
    Histo_Map_otec["h_otec_radius"]->Draw("histsame"); 
    Histo_Map_vbc["h_vbc_radius"]->Draw("histsame");  
    Histo_Map_vec["h_vec_radius"]->Draw("histsame");
    Histo_Map_total["h_total_radius"]->Draw("histsame");

    TLegend* l = new TLegend(0.5, 0.65, 0.75, 0.9);
    l->AddEntry(Histo_Map_itec["h_itec_radius"], "Inner Tracker Endcap");
    l->AddEntry(Histo_Map_itbc["h_itbc_radius"], "Inner Tracker Barrel");
    l->AddEntry(Histo_Map_otbc["h_otbc_radius"], "Outer Tracker Barrel");
    l->AddEntry(Histo_Map_otec["h_otec_radius"], "Outer Tracker Endcap");
    l->AddEntry(Histo_Map_vbc["h_vbc_radius"], "Vertex Barrel");
    l->AddEntry(Histo_Map_vec["h_vec_radius"], "Vertex Endcap");
    l->AddEntry(Histo_Map_total["h_total_radius"], "Total");
    l->Draw(); 

    t.DrawLatexNDC(0.5, 0.4, Form("Single %s", name.c_str()));
    t.DrawLatexNDC(0.15, 0.93935, "CLD #font[52]{work in progress}");

    c_hits_radius->SaveAs("plots/trk_hits_wogun/hits_r.pdf");

    TCanvas* c_hits_z = new TCanvas("c2", "comp", 200, 10, 700, 500);

    Histo_Map_itec["h_itec_z"]->Draw("hist");
    Histo_Map_itbc["h_itbc_z"]->Draw("histsame"); 
    Histo_Map_otbc["h_otbc_z"]->Draw("histsame");  
    Histo_Map_otec["h_otec_z"]->Draw("histsame"); 
    Histo_Map_vbc["h_vbc_z"]->Draw("histsame");  
    Histo_Map_vec["h_vec_z"]->Draw("histsame");
    Histo_Map_total["h_total_z"]->Draw("histsame");

    TLegend* l2 = new TLegend(0.5, 0.65, 0.75, 0.9);
    l2->AddEntry(Histo_Map_itec["h_itec_z"], "Inner Tracker Endcap");
    l2->AddEntry(Histo_Map_itbc["h_itbc_z"], "Inner Tracker Barrel");
    l2->AddEntry(Histo_Map_otbc["h_otbc_z"], "Outer Tracker Barrel");
    l2->AddEntry(Histo_Map_otec["h_otec_z"], "Outer Tracker Endcap");
    l2->AddEntry(Histo_Map_vbc["h_vbc_z"], "Vertex Barrel");
    l2->AddEntry(Histo_Map_vec["h_vec_z"], "Vertex Endcap");
    l2->AddEntry(Histo_Map_total["h_total_z"], "Total");
    l2->Draw(); 

    t.DrawLatexNDC(0.5, 0.4, Form("Single %s", name.c_str()));
    t.DrawLatexNDC(0.15, 0.93935, "CLD #font[52]{work in progress}");

    c_hits_z->SaveAs("plots/trk_hits_wogun/hits_z.pdf");

    TCanvas* c_hits_path = new TCanvas("c3", "comp", 200, 10, 700, 500);

    Histo_Map_itec["h_itec_path"]->Draw("hist");
    Histo_Map_itbc["h_itbc_path"]->Draw("histsame"); 
    Histo_Map_otbc["h_otbc_path"]->Draw("histsame");  
    Histo_Map_otec["h_otec_path"]->Draw("histsame"); 
    Histo_Map_vbc["h_vbc_path"]->Draw("histsame");  
    Histo_Map_vec["h_vec_path"]->Draw("histsame");
    Histo_Map_total["h_total_path"]->Draw("histsame");

    TLegend* l3 = new TLegend(0.5, 0.65, 0.75, 0.9);
    l3->AddEntry(Histo_Map_itec["h_itec_path"], "Inner Tracker Endcap");
    l3->AddEntry(Histo_Map_itbc["h_itbc_path"], "Inner Tracker Barrel");
    l3->AddEntry(Histo_Map_otbc["h_otbc_path"], "Outer Tracker Barrel");
    l3->AddEntry(Histo_Map_otec["h_otec_path"], "Outer Tracker Endcap");
    l3->AddEntry(Histo_Map_vbc["h_vbc_path"], "Vertex Barrel");
    l3->AddEntry(Histo_Map_vec["h_vec_path"], "Vertex Endcap");
    l3->AddEntry(Histo_Map_total["h_total_path"], "Total");
    l3->Draw(); 

    t.DrawLatexNDC(0.5, 0.4, Form("Single %s", name.c_str()));
    t.DrawLatexNDC(0.15, 0.93935, "CLD #font[52]{work in progress}");

    c_hits_path->SaveAs("plots/trk_hits_wogun/hits_path.pdf");

    TCanvas* c_hits_time = new TCanvas("c4", "comp", 200, 10, 700, 500);

    Histo_Map_itec["h_itec_time"]->Draw("hist");
    Histo_Map_itbc["h_itbc_time"]->Draw("histsame"); 
    Histo_Map_otbc["h_otbc_time"]->Draw("histsame");  
    Histo_Map_otec["h_otec_time"]->Draw("histsame"); 
    Histo_Map_vbc["h_vbc_time"]->Draw("histsame");  
    Histo_Map_vec["h_vec_time"]->Draw("histsame");
    Histo_Map_total["h_total_time"]->Draw("histsame");

    TLegend* l4 = new TLegend(0.5, 0.65, 0.75, 0.9);
    l4->AddEntry(Histo_Map_itec["h_itec_time"], "Inner Tracker Endcap");
    l4->AddEntry(Histo_Map_itbc["h_itbc_time"], "Inner Tracker Barrel");
    l4->AddEntry(Histo_Map_otbc["h_otbc_time"], "Outer Tracker Barrel");
    l4->AddEntry(Histo_Map_otec["h_otec_time"], "Outer Tracker Endcap");
    l4->AddEntry(Histo_Map_vbc["h_vbc_time"], "Vertex Barrel");
    l4->AddEntry(Histo_Map_vec["h_vec_time"], "Vertex Endcap");
    l4->AddEntry(Histo_Map_total["h_total_time"], "Total");
    l4->Draw(); 

    t.DrawLatexNDC(0.5, 0.4, Form("Single %s", name.c_str()));
    t.DrawLatexNDC(0.15, 0.93935, "CLD #font[52]{work in progress}");

    c_hits_time->SaveAs("plots/trk_hits_wogun/hits_time.pdf");

    outputfile->Write();
    outputfile->Close(); 
}

