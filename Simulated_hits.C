
#include <iostream>
#include <cmath>
#include <tuple>
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
#include "edm4hep/MCRecoParticleAssociationCollection.h"
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

void Simulated_hits(const char* input, const char* output, int pdg_number) {
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
    Histo_Map_theta["h_itec_hits"] = new TH1D("h_itec_hits", ";#theta [mrad]", 30, 0, 250);
    Histo_Map_theta["h_itbc_hits"] = new TH1D("h_itbc_hits", ";#theta [mrad]", 30, 0, 250);
    Histo_Map_theta["h_otbc_hits"] = new TH1D("h_otbc_hits", ";#theta [mrad]", 30, 0, 250);
    Histo_Map_theta["h_otec_hits"] = new TH1D("h_otec_hits", ";#theta [mrad]", 30, 0, 250);
    Histo_Map_theta["h_vbc_hits"] = new TH1D("h_vbc_hits", ";#theta [mrad]", 30, 0, 250);
    Histo_Map_theta["h_vec_hits"] = new TH1D("h_vec_hits", ";#theta [mrad]", 30, 0, 250);
    Histo_Map_theta["h_total_hits"] = new TH1D("h_total_hits", ";#theta [mrad]", 30, 0, 250);

    TH1D* h_pass_theta = new TH1D("h_pass_theta", ";#theta [mrad]", 30, 0, 250);
    TH1D* h_total_theta = new TH1D("h_total_theta", ";#theta [mrad]", 30, 0, 250);

    int all = 0;

    for (size_t i = 0; i < reader.getEntries("events"); ++i) {
        if (all >= 10000) break;
        all++;

        const auto event = podio::Frame(reader.readNextEntry("events"));
        const auto &mcparticles = event.get<edm4hep::MCParticleCollection>("MCParticles");
        const auto &mctruthrecolinks = event.get<edm4hep::TrackMCParticleLinkCollection>("SiTracksMCTruthLink"); //before MCRecoTrackParticleAssociation
        
        const auto &itec = event.get<edm4hep::SimTrackerHitCollection>("InnerTrackerEndcapCollection");
        const auto &itbc = event.get<edm4hep::SimTrackerHitCollection>("InnerTrackerBarrelCollection");
        const auto &otbc = event.get<edm4hep::SimTrackerHitCollection>("OuterTrackerBarrelCollection");
        const auto &otec = event.get<edm4hep::SimTrackerHitCollection>("OuterTrackerEndcapCollection");
        const auto &vbc = event.get<edm4hep::SimTrackerHitCollection>("VertexBarrelCollection");
        const auto &vec = event.get<edm4hep::SimTrackerHitCollection>("VertexEndcapCollection");

        // Map for the SimTrackerHit
        std::vector<std::tuple<std::string,const edm4hep::SimTrackerHitCollection*>> SimTrackerHit = {
            {"itec", &itec},
            {"itbc", &itbc},
            {"otbc", &otbc},
            {"otec", &otec},
            {"vbc", &vbc},
            {"vec", &vec}
        };
        
        
        // Search for first genStat 1 electron
        for (const auto& mcp : mcparticles) {
            MC_Particle particle;
            particle = getStableMCParticle(mcp, particle, pdg_number);
            if (particle.flag == true) { 
                h_total_theta->Fill(particle.theta *1000);  
                for (const auto& simtracker : SimTrackerHit) {
                    std::string name_sim;
                    const edm4hep::SimTrackerHitCollection*  hit_collection;   
                    std::tie(name_sim, hit_collection) = simtracker;                  
                    for (const auto& hit : *hit_collection) {
                        if (hit.getParticle() == mcp) { 
                            std::string histo_si_name = "h_" + name_sim + "_hits";                       
                            Histo_Map_theta[histo_si_name]->Fill(particle.theta *1000);
                            Histo_Map_theta["h_total_hits"]->Fill(particle.theta *1000);                                        
                        }
                    }
                }
                for (const auto& link : mctruthrecolinks) {
                    if (link.getSim() == mcp && link.getWeight() > 0.99) {
                        h_pass_theta->Fill(particle.theta *1000);
                    }
                }
            }
        }            
    }

    // TEfficency
    TEfficiency* eff_theta = new TEfficiency(*h_pass_theta, *h_total_theta);
    eff_theta->CreateGraph();

    // Average amount of hits
    for (auto& histo: Histo_Map_theta) {   
        histo.second->Divide(h_total_theta);
    }

    // Create canvas and draw histograms
    TLatex t;
    t.SetTextSize(30);
    t.SetTextFont(63);
    gStyle->Reset("Modern");
    TCanvas* c_hits = new TCanvas("c", "comp", 200, 10, 700, 500);
    
    TPad *p1 = new TPad("p1", "", 0, 0, 1, 1);
    TPad *p2 = new TPad("p2", "", 0, 0, 1, 1);
    p1->SetFillStyle(4000); // will be transparent
    p2->SetFillStyle(4000); 

    p1->Draw();
    p1->cd(); 
     
    // Create the first graph
    for (auto& histo: Histo_Map_theta) {   
        histo.second->SetLineWidth(2);
        histo.second->GetYaxis()->SetRangeUser(0,25);
    }
    
    Histo_Map_theta["h_itec_hits"]->GetXaxis()->SetRangeUser(0, 250);
    Histo_Map_theta["h_itec_hits"]->GetYaxis()->SetTitle(Form("Number of hits per %s", name.c_str()));
    Histo_Map_theta["h_itec_hits"]->GetYaxis()->SetTitleSize(0.05);
    Histo_Map_theta["h_itec_hits"]->GetXaxis()->SetTitleSize(0.05);
    Histo_Map_theta["h_itec_hits"]->GetXaxis()->SetTitleOffset(0.8);
    Histo_Map_theta["h_itec_hits"]->GetYaxis()->SetTitleOffset(0.8);
    Histo_Map_theta["h_itec_hits"]->SetLineColor(kCyan+1);
    Histo_Map_theta["h_itbc_hits"]->SetLineColor(kGreen+2);
    Histo_Map_theta["h_otbc_hits"]->SetLineColor(kRed+2);
    Histo_Map_theta["h_otec_hits"]->SetLineColor(kYellow+1);
    Histo_Map_theta["h_vbc_hits"]->SetLineColor(kRed-7);
    Histo_Map_theta["h_vec_hits"]->SetLineColor(kBlue+1);
    Histo_Map_theta["h_total_hits"]->SetLineColor(kBlack);
    Histo_Map_theta["h_itec_hits"]->Draw("hist");
    Histo_Map_theta["h_itbc_hits"]->Draw("histsame");
    Histo_Map_theta["h_otbc_hits"]->Draw("histsame");
    Histo_Map_theta["h_otec_hits"]->Draw("histsame");
    Histo_Map_theta["h_vbc_hits"]->Draw("histsame");  
    Histo_Map_theta["h_vec_hits"]->Draw("histsame");
    Histo_Map_theta["h_total_hits"]->Draw("histsame");

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

    // Draw an axis on the right side
    TGaxis *axis = new TGaxis(xmax, ymin, xmax, ymax, ymin, ymax, 510, "+L");
    axis->SetLineColor(kViolet-7);
    axis->SetLabelColor(kViolet-7);
    axis->SetTitleColor(kViolet-7);
    axis->SetTitle("Efficiency");
    axis->Draw();

    TLegend* l = new TLegend(0.1, 0.3, 0.35, 0.55);
    l->AddEntry(eff_theta, "Efficiency");
    l->AddEntry(Histo_Map_theta["h_itec_hits"], "Inner Tracker Endcap");
    l->AddEntry(Histo_Map_theta["h_itbc_hits"], "Inner Tracker Barrel");
    l->AddEntry(Histo_Map_theta["h_otbc_hits"], "Outer Tracker Barrel");
    l->AddEntry(Histo_Map_theta["h_otec_hits"], "Outer Tracker Endcap");
    l->AddEntry(Histo_Map_theta["h_vbc_hits"], "Vertex Barrel");
    l->AddEntry(Histo_Map_theta["h_vec_hits"], "Vertex Endcap");
    l->AddEntry(Histo_Map_theta["h_total_hits"], "Total");
    l->Draw(); 

    t.DrawLatexNDC(0.15, 0.6, Form("Single %s", name.c_str()));
    t.DrawLatexNDC(0.15, 0.93935, "CLD #font[52]{work in progress}");

    //c_hits->Draw();
    c_hits->SaveAs("plots/Tracking_efficiency/hits.pdf");

    outputfile->Write();
    outputfile->Close(); 
}
