
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
#include "edm4hep/MCRecoParticleAssociation.h"
#include "edm4hep/ReconstructedParticleCollection.h"
#include "edm4hep/TrackerHitSimTrackerHitLink.h"
#include "edm4hep/TrackMCParticleLink.h"
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
    std::map<std::string,std::map<std::string,TH1D*>> Histo_Maps;
    Histo_Maps["itec"] = {
    {"h_itec_radius", new TH1D("h_itec_radius", ";r [m]", 30, 0.127, 2.1)},
    {"h_itec_z",  new TH1D("h_itec_z", ";z [m]", 30, 0., 2.2)},
    {"h_itec_path", new TH1D("h_itec_path", ";Path lenght", 30, 0., 0.35)},
    {"h_itec_eDep", new TH1D("h_itec_eDep", ";E dep [MeV]", 30, 0., 0.25)},
    {"h_itec_pt", new TH1D("h_itec_pt", "; p_{T} [GeV]", 30, 0., 100)},
    {"h_itec_p", new TH1D("h_itec_p", ";p [GeV]", 30, 0., 100)},
    {"h_itec_time",new TH1D("h_itec_time", ";Time [ns]", 30, 0., 30)}
    };
    Histo_Maps["itbc"] = {
    {"h_itbc_radius", new TH1D("h_itbc_radius", ";r [m]", 30, 0.127, 2.1)},
    {"h_itbc_z",  new TH1D("h_itbc_z", ";z [m]", 30, 0., 2.2)},
    {"h_itbc_path", new TH1D("h_itbc_path", ";Path lenght", 30, 0., 0.35)},
    {"h_itbc_eDep", new TH1D("h_itbc_eDep", ";E dep [MeV]", 30, 0., 0.25)},
    {"h_itbc_pt", new TH1D("h_itbc_pt", "; p_{T} [GeV]", 30, 0., 100)},
    {"h_itbc_p", new TH1D("h_itbc_p", ";p [GeV]", 30, 0., 100)},
    {"h_itbc_time",new TH1D("h_itbc_time", ";Time [ns]", 30, 0., 30)}
    };
    Histo_Maps["otbc"] = {
    {"h_otbc_radius", new TH1D("h_otbc_radius", ";r [m]", 30, 0.127, 2.1)},
    {"h_otbc_z",  new TH1D("h_otbc_z", ";z [m]", 30, 0., 2.2)},
    {"h_otbc_path", new TH1D("h_otbc_path", ";Path lenght", 30, 0., 0.35)},
    {"h_otbc_eDep", new TH1D("h_otbc_eDep", ";E dep [MeV]", 30, 0., 0.25)},
    {"h_otbc_pt", new TH1D("h_otbc_pt", "; p_{T} [GeV]", 30, 0., 100)},
    {"h_otbc_p", new TH1D("h_otbc_p", ";p [GeV]", 30, 0., 100)},
    {"h_otbc_time",new TH1D("h_otbc_time", ";Time [ns]", 30, 0., 30)}
    };
    Histo_Maps["otec"] = {
    {"h_otec_radius", new TH1D("h_otec_radius", ";r [m]", 30, 0.127, 2.1)},
    {"h_otec_z",  new TH1D("h_otec_z", ";z [m]", 30, 0., 2.2)},
    {"h_otec_path", new TH1D("h_otec_path", ";Path lenght", 30, 0., 0.35)},
    {"h_otec_eDep", new TH1D("h_otec_eDep", ";E dep [MeV]", 30, 0., 0.25)},
    {"h_otec_pt", new TH1D("h_otec_pt", "; p_{T} [GeV]", 30, 0., 100)},
    {"h_otec_p", new TH1D("h_otec_p", ";p [GeV]", 30, 0., 100)},
    {"h_otec_time",new TH1D("h_otec_time", ";Time [ns]", 30, 0., 30)}
    };
    Histo_Maps["vbc"] = {
    {"h_vbc_radius", new TH1D("h_vbc_radius", ";r [m]", 30, 0.127, 2.1)},
    {"h_vbc_z",  new TH1D("h_vbc_z", ";z [m]", 30, 0., 2.2)},
    {"h_vbc_path", new TH1D("h_vbc_path", ";Path lenght", 30, 0., 0.35)},
    {"h_vbc_eDep", new TH1D("h_vbc_eDep", ";E dep [MeV]", 30, 0., 0.25)},
    {"h_vbc_pt", new TH1D("h_vbc_pt", "; p_{T} [GeV]", 30, 0., 100)},
    {"h_vbc_p", new TH1D("h_vbc_p", ";p [GeV]", 30, 0., 100)},
    {"h_vbc_time",new TH1D("h_vbc_time", ";Time [ns]", 30, 0., 30)}
    };
    Histo_Maps["vec"] = {
    {"h_vec_radius", new TH1D("h_vec_radius", ";r [m]", 30, 0.127, 2.1)},
    {"h_vec_z",  new TH1D("h_vec_z", ";z [m]", 30, 0., 2.2)},
    {"h_vec_path", new TH1D("h_vec_path", ";Path lenght", 30, 0., 0.35)},
    {"h_vec_eDep", new TH1D("h_vec_eDep", ";E dep [MeV]", 30, 0., 0.25)},
    {"h_vec_pt", new TH1D("h_vec_pt", "; p_{T} [GeV]", 30, 0., 100)},
    {"h_vec_p", new TH1D("h_vec_p", ";p [GeV]", 30, 0., 100)},
    {"h_vec_time",new TH1D("h_vec_time", ";Time [ns]", 30, 0., 30)}
    };
  
    std::map<std::string,TH1D*> Histo_Map_total = {
    {"h_total_radius", new TH1D("h_total_hits", ";r [m]", 30, 0.127, 2.1)},
    {"h_total_z", new TH1D("h_total_z", ";z [m]", 30, 0., 2.2)},
    {"h_total_path",new TH1D("h_total_path", ";Path lenght", 30, 0., 0.35)},
    {"h_total_eDep", new TH1D("h_total_eDep", ";E dep [MeV]", 30, 0., 0.25)},
    {"h_total_pt", new TH1D("h_total_pt", "; p_{T} [GeV]", 30, 0., 100)},
    {"h_total_p", new TH1D("h_total_p", ";p [GeV]", 30, 0., 100)},
    {"h_total_time", new TH1D("h_total_time", ";Time [ns]", 30, 0., 30)}
    };

    TH1D* h_pass_pt = new TH1D("h_pass_pt", ";p_{T} [GeV]", 30, 0, 100);
    TH1D* h_eff_pt = new TH1D("h_eff_pt", ";p_{T} [GeV]", 30, 0, 100);

    TH1D* h_pass_p = new TH1D("h_pass_p", ";p [GeV]", 30, 0, 100);
    TH1D* h_eff_p = new TH1D("h_eff_p", ";p [GeV]", 30, 0, 100);

    int all = 0;

    for (size_t i = 0; i < reader.getEntries("events"); ++i) {
        if (all >= 100000) break;
        all++;

        const auto event = podio::Frame(reader.readNextEntry("events"));
        const auto &mcparticles = event.get<edm4hep::MCParticleCollection>("MCParticles");
        const auto& mclinks = event.get<edm4hep::MCRecoParticleAssociationCollection>("MCTruthRecoLink");
        const auto &mctruthrecolinks = event.get<edm4hep::TrackMCParticleLinkCollection>("SiTracksMCTruthLink");

        const auto &itec = event.get<edm4hep::SimTrackerHitCollection>("InnerTrackerEndcapCollection");
        const auto &itbc = event.get<edm4hep::SimTrackerHitCollection>("InnerTrackerBarrelCollection");
        const auto &otbc = event.get<edm4hep::SimTrackerHitCollection>("OuterTrackerBarrelCollection");
        const auto &otec = event.get<edm4hep::SimTrackerHitCollection>("OuterTrackerEndcapCollection");
        const auto &vbc = event.get<edm4hep::SimTrackerHitCollection>("VertexBarrelCollection");
        const auto &vec = event.get<edm4hep::SimTrackerHitCollection>("VertexEndcapCollection");

        const auto &itec_link = event.get<edm4hep::TrackerHitSimTrackerHitLinkCollection>("InnerTrackerEndcapHitsRelations");
        const auto &itbc_link = event.get<edm4hep::TrackerHitSimTrackerHitLinkCollection>("InnerTrackerBarrelHitsRelations");
        const auto &otbc_link = event.get<edm4hep::TrackerHitSimTrackerHitLinkCollection>("OuterTrackerBarrelHitsRelations");
        const auto &otec_link = event.get<edm4hep::TrackerHitSimTrackerHitLinkCollection>("OuterTrackerEndcapHitsRelations");
        const auto &vbc_link = event.get<edm4hep::TrackerHitSimTrackerHitLinkCollection>("VXDTrackerHitRelations");
        const auto &vec_link = event.get<edm4hep::TrackerHitSimTrackerHitLinkCollection>("VXDEndcapTrackerHitRelations"); 

        // Map for the SimTrackerHit
        // Map for the SimTrackerHit
        std::vector<std::tuple<std::string,const edm4hep::SimTrackerHitCollection*,const edm4hep::TrackerHitSimTrackerHitLinkCollection*>> SimTrackerHit = {
            {"itec", &itec, &itec_link},
            {"itbc", &itbc, &itbc_link},
            {"otbc", &otbc, &otbc_link},
            {"otec", &otec, &otec_link},
            {"vbc", &vbc, &vbc_link},
            {"vec", &vec, &vec_link}
        };
        
        // Search for first genStat 1 electron
        for (const auto& mcp : mcparticles) {
            MC_Particle particle;
            particle = getStableMCParticle(mcp, particle, pdg_number);
            if (particle.flag == true) { 
                h_eff_pt->Fill(particle.pt);
                h_eff_p->Fill(particle.p);
                for (const auto& link : mctruthrecolinks) {
                    if (link.getSim() == mcp && link.getWeight() > 0.99) {
                        h_pass_pt->Fill(particle.pt);
                        h_pass_p->Fill(particle.p);
                    }
                }
                for (const auto& simtracker : SimTrackerHit) {
                    std::string name_sim;
                    const edm4hep::SimTrackerHitCollection*  hit_collection;  
                    const edm4hep::TrackerHitSimTrackerHitLinkCollection* hit_link;  
                    std::tie(name_sim, hit_collection, hit_link) = simtracker;    
                    for (const auto& hit : *hit_collection) {
                        if (hit.getParticle() == mcp) { 
                            Tracker_Hit tracker_hit;
                            tracker_hit = getHit(mcp, hit, tracker_hit);
                            Reconstructed_Particle reco_p;
                            reco_p = getReco(mcp, mclinks, reco_p); 
                            auto tracks = reco_p.reco.first.getTracks();
                            bool flag;
                            flag = getNumberofHitReco(tracks, hit, *hit_link); 
                            bool flag_secondary = hit.isProducedBySecondary();
                            if (flag == false && flag_secondary == false){    
                                auto& Histo_Map_current = Histo_Maps[name_sim];  
                                Histo_Map_current["h_" + name_sim + "_radius"]->Fill(tracker_hit.radius/1000);
                                Histo_Map_current["h_" + name_sim + "_z"]->Fill(tracker_hit.z/1000);
                                Histo_Map_current["h_" + name_sim + "_path"]->Fill(tracker_hit.path);
                                Histo_Map_current["h_" + name_sim + "_eDep"]->Fill(tracker_hit.eDep*1000);
                                Histo_Map_current["h_" + name_sim + "_pt"]->Fill(tracker_hit.pt);
                                Histo_Map_current["h_" + name_sim + "_p"]->Fill(tracker_hit.p);
                                Histo_Map_current["h_" + name_sim + "_time"]->Fill(tracker_hit.time);
                                Histo_Map_total["h_total_radius"]->Fill(tracker_hit.radius/1000);  
                                Histo_Map_total["h_total_z"]->Fill(tracker_hit.z/1000);
                                Histo_Map_total["h_total_path"]->Fill(tracker_hit.path);
                                Histo_Map_total["h_total_eDep"]->Fill(tracker_hit.eDep*1000);
                                Histo_Map_total["h_total_pt"]->Fill(tracker_hit.pt);
                                Histo_Map_total["h_total_p"]->Fill(tracker_hit.p);
                                Histo_Map_total["h_total_time"]->Fill(tracker_hit.time);
                            }
                        }     
                    }                                                    
                }              
            }
        }            
    }
    
    // TEfficency
    TEfficiency* eff_pt = new TEfficiency(*h_pass_pt, *h_eff_pt);
    eff_pt->CreateGraph();

    TEfficiency* eff_p = new TEfficiency(*h_pass_p, *h_eff_p);
    eff_p->CreateGraph();

    // Create graphs
    for(auto& histo_map_pair : Histo_Maps) {
        auto& histo_map = histo_map_pair.second;
        int color = kBlack; 
        if (histo_map_pair.first == "itec") color = kCyan + 1;
        else if (histo_map_pair.first == "itbc") color = kGreen + 2;
        else if (histo_map_pair.first == "otbc") color = kRed + 2;
        else if (histo_map_pair.first == "otec") color = kYellow + 1;
        else if (histo_map_pair.first == "vbc") color = kRed - 7;
        else if (histo_map_pair.first == "vec") color = kBlue + 1;
        for (auto& histo: histo_map) {  
            histo.second->SetLineWidth(2); 
            histo.second->SetLineColor(color);
            histo.second->GetYaxis()->SetTitle("Total number of hits");
            histo.second->GetYaxis()->SetTitleSize(0.05);
            histo.second->GetXaxis()->SetTitleSize(0.05);
            histo.second->GetXaxis()->SetTitleOffset(0.8);
            histo.second->GetYaxis()->SetTitleOffset(0.85);
        }
    }
    Histo_Maps["itec"]["h_itec_radius"]->GetXaxis()->SetRangeUser(0, 2.1);
    Histo_Maps["itec"]["h_itec_z"]->GetXaxis()->SetRangeUser(0, 2.2);
    Histo_Maps["itec"]["h_itec_path"]->GetXaxis()->SetRangeUser(0, 0.35);
    Histo_Maps["itec"]["h_itec_eDep"]->GetXaxis()->SetRangeUser(0, 0.25);
    Histo_Maps["itec"]["h_itec_pt"]->GetXaxis()->SetRangeUser(0, 100);
    Histo_Maps["itec"]["h_itec_p"]->GetXaxis()->SetRangeUser(0, 100);
    Histo_Maps["itec"]["h_itec_time"]->GetXaxis()->SetRangeUser(0, 30);
    Histo_Maps["itec"]["h_itec_radius"]->GetYaxis()->SetRangeUser(0, 3500);
    Histo_Maps["itec"]["h_itec_z"]->GetYaxis()->SetRangeUser(0, 2000);
    Histo_Maps["itec"]["h_itec_path"]->GetYaxis()->SetRangeUser(0, 3000);
    Histo_Maps["itec"]["h_itec_eDep"]->GetYaxis()->SetRangeUser(0, 3000);
    Histo_Maps["itec"]["h_itec_pt"]->GetYaxis()->SetRangeUser(0, 3000);
    Histo_Maps["itec"]["h_itec_p"]->GetYaxis()->SetRangeUser(0, 3000);
    Histo_Maps["itec"]["h_itec_time"]->GetYaxis()->SetRangeUser(0, 3500);

    for (auto& histo : Histo_Map_total) {
        histo.second->SetLineWidth(2);
        histo.second->SetLineColor(kBlack);
    }
    Histo_Map_total["h_total_radius"]->GetXaxis()->SetRangeUser(0, 2.1);
    Histo_Map_total["h_total_z"]->GetXaxis()->SetRangeUser(0, 2.2);
    Histo_Map_total["h_total_path"]->GetXaxis()->SetRangeUser(0, 0.35);
    Histo_Map_total["h_total_eDep"]->GetXaxis()->SetRangeUser(0, 0.25);
    Histo_Map_total["h_total_pt"]->GetXaxis()->SetRangeUser(0, 100);
    Histo_Map_total["h_total_p"]->GetXaxis()->SetRangeUser(0, 100);
    Histo_Map_total["h_total_time"]->GetXaxis()->SetRangeUser(0, 30);

    // Create canvas and draw histograms
    TLatex t;
    t.SetTextSize(30);
    t.SetTextFont(63);
    gStyle->Reset("Modern");

    TCanvas* c_hits_radius = new TCanvas("c_hits_radius", "comp", 200, 10, 700, 500);
    Histo_Maps["itec"]["h_itec_radius"]->Draw("hist");
    Histo_Maps["itbc"]["h_itbc_radius"]->Draw("histsame");
    Histo_Maps["otbc"]["h_otbc_radius"]->Draw("histsame");
    Histo_Maps["otec"]["h_otec_radius"]->Draw("histsame");
    Histo_Maps["vbc"]["h_vbc_radius"]->Draw("histsame");
    Histo_Maps["vec"]["h_vec_radius"]->Draw("histsame");
    Histo_Map_total["h_total_radius"]->Draw("histsame");

    TLegend* l = new TLegend(0.5, 0.65, 0.75, 0.9);
    l->AddEntry(Histo_Maps["itec"]["h_itec_radius"], "Inner Tracker Endcap");
    l->AddEntry(Histo_Maps["itbc"]["h_itbc_radius"], "Inner Tracker Barrel");
    l->AddEntry(Histo_Maps["otbc"]["h_otbc_radius"], "Outer Tracker Barrel");
    l->AddEntry(Histo_Maps["otec"]["h_otec_radius"], "Outer Tracker Endcap");
    l->AddEntry(Histo_Maps["vbc"]["h_vbc_radius"], "Vertex Barrel");
    l->AddEntry(Histo_Maps["vec"]["h_vec_radius"], "Vertex Endcap");
    l->AddEntry(Histo_Map_total["h_total_radius"], "Total");
    l->Draw();

    t.DrawLatexNDC(0.5, 0.4, Form("Single %s", name.c_str()));
    t.DrawLatexNDC(0.15, 0.93935, "CLD #font[52]{work in progress}");
    t.DrawLatexNDC(0.58, 0.93935, "#scale[0.8]{Missing hits}");
    c_hits_radius->SaveAs("plots/Tracking_efficiency/hits_r_missing.pdf");

    TCanvas* c_hits_z = new TCanvas("c_hits_z", "comp", 200, 10, 700, 500);
    Histo_Maps["itec"]["h_itec_z"]->Draw("hist");
    Histo_Maps["itbc"]["h_itbc_z"]->Draw("histsame");
    Histo_Maps["otbc"]["h_otbc_z"]->Draw("histsame");
    Histo_Maps["otec"]["h_otec_z"]->Draw("histsame");
    Histo_Maps["vbc"]["h_vbc_z"]->Draw("histsame");
    Histo_Maps["vec"]["h_vec_z"]->Draw("histsame");
    Histo_Map_total["h_total_z"]->Draw("histsame");

    TLegend* l2 = new TLegend(0.5, 0.65, 0.75, 0.9);
    l2->AddEntry(Histo_Maps["itec"]["h_itec_z"], "Inner Tracker Endcap");
    l2->AddEntry(Histo_Maps["itbc"]["h_itbc_z"], "Inner Tracker Barrel");
    l2->AddEntry(Histo_Maps["otbc"]["h_otbc_z"], "Outer Tracker Barrel");
    l2->AddEntry(Histo_Maps["otec"]["h_otec_z"], "Outer Tracker Endcap");
    l2->AddEntry(Histo_Maps["vbc"]["h_vbc_z"], "Vertex Barrel");
    l2->AddEntry(Histo_Maps["vec"]["h_vec_z"], "Vertex Endcap");
    l2->AddEntry(Histo_Map_total["h_total_z"], "Total");
    l2->Draw();

    t.DrawLatexNDC(0.5, 0.4, Form("Single %s", name.c_str()));
    t.DrawLatexNDC(0.15, 0.93935, "CLD #font[52]{work in progress}");
    t.DrawLatexNDC(0.58, 0.93935, "#scale[0.8]{Missing hits}");
    c_hits_z->SaveAs("plots/Tracking_efficiency/hits_z_missing.pdf");

    TCanvas* c_hits_path = new TCanvas("c_hits_path", "comp", 200, 10, 700, 500);
    Histo_Maps["itec"]["h_itec_path"]->Draw("hist");
    Histo_Maps["itbc"]["h_itbc_path"]->Draw("histsame");
    Histo_Maps["otbc"]["h_otbc_path"]->Draw("histsame");
    Histo_Maps["otec"]["h_otec_path"]->Draw("histsame");
    Histo_Maps["vbc"]["h_vbc_path"]->Draw("histsame");
    Histo_Maps["vec"]["h_vec_path"]->Draw("histsame");
    Histo_Map_total["h_total_path"]->Draw("histsame");

    TLegend* l3 = new TLegend(0.5, 0.65, 0.75, 0.9);
    l3->AddEntry(Histo_Maps["itec"]["h_itec_path"], "Inner Tracker Endcap");
    l3->AddEntry(Histo_Maps["itbc"]["h_itbc_path"], "Inner Tracker Barrel");
    l3->AddEntry(Histo_Maps["otbc"]["h_otbc_path"], "Outer Tracker Barrel");
    l3->AddEntry(Histo_Maps["otec"]["h_otec_path"], "Outer Tracker Endcap");
    l3->AddEntry(Histo_Maps["vbc"]["h_vbc_path"], "Vertex Barrel");
    l3->AddEntry(Histo_Maps["vec"]["h_vec_path"], "Vertex Endcap");
    l3->AddEntry(Histo_Map_total["h_total_path"], "Total");
    l3->Draw();

    t.DrawLatexNDC(0.5, 0.4, Form("Single %s", name.c_str()));
    t.DrawLatexNDC(0.15, 0.93935, "CLD #font[52]{work in progress}");
    t.DrawLatexNDC(0.58, 0.93935, "#scale[0.8]{Missing hits}");
    c_hits_path->SaveAs("plots/Tracking_efficiency/hits_path_missing.pdf");

    TCanvas* c_hits_time = new TCanvas("c_hits_time", "comp", 200, 10, 700, 500);
    Histo_Maps["itec"]["h_itec_time"]->Draw("hist");
    Histo_Maps["itbc"]["h_itbc_time"]->Draw("histsame");
    Histo_Maps["otbc"]["h_otbc_time"]->Draw("histsame");
    Histo_Maps["otec"]["h_otec_time"]->Draw("histsame");
    Histo_Maps["vbc"]["h_vbc_time"]->Draw("histsame");
    Histo_Maps["vec"]["h_vec_time"]->Draw("histsame");
    Histo_Map_total["h_total_time"]->Draw("histsame");

    TLegend* l4 = new TLegend(0.5, 0.65, 0.75, 0.9);
    l4->AddEntry(Histo_Maps["itec"]["h_itec_time"], "Inner Tracker Endcap");
    l4->AddEntry(Histo_Maps["itbc"]["h_itbc_time"], "Inner Tracker Barrel");
    l4->AddEntry(Histo_Maps["otbc"]["h_otbc_time"], "Outer Tracker Barrel");
    l4->AddEntry(Histo_Maps["otec"]["h_otec_time"], "Outer Tracker Endcap");
    l4->AddEntry(Histo_Maps["vbc"]["h_vbc_time"], "Vertex Barrel");
    l4->AddEntry(Histo_Maps["vec"]["h_vec_time"], "Vertex Endcap");
    l4->AddEntry(Histo_Map_total["h_total_time"], "Total");
    l4->Draw();

    t.DrawLatexNDC(0.5, 0.4, Form("Single %s", name.c_str()));
    t.DrawLatexNDC(0.15, 0.93935, "CLD #font[52]{work in progress}");
    t.DrawLatexNDC(0.58, 0.93935, "#scale[0.8]{Missing hits}");

    c_hits_time->SaveAs("plots/Tracking_efficiency/hits_time_missing.pdf");

    TCanvas* c_hits_eDep = new TCanvas("c_hits_eDep", "comp", 200, 10, 700, 500);
    Histo_Maps["itec"]["h_itec_eDep"]->Draw("hist");
    Histo_Maps["itbc"]["h_itbc_eDep"]->Draw("histsame");
    Histo_Maps["otbc"]["h_otbc_eDep"]->Draw("histsame");
    Histo_Maps["otec"]["h_otec_eDep"]->Draw("histsame");
    Histo_Maps["vbc"]["h_vbc_eDep"]->Draw("histsame");
    Histo_Maps["vec"]["h_vec_eDep"]->Draw("histsame");
    Histo_Map_total["h_total_eDep"]->Draw("histsame");

    TLegend* l5 = new TLegend(0.5, 0.65, 0.75, 0.9);
    l5->AddEntry(Histo_Maps["itec"]["h_itec_eDep"], "Inner Tracker Endcap");
    l5->AddEntry(Histo_Maps["itbc"]["h_itbc_eDep"], "Inner Tracker Barrel");
    l5->AddEntry(Histo_Maps["otbc"]["h_otbc_eDep"], "Outer Tracker Barrel");
    l5->AddEntry(Histo_Maps["otec"]["h_otec_eDep"], "Outer Tracker Endcap");
    l5->AddEntry(Histo_Maps["vbc"]["h_vbc_eDep"], "Vertex Barrel");
    l5->AddEntry(Histo_Maps["vec"]["h_vec_eDep"], "Vertex Endcap");
    l5->AddEntry(Histo_Map_total["h_total_eDep"], "Total");
    l5->Draw();

    t.DrawLatexNDC(0.5, 0.4, Form("Single %s", name.c_str()));
    t.DrawLatexNDC(0.15, 0.93935, "CLD #font[52]{work in progress}");
    t.DrawLatexNDC(0.58, 0.93935, "#scale[0.8]{Missing hits}");

    c_hits_eDep->SaveAs("plots/Tracking_efficiency/hits_eDep_missing.pdf");

    TCanvas* c_hits_pt = new TCanvas("c_hits_pt", "comp", 200, 10, 700, 500);
    TPad *p1 = new TPad("p1", "", 0, 0, 1, 1);
    TPad *p2 = new TPad("p2", "", 0, 0, 1, 1);
    p1->SetFillStyle(4000); // will be transparent
    p2->SetFillStyle(4000); 

    p1->Draw();
    p1->cd(); 

    Histo_Maps["itec"]["h_itec_pt"]->Draw("hist");
    Histo_Maps["itbc"]["h_itbc_pt"]->Draw("histsame");
    Histo_Maps["otbc"]["h_otbc_pt"]->Draw("histsame");
    Histo_Maps["otec"]["h_otec_pt"]->Draw("histsame");
    Histo_Maps["vbc"]["h_vbc_pt"]->Draw("histsame");
    Histo_Maps["vec"]["h_vec_pt"]->Draw("histsame");
    Histo_Map_total["h_total_pt"]->Draw("histsame");

    // Create the efficiency graph
    Double_t xmin = 0;
    Double_t xmax = 100;
    Double_t dx = (xmax - xmin) / 0.8; // 10 percent margins left and right
    Double_t ymin = 0.;
    Double_t ymax = 1.05;
    Double_t dy = (ymax - ymin) / 0.8; // 10 percent margins top and bottom
    p2->Range(xmin-0.1*dx, ymin-0.1*dy, xmax+0.1*dx, ymax+0.1*dy);
    p2->Draw();
    p2->cd();
    eff_pt->SetLineColor(kViolet-7);
    eff_pt->SetMarkerColor(kViolet-7);
    eff_pt->SetMarkerSize(2);
    eff_pt->SetLineWidth(2);
    eff_pt->Draw("psame");

    // Draw an axis on the right side
    TGaxis *axis = new TGaxis(xmax, ymin, xmax, ymax, ymin, ymax, 510, "+L");
    axis->SetLineColor(kViolet-7);
    axis->SetLabelColor(kViolet-7);
    axis->SetTitleColor(kViolet-7);
    axis->SetTitle("Efficiency");
    axis->Draw();


    TLegend* l6 = new TLegend(0.5, 0.55, 0.75, 0.8);
    l6->AddEntry(Histo_Maps["itec"]["h_itec_pt"], "Inner Tracker Endcap");
    l6->AddEntry(Histo_Maps["itbc"]["h_itbc_pt"], "Inner Tracker Barrel");
    l6->AddEntry(Histo_Maps["otbc"]["h_otbc_pt"], "Outer Tracker Barrel");
    l6->AddEntry(Histo_Maps["otec"]["h_otec_pt"], "Outer Tracker Endcap");
    l6->AddEntry(Histo_Maps["vbc"]["h_vbc_pt"], "Vertex Barrel");
    l6->AddEntry(Histo_Maps["vec"]["h_vec_pt"], "Vertex Endcap");
    l6->AddEntry(Histo_Map_total["h_total_pt"], "Total");
    l6->Draw();

    t.DrawLatexNDC(0.5, 0.4, Form("Single %s", name.c_str()));
    t.DrawLatexNDC(0.15, 0.93935, "CLD #font[52]{work in progress}");
    t.DrawLatexNDC(0.58, 0.93935, "#scale[0.8]{Missing hits}");

    c_hits_pt->SaveAs("plots/Tracking_efficiency/hits_pt_missing_eff.pdf");

    TCanvas* c_hits_p = new TCanvas("c_hits_p", "comp", 200, 10, 700, 500);
    TPad *p3 = new TPad("p3", "", 0, 0, 1, 1);
    TPad *p4 = new TPad("p4", "", 0, 0, 1, 1);
    p3->SetFillStyle(4000); // will be transparent
    p4->SetFillStyle(4000); 

    p3->Draw();
    p3->cd();

    Histo_Maps["itec"]["h_itec_p"]->Draw("hist");
    Histo_Maps["itbc"]["h_itbc_p"]->Draw("histsame");
    Histo_Maps["otbc"]["h_otbc_p"]->Draw("histsame");
    Histo_Maps["otec"]["h_otec_p"]->Draw("histsame");
    Histo_Maps["vbc"]["h_vbc_p"]->Draw("histsame");
    Histo_Maps["vec"]["h_vec_p"]->Draw("histsame");
    Histo_Map_total["h_total_p"]->Draw("histsame");

    // Create the second graph
    p4->Range(xmin-0.1*dx, ymin-0.1*dy, xmax+0.1*dx, ymax+0.1*dy);
    p4->Draw();
    p4->cd();
    eff_p->SetLineColor(kViolet-7);
    eff_p->SetMarkerColor(kViolet-7);
    eff_p->SetMarkerSize(2);
    eff_p->SetLineWidth(2);
    eff_p->Draw("psame");

    // Draw an axis on the right side
    TGaxis *axis2 = new TGaxis(xmax, ymin, xmax, ymax, ymin, ymax, 510, "+L");
    axis2->SetLineColor(kViolet-7);
    axis2->SetLabelColor(kViolet-7);
    axis2->SetTitleColor(kViolet-7);
    axis2->SetTitle("Efficiency");
    axis2->Draw();

    TLegend* l7 = new TLegend(0.5, 0.55, 0.75, 0.8);
    l7->AddEntry(Histo_Maps["itec"]["h_itec_p"], "Inner Tracker Endcap");
    l7->AddEntry(Histo_Maps["itbc"]["h_itbc_p"], "Inner Tracker Barrel");
    l7->AddEntry(Histo_Maps["otbc"]["h_otbc_p"], "Outer Tracker Barrel");
    l7->AddEntry(Histo_Maps["otec"]["h_otec_p"], "Outer Tracker Endcap");
    l7->AddEntry(Histo_Maps["vbc"]["h_vbc_p"], "Vertex Barrel");
    l7->AddEntry(Histo_Maps["vec"]["h_vec_p"], "Vertex Endcap");
    l7->AddEntry(Histo_Map_total["h_total_p"], "Total");
    l7->Draw();

    t.DrawLatexNDC(0.5, 0.4, Form("Single %s", name.c_str()));
    t.DrawLatexNDC(0.15, 0.93935, "CLD #font[52]{work in progress}");
    t.DrawLatexNDC(0.58, 0.93935, "#scale[0.8]{Missing hits}");

    c_hits_p->SaveAs("plots/Tracking_efficiency/hits_p_missing_eff.pdf");

    outputfile->Write();
    outputfile->Close(); 
}

