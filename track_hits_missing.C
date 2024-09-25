
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
#include "podio/RelationRange.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/MCRecoTrackParticleAssociation.h"
#include "edm4hep/TrackerHitSimTrackerHitLink.h"
#include "edm4hep/ReconstructedParticleCollection.h"
#include "edm4hep/MCRecoParticleAssociationCollection.h"
#include "edm4hep/TrackCollection.h"
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

void track_hits_missing(const char* input, const char* output, int pdg_number) {
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
    std::map<std::string,TH1D*> Histo_Map_si;
    Histo_Map_si["h_itec_hits"] = new TH1D("h_itec_hits", ";#theta [mrad]", 30, 0, TMath::Pi()*1000);
    Histo_Map_si["h_itbc_hits"] = new TH1D("h_itbc_hits", ";#theta [mrad]", 30, 0, TMath::Pi()*1000);
    Histo_Map_si["h_otbc_hits"] = new TH1D("h_otbc_hits", ";#theta [mrad]", 30, 0, TMath::Pi()*1000);
    Histo_Map_si["h_otec_hits"] = new TH1D("h_otec_hits", ";#theta [mrad]", 30, 0, TMath::Pi()*1000);
    Histo_Map_si["h_vbc_hits"] = new TH1D("h_vbc_hits", ";#theta [mrad]", 30, 0, TMath::Pi()*1000);
    Histo_Map_si["h_vec_hits"] = new TH1D("h_vec_hits", ";#theta [mrad]", 30, 0, TMath::Pi()*1000);
    Histo_Map_si["h_total_hits"] = new TH1D("h_total_hits", ";#theta [mrad]", 30, 0, TMath::Pi()*1000);

    TH1D* h_total_theta = new TH1D("h_total_theta", ";#theta [mrad]", 30, 0, TMath::Pi()*1000);

    int all = 0;

    for (size_t i = 0; i < reader.getEntries("events"); ++i) {
        if (all >= 10000) break;
        all++;

        const auto event = podio::Frame(reader.readNextEntry("events"));
        const auto &mcparticles = event.get<edm4hep::MCParticleCollection>("MCParticles");
        const auto &mctruthrecolinks = event.get<edm4hep::TrackMCParticleLinkCollection>("SiTracksMCTruthLink"); //before MCRecoTrackParticleAssociation
        const auto& mclinks = event.get<edm4hep::MCRecoParticleAssociationCollection>("MCTruthRecoLink");

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
        std::vector<std::tuple<std::string,const edm4hep::SimTrackerHitCollection*,const edm4hep::TrackerHitSimTrackerHitLinkCollection*>> SimTrackerHit = {
            {"itec", &itec, &itec_link},
            {"itbc", &itbc, &itbc_link},
            {"otbc", &otbc, &otbc_link},
            {"otec", &otec, &otec_link},
            {"vbc", &vbc, &vbc_link},
            {"vec", &vec, &vec_link}
        };
        //std::cout << i << std::endl;
        for (const auto& mcp : mcparticles) {
            MC_Particle particle;
            particle = getStableMCParticle(mcp, particle, pdg_number);
            if (particle.flag == true) { 
                h_total_theta->Fill(particle.theta *1000);
                for (const auto& simtracker : SimTrackerHit) {
                    std::string name_sim;
                    const edm4hep::SimTrackerHitCollection*  hit_collection;
                    const edm4hep::TrackerHitSimTrackerHitLinkCollection* hit_link;    
                    std::tie(name_sim, hit_collection,  hit_link) = simtracker;    
                    for (const auto& hit : *hit_collection) {
                        if (hit.getParticle() == mcp) { 
                            Tracker_Hit tracker_hit;
                            tracker_hit = getHit(mcp, hit, tracker_hit);
                            std::string histo_si_name = "h_" + name_sim + "_hits";
                            Reconstructed_Particle reco_p;
                            reco_p = getReco(mcp, mclinks, reco_p); 
                            auto tracks = reco_p.reco.first.getTracks(); //tracks that have been used for this particle
                            bool flag;
                            flag = getNumberofHitReco(tracks, hit, *hit_link); 
                            bool flag_secondary = hit.isProducedBySecondary();
                            if (flag == false && flag_secondary == false){
                                Histo_Map_si[histo_si_name]->Fill(particle.theta *1000);
                                Histo_Map_si["h_total_hits"]->Fill(particle.theta *1000);                        
                            }                                                                                                                
                        }
                    }
                } 
            }
        } 
    } 
          
    // Average amount of hits
    for (auto& histo: Histo_Map_si) {   
        histo.second->Divide(h_total_theta);
    } 

    // Create canvas and draw histograms
    TLatex t;
    t.SetTextSize(30);
    t.SetTextFont(63);
    gStyle->Reset("Modern");
    TCanvas* c_hits = new TCanvas("c", "comp", 200, 10, 700, 500);
     
    // Create the first graph
    for (auto& histo: Histo_Map_si) {   
        histo.second->SetLineWidth(2);
        histo.second->GetYaxis()->SetRangeUser(0,6);
    }
    
    Histo_Map_si["h_itec_hits"]->GetXaxis()->SetRangeUser(0, TMath::Pi()*1000);
    Histo_Map_si["h_itec_hits"]->GetYaxis()->SetTitle(Form("Number of hits per %s missing in the track", name.c_str()));
    Histo_Map_si["h_itec_hits"]->GetYaxis()->SetTitleSize(0.05);
    Histo_Map_si["h_itec_hits"]->GetXaxis()->SetTitleSize(0.05);
    Histo_Map_si["h_itec_hits"]->GetXaxis()->SetTitleOffset(0.8);
    Histo_Map_si["h_itec_hits"]->GetYaxis()->SetTitleOffset(0.8);
    Histo_Map_si["h_itec_hits"]->SetLineColor(kCyan+1);
    Histo_Map_si["h_itbc_hits"]->SetLineColor(kGreen+2);
    Histo_Map_si["h_otbc_hits"]->SetLineColor(kRed+2);
    Histo_Map_si["h_otec_hits"]->SetLineColor(kYellow+1);
    Histo_Map_si["h_vbc_hits"]->SetLineColor(kRed-7);
    Histo_Map_si["h_vec_hits"]->SetLineColor(kBlue+1);
    Histo_Map_si["h_total_hits"]->SetLineColor(kBlack);
    Histo_Map_si["h_itec_hits"]->Draw("hist");
    Histo_Map_si["h_itbc_hits"]->Draw("histsame");
    Histo_Map_si["h_otbc_hits"]->Draw("histsame");
    Histo_Map_si["h_otec_hits"]->Draw("histsame");
    Histo_Map_si["h_vbc_hits"]->Draw("histsame");  
    Histo_Map_si["h_vec_hits"]->Draw("histsame");
    Histo_Map_si["h_total_hits"]->Draw("histsame");

    TLegend* l = new TLegend(0.1, 0.6, 0.35, 0.9);
    l->AddEntry(Histo_Map_si["h_itec_hits"], "Inner Tracker Endcap");
    l->AddEntry(Histo_Map_si["h_itbc_hits"], "Inner Tracker Barrel");
    l->AddEntry(Histo_Map_si["h_otbc_hits"], "Outer Tracker Barrel");
    l->AddEntry(Histo_Map_si["h_otec_hits"], "Outer Tracker Endcap");
    l->AddEntry(Histo_Map_si["h_vbc_hits"], "Vertex Barrel");
    l->AddEntry(Histo_Map_si["h_vec_hits"], "Vertex Endcap");
    l->AddEntry(Histo_Map_si["h_total_hits"], "Total");
    l->Draw(); 

    t.DrawLatexNDC(0.4, 0.55, Form("Single %s", name.c_str()));
    t.DrawLatexNDC(0.15, 0.93935, "CLD #font[52]{work in progress}");
    //t.DrawLatexNDC(0.68, 0.93935, "#scale[0.8]{p_{T} > 10GeV}");

    c_hits->SaveAs("plots/Missing_hits/hits_missing.pdf");

    outputfile->Write();
    outputfile->Close(); 
}

