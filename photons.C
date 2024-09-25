
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
#include "edm4hep/CaloHitMCParticleLinkCollection.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/CaloHitContribution.h"
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

void photons(const char* input, const char* output, int pdg_number) {
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
    
    // Histograms
    TH1D* h_found_radius = new TH1D("h_found_radius", ";#theta [mrad]", 30, 0, TMath::Pi()*1000);
    TH1D* h_daughter = new TH1D("h_daughter", ";#theta [mrad]", 30, 0, TMath::Pi()*1000); //ph rec as e
    TH1D* h_photons = new TH1D("h_photons", ";#theta [mrad]", 30, 0, TMath::Pi()*1000); // ph rec as photons
    TH1D* h_others = new TH1D("h_others", ";#theta [mrad]", 30, 0, TMath::Pi()*1000); //ph rec as other particles

    std::map<std::string,TH1D*> Histo_Map_si;
    Histo_Map_si["h_ecb_hits"] = new TH1D("h_ecb_hits", ";#theta [mrad]", 30, 0, TMath::Pi()*1000);
    Histo_Map_si["h_ece_hits"] = new TH1D("h_ece_hits", ";#theta [mrad]", 30, 0, TMath::Pi()*1000);
    Histo_Map_si["h_hcb_hits"] = new TH1D("h_hcb_hits", ";#theta [mrad]", 30, 0, TMath::Pi()*1000);
    Histo_Map_si["h_hce_hits"] = new TH1D("h_hce_hits", ";#theta [mrad]", 30, 0, TMath::Pi()*1000);
    Histo_Map_si["h_hcr_hits"] = new TH1D("h_hcr_hits", ";#theta [mrad]", 30, 0, TMath::Pi()*1000);
    Histo_Map_si["h_lumical_hits"] = new TH1D("h_lumical_hits", ";#theta [mrad]", 30, 0, TMath::Pi()*1000);
    Histo_Map_si["h_yb_hits"] = new TH1D("h_yb_hits", ";#theta [mrad]", 30, 0, TMath::Pi()*1000);
    Histo_Map_si["h_ye_hits"] = new TH1D("h_ye_hits", ";#theta [mrad]", 30, 0, TMath::Pi()*1000);

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

        //Calorimiter
        const auto& calolinks = event.get<edm4hep::CaloHitMCParticleLinkCollection>("CalohitMCTruthLink"); 
        const auto& calohit_cont = event.get<edm4hep::CaloHitContributionCollection>("ToolSvc.lcio2EDM4hep_CaloHitContributions");
        const auto& ecb = event.get<edm4hep::SimCalorimeterHitCollection>("ECalBarrelCollection");
        const auto& ece = event.get<edm4hep::SimCalorimeterHitCollection>("ECalEndcapCollection");
        const auto& hcb = event.get<edm4hep::SimCalorimeterHitCollection>("HCalBarrelCollection");
        const auto& hce = event.get<edm4hep::SimCalorimeterHitCollection>("HCalEndcapCollection");
        const auto& hcr = event.get<edm4hep::SimCalorimeterHitCollection>("HCalRingCollection");
        const auto& lumical = event.get<edm4hep::SimCalorimeterHitCollection>("LumiCalCollection ");
        const auto& yb = event.get<edm4hep::SimCalorimeterHitCollection>("YokeBarrelCollection");
        const auto& ye = event.get<edm4hep::SimCalorimeterHitCollection>("YokeEndcapCollection");


        // Map for the SimTrackerHit
        std::vector<std::tuple<std::string,const edm4hep::SimTrackerHitCollection*,const edm4hep::TrackerHitSimTrackerHitLinkCollection*>> SimTrackerHit = {
            {"itec", &itec, &itec_link},
            {"itbc", &itbc, &itbc_link},
            {"otbc", &otbc, &otbc_link},
            {"otec", &otec, &otec_link},
            {"vbc", &vbc, &vbc_link},
            {"vec", &vec, &vec_link}
        };

        // Map for the CalorimiterHit
        std::vector<std::tuple<std::string,const edm4hep::SimCalorimeterHitCollection*>> SimCalorimiterHit = {
            {"ecb", &ecb},
            {"ece", &ece},
            {"hcb", &hcb},
            {"hce", &hce},
            {"hcr", &hcr},
            {"lumical", &lumical},
            {"yb", &yb},
            {"ye", &ye},
        };

        for (const auto& mcp : mcparticles) {
            MC_Particle particle;
            particle = getStableMCParticle(mcp, particle, pdg_number);
            if (particle.flag == true ) { 
                minimum min;
                maximum max;
                min.radius_min = 10000.; 
                max.radius_max = 0.; 
                bool is_missing = false; 
                for (const auto& simtracker : SimTrackerHit) {
                    std::string name_sim;
                    const edm4hep::SimTrackerHitCollection*  hit_collection;
                    const edm4hep::TrackerHitSimTrackerHitLinkCollection* hit_link;    
                    std::tie(name_sim, hit_collection,  hit_link) = simtracker;    
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
                            if (flag == false && flag_secondary == false){ // primary missing  
                                is_missing = true;  
                                min = getMinRadius(tracker_hit, min);                                             
                            } 
                            if (flag == true && flag_secondary == false){ // primary not missing  
                                max = getMaxRadius(tracker_hit, max, hit, mcp);                                                         
                            } 
                        }
                    }
                } 
                if(is_missing == true) {
                    if( max.radius_max < min.radius_min && abs(max.z_max) < abs(min.z_min) ) { 
                        MC_Particle daughter =  getDaughters(mcp, max, min);
                        float total_energy = getTotalEnergy(mcp, max, min);
                        if (daughter.flag == true && abs(daughter.pdg) == 22) {
                            for (const auto& link : mclinks) {
                                if (link.getSim() == daughter.mc_particle ) {
                                    h_found_radius->Fill(particle.theta*1000);
                                    if (abs(link.getRec().getPDG()) == 11) {
                                        h_daughter->Fill(particle.theta*1000);
                                    }
                                    else if (abs(link.getRec().getPDG()) == 22) {
                                        h_photons->Fill(particle.theta*1000);
                                    }
                                    else h_others->Fill(particle.theta*1000);
                                }
                            }
                            for (const auto& simcalorimiter : SimCalorimiterHit) {
                                std::string name_calo;
                                const edm4hep::SimCalorimeterHitCollection*  calo_hit_collection;   
                                std::tie(name_calo, calo_hit_collection) = simcalorimiter;    
                                for (const auto& hit : *calo_hit_collection) {
                                    auto contributions = hit.getContributions();
                                    for (const auto& contribution : contributions) {
                                        if (contribution.getParticle() == daughter.mc_particle) {
                                            std::string histo_si_name = "h_" + name_calo + "_hits";
                                            Histo_Map_si[histo_si_name]->Fill(particle.theta *1000);
                                        }
                                    }
                                }
                            }                                  
                        }
                    } 
                }
            }                               
        } 
    } 

    // Percentage of right missing
    h_daughter->Divide(h_found_radius);
    h_photons->Divide(h_found_radius);
    h_others->Divide(h_found_radius);
    for (auto& histo: Histo_Map_si) {   
        histo.second->Divide(h_found_radius);
    } 
         
    // Create canvas and draw histograms
    TLatex t;
    t.SetTextSize(30);
    t.SetTextFont(63);
    gStyle->Reset("Modern");
    TCanvas* c_reco = new TCanvas("c", "comp", 200, 10, 700, 500);
    
    // Create the first graph   
    h_daughter->SetLineWidth(2);
    h_daughter->GetYaxis()->SetRangeUser(0, 1.05);
    h_photons->SetLineWidth(2);
    h_photons->GetYaxis()->SetRangeUser(0, 1.05);
    h_others->SetLineWidth(2);
    h_others->GetYaxis()->SetRangeUser(0, 1.05);

    h_daughter->GetXaxis()->SetRangeUser(0, TMath::Pi()*1000);
    h_daughter->GetYaxis()->SetTitle(Form("Secondary #gamma [%]", name.c_str())); 
    h_daughter->GetYaxis()->SetTitleSize(0.05);
    h_daughter->GetXaxis()->SetTitleSize(0.05);
    h_daughter->GetXaxis()->SetTitleOffset(0.8);
    h_daughter->GetYaxis()->SetTitleOffset(0.8);
    h_daughter->SetLineColor(kRed);
    h_photons->SetLineColor(kYellow+1);
    h_others->SetLineColor(kBlack);
    h_daughter->Draw("hist");
    h_photons->Draw("histsame"); 
    h_others->Draw("histsame");

    TLegend* l = new TLegend(0.2, 0.7, 0.8, 0.9);
    l->AddEntry(h_daughter, "Reconstructed as e^{#pm}");
    l->AddEntry(h_photons, "Reconstructed as #gamma");
    l->AddEntry(h_others, "Reconstructed as other particles");
    l->Draw(); 

    t.DrawLatexNDC(0.45, 0.3, Form("Single %s", name.c_str()));
    t.DrawLatexNDC(0.15, 0.93935, "CLD #font[52]{work in progress}");

    c_reco->SaveAs("plots/Missing_hits/Reco_photons.pdf");

    // Create the second graph
    TCanvas* c_calo = new TCanvas("c2", "comp", 200, 10, 700, 500);
    for (auto& histo: Histo_Map_si) {   
        histo.second->SetLineWidth(2);
        histo.second->GetYaxis()->SetRangeUser(0,20000);
    }
    
    Histo_Map_si["h_ecb_hits"]->GetXaxis()->SetRangeUser(0, TMath::Pi()*1000);
    Histo_Map_si["h_ecb_hits"]->GetYaxis()->SetTitle("Number of hits per secondary #gamma");
    Histo_Map_si["h_ecb_hits"]->GetYaxis()->SetTitleSize(0.05);
    Histo_Map_si["h_ecb_hits"]->GetXaxis()->SetTitleSize(0.05);
    Histo_Map_si["h_ecb_hits"]->GetXaxis()->SetTitleOffset(0.8);
    Histo_Map_si["h_ecb_hits"]->GetYaxis()->SetTitleOffset(0.8);
    Histo_Map_si["h_ecb_hits"]->SetLineColor(kCyan+1);
    Histo_Map_si["h_ece_hits"]->SetLineColor(kGreen+2);
    Histo_Map_si["h_hcb_hits"]->SetLineColor(kRed+2);
    Histo_Map_si["h_hce_hits"]->SetLineColor(kYellow+1);
    Histo_Map_si["h_hcr_hits"]->SetLineColor(kRed-7);
    Histo_Map_si["h_lumical_hits"]->SetLineColor(kBlue+1);
    Histo_Map_si["h_yb_hits"]->SetLineColor(kBlack);
    Histo_Map_si["h_ye_hits"]->SetLineColor(kYellow-7);
    Histo_Map_si["h_ecb_hits"]->Draw("hist");
    Histo_Map_si["h_ece_hits"]->Draw("histsame");
    Histo_Map_si["h_hcb_hits"]->Draw("histsame");
    Histo_Map_si["h_hce_hits"]->Draw("histsame");
    Histo_Map_si["h_hcr_hits"]->Draw("histsame");  
    Histo_Map_si["h_lumical_hits"]->Draw("histsame");
    Histo_Map_si["h_yb_hits"]->Draw("histsame");
    Histo_Map_si["h_ye_hits"]->Draw("histsame");

    TLegend* l = new TLegend(0.2, 0.6, 0.45, 0.9);
    l->AddEntry(Histo_Map_si["h_ecb_hits"], "ECAL Barrel");
    l->AddEntry(Histo_Map_si["h_ece_hits"], "ECAL Endcup");
    l->AddEntry(Histo_Map_si["h_hcb_hits"], "HCAL Barrel");
    l->AddEntry(Histo_Map_si["h_hce_hits"], "HCAL Endcup");
    l->AddEntry(Histo_Map_si["h_hcr_hits"], "HCAL Ring");
    l->AddEntry(Histo_Map_si["h_lumical_hits"], "Lumical");
    l->AddEntry(Histo_Map_si["h_yb_hits"], "Yoke Barrel");
    l->AddEntry(Histo_Map_si["h_ye_hits"], "yoke Endcup");
    l->Draw(); 


    t.DrawLatexNDC(0.45, 0.3, Form("Single %s", name.c_str()));
    t.DrawLatexNDC(0.15, 0.93935, "CLD #font[52]{work in progress}");

    c_calo->SaveAs("plots/Missing_hits/Calorimiter.pdf");


    outputfile->Write();
    outputfile->Close(); 
}

