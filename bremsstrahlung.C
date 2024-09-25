
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

void bremsstrahlung(const char* input, const char* output, int pdg_number) {
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
    TH1D* h_missing_radius = new TH1D("h_missing_radius", ";#theta [mrad]", 30, 0, TMath::Pi()*1000);
    TH1D* h_found_radius = new TH1D("h_found_radius", ";#theta [mrad]", 30, 0, TMath::Pi()*1000);
    TH1D* h_daughter_radius = new TH1D("h_daughter_radius", ";#theta [mrad]", 30, 0, TMath::Pi()*1000);
    TH1D* h_photon_radius = new TH1D("h_photon_radius", ";#theta [mrad]", 30, 0, TMath::Pi()*1000);
    TH1D* h_electron_radius = new TH1D("h_electron_radius", ";#theta [mrad]", 30, 0, TMath::Pi()*1000);
    TH1D* h_other_radius = new TH1D("h_other_radius", ";#theta [mrad]", 30, 0, TMath::Pi()*1000);

     Int_t nbins = 30;

     Double_t xmin = 0.001;  // 10^(-3)
     Double_t xmax = 100; // 10^3

    Double_t bins[nbins+1];
    for (int i = 0; i <= nbins; ++i) {
    bins[i] = TMath::Power(10, TMath::Log10(xmin) + i*(TMath::Log10(xmax) - TMath::Log10(xmin))/nbins);
    }

    TH1D* h_daughter_energy = new TH1D("h_daughter_energy", ";Energy sec. particle [GeV]", nbins, bins);
    TH1D* h_photon_energy = new TH1D("h_photon_energy", ";Energy sec. particle [GeV]", nbins, bins);
    TH1D* h_electron_energy = new TH1D("h_electron_energy", ";Energy sec. particle [GeV]", nbins, bins);
    TH1D* h_other_energy = new TH1D("h_other_energy", ";Energy sec. particle [GeV]", nbins, bins); 

    

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
                    h_found_radius->Fill(particle.theta*1000);
                    if( max.radius_max < min.radius_min && abs(max.z_max) < abs(min.z_min) ) { 
                        h_missing_radius->Fill(particle.theta*1000);
                        MC_Particle daughter =  getDaughters(mcp, max, min);
                        if (daughter.flag == true) {
                            h_daughter_radius->Fill(particle.theta*1000);
                            h_daughter_energy->Fill(daughter.energy);
                            if(abs(daughter.pdg) == 22) {
                                 h_photon_radius->Fill(particle.theta*1000);
                                 h_photon_energy->Fill(daughter.energy);
                            }
                            else if(abs(daughter.pdg) == 11) {
                                h_electron_radius->Fill(particle.theta*1000);
                                h_electron_energy->Fill(daughter.energy);

                            }
                            else {
                                h_other_radius->Fill(particle.theta*1000);
                                h_other_energy->Fill(daughter.energy);  
                            }                        
                        }
                    } 
                }
            }                               
        } 
    } 

    // Percentage of right missing
    h_missing_radius->Divide(h_found_radius);
    h_daughter_radius->Divide(h_found_radius);
    h_photon_radius->Divide(h_found_radius);
    h_electron_radius->Divide(h_found_radius);
    h_other_radius->Divide(h_found_radius);
        
    // Create canvas and draw histograms
    TLatex t;
    t.SetTextSize(30);
    t.SetTextFont(63);
    gStyle->Reset("Modern");
    TCanvas* c_hits = new TCanvas("c", "comp", 200, 10, 700, 500);
     
    // Create the first graph   
    h_missing_radius->SetLineWidth(2);
    h_missing_radius->GetYaxis()->SetRangeUser(0, 1.05);
    h_daughter_radius->SetLineWidth(2);
    h_daughter_radius->GetYaxis()->SetRangeUser(0, 1.05);

    h_photon_radius->SetLineWidth(2);
    h_photon_radius->GetYaxis()->SetRangeUser(0, 1.05);
    h_electron_radius->SetLineWidth(2);
    h_electron_radius->GetYaxis()->SetRangeUser(0, 1.05);
    h_other_radius->SetLineWidth(2);
    h_other_radius->GetYaxis()->SetRangeUser(0, 1.05);
 
    h_missing_radius->GetXaxis()->SetRangeUser(0, TMath::Pi()*1000);
    h_missing_radius->GetYaxis()->SetTitle(Form("Percentage of primary %s", name.c_str())); 
    h_missing_radius->GetYaxis()->SetTitleSize(0.05);
    h_missing_radius->GetXaxis()->SetTitleSize(0.05);
    h_missing_radius->GetXaxis()->SetTitleOffset(0.8);
    h_missing_radius->GetYaxis()->SetTitleOffset(0.8);

    h_missing_radius->SetLineColor(kBlue+1);
    h_missing_radius->Draw("hist");
    h_daughter_radius->SetLineColor(kRed);
    h_daughter_radius->Draw("histsame");
    h_photon_radius->SetLineColor(kYellow +1);
    h_photon_radius->Draw("histsame");
    h_electron_radius->SetLineColor(kGreen+2);
    h_electron_radius->Draw("histsame");
    h_other_radius->SetLineColor(kBlack);
    h_other_radius->Draw("histsame");

    TLegend* l = new TLegend(0.2, 0.65, 0.8, 0.9);
    l->AddEntry(h_missing_radius, "radius (and z) max in the track < radius (and z) min missing");
    l->AddEntry(h_daughter_radius, "Vertex of one daughter between radius (and z) max and min");
    l->AddEntry(h_photon_radius, "Photons");
    l->AddEntry(h_electron_radius, "Electrons/positrons");
    l->AddEntry(h_other_radius, "Other particles");
    l->Draw();

    t.DrawLatexNDC(0.4, 0.55, Form("Single %s", name.c_str()));
    t.DrawLatexNDC(0.15, 0.93935, "CLD #font[52]{work in progress}");

    c_hits->SaveAs("plots/Missing_hits/Bremsstrahlung.pdf");


    TCanvas* c_energy = new TCanvas("c2", "comp", 800,600);
     
    // Create the second graph   
    gPad->SetLogx();
    h_daughter_energy->SetLineWidth(2);
    h_daughter_energy->GetYaxis()->SetRangeUser(0, 100);
    h_photon_energy->SetLineWidth(2);
    h_photon_energy->GetYaxis()->SetRangeUser(0, 100);
    h_electron_energy->SetLineWidth(2);
    h_electron_energy->GetYaxis()->SetRangeUser(0, 100);
    h_other_energy->SetLineWidth(2);
    h_other_energy->GetYaxis()->SetRangeUser(0, 100);
    h_daughter_energy->GetXaxis()->SetRangeUser(0, 100);
    h_daughter_energy->GetYaxis()->SetTitle("Number of secondary particles");
    h_daughter_energy->GetYaxis()->SetTitleSize(0.05);
    h_daughter_energy->GetXaxis()->SetTitleSize(0.05);
    h_daughter_energy->GetXaxis()->SetTitleOffset(0.8);
    h_daughter_energy->GetYaxis()->SetTitleOffset(0.8);

    h_daughter_energy->SetLineColor(kRed);
    h_daughter_energy->Draw("hist");
    h_photon_energy->SetLineColor(kYellow +1);
    h_photon_energy->Draw("histsame");
    h_electron_energy->SetLineColor(kGreen+2);
    h_electron_energy->Draw("histsame");
    h_other_energy->SetLineColor(kBlack);
    h_other_energy->Draw("histsame");


    TLegend* l2 = new TLegend(0.2, 0.65, 0.8, 0.9);
    l2->AddEntry(h_daughter_energy, "Total daughters");
    l2->AddEntry(h_photon_energy, "Photons");
    l2->AddEntry(h_electron_energy, "Electrons/positrons");
    l2->AddEntry(h_other_energy, "Other particles");
    l2->Draw();
    t.DrawLatexNDC(0.5, 0.55, Form("Single %s", name.c_str()));
    t.DrawLatexNDC(0.15, 0.93935, "CLD #font[52]{work in progress}");

    c_energy->SaveAs("plots/Missing_hits/Energy_daughters.pdf");

    outputfile->Write();
    outputfile->Close(); 
}

