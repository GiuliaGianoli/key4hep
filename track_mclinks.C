
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

// Struct to hold particle data
struct MC_Particle {
    float theta;
    float pt;
    float p;  
    bool flag; 
};

// Struct for reconstructed particle
struct Reconstructed_Particle {
    float trwgt;
    edm4hep::ReconstructedParticle reco;
    bool flag;
};

// Function to get electron data from MCParticleCollection
MC_Particle getElectron(const edm4hep::MCParticle& mcp, MC_Particle particle, int pdg_number) {  
    if (abs(mcp.getPDG()) == pdg_number && mcp.getGeneratorStatus() == 1) {  
        auto mom = mcp.getMomentum(); 
        TVector3 mom_v(mom.x, mom.y, mom.z); 
        particle.theta = mom_v.Theta();
        particle.pt = mom_v.Pt();
        particle.p = sqrt(mom_v.Px()*mom_v.Px()+mom_v.Py()*mom_v.Py()+mom_v.Pz()*mom_v.Pz());
        particle.flag = true;
    }
    else particle.flag = false;
    return particle;
}
    
// Function to get reconstructed particle
Reconstructed_Particle getReco(const edm4hep::MCParticle& mcp, const edm4hep::MCRecoParticleAssociation& link, Reconstructed_Particle reco_p) {       
    reco_p.trwgt = (int(link.getWeight())%10000)/1000.;
    double clwgt = (int(link.getWeight())/10000)/1000.;
    double weight = reco_p.trwgt > 0.5 ? reco_p.trwgt : clwgt;
    if (link.getSim() == mcp && weight > 0.5)  { 
        reco_p.reco = link.getRec();
        reco_p.flag = true;
    }
    else reco_p.flag = false;
    return reco_p;   
}   

// Main function
int track_mclinks(const char* input, const char* output, int pdg_number) {
    gStyle->SetOptStat(0000);
    double Bz = 2;
    std::cout << Bz << std::endl;

    TFile* outputfile = new TFile(output, "RECREATE");

    // To read the generated file
    podio::ROOTReader reader;
    reader.openFile(input);
    
    // Create histograms
    TH1D* void_theta = new TH1D("void_theta", ";#theta [mrad]", 30, 0, 250);
    TH1D* void_pt = new TH1D("void_pt", ";p_{T} [GeV]", 30, 0, 250);
    TH1D* void_p = new TH1D("void_p", ";p [GeV]", 30, 0, 250);
    TH1D* h_pass_theta = new TH1D("h_pass_theta", ";#theta [mrad]", 30, 0, 250);
    TH1D* h_total_theta = new TH1D("h_total_theta", ";#theta [mrad]", 30, 0, 250);
    TH1D* h_pass_pt = new TH1D("h_pass_pt", ";p_{T} [GeV]", 30, 0, 20);
    TH1D* h_total_pt = new TH1D("h_total_pt", ";p_{T} [GeV]", 30, 0, 20);
    TH1D* h_pass_p = new TH1D("h_pass_p", ";p [GeV]", 30, 0, 20);
    TH1D* h_total_p = new TH1D("h_total_p", ";p [GeV]", 30, 0, 20);

    TH1D* h_photon_theta = new TH1D("h_photon_theta", ";#theta [mrad]", 30, 0, 250);
    TH1D* h_proton_theta = new TH1D("h_proton_theta", ";#theta [mrad]", 30, 0, 250);
    TH1D* h_neutron_theta = new TH1D("h_neutron_theta", ";#theta [mrad]", 30, 0, 250);
    TH1D* h_electron_theta = new TH1D("h_electron_theta", ";#theta [mrad]", 30, 0, 250);
    TH1D* h_pion_theta = new TH1D("h_pion_theta", ";#theta [mrad]", 30, 0, 250);
    TH1D* h_other_theta = new TH1D("h_other_theta", ";#theta [mrad]", 30, 0, 250);
    TH1D* h_photon_pt = new TH1D("h_photon_pt", ";p_{T} [GeV]", 30, 0, 20);
    TH1D* h_proton_pt = new TH1D("h_proton_pt", ";p_{T} [GeV]", 30, 0, 20);
    TH1D* h_neutron_pt = new TH1D("h_neutron_pt", ";p_{T} [GeV]", 30, 0, 20);
    TH1D* h_electron_pt = new TH1D("h_electron_pt", ";p_{T} [GeV]", 30, 0, 20);
    TH1D* h_pion_pt = new TH1D("h_pion_pt", ";p_{T} [GeV]", 30, 0, 20);
    TH1D* h_other_pt = new TH1D("h_other_pt", ";p_{T} [GeV]", 30, 0, 20);
    TH1D* h_photon_p = new TH1D("h_photon_p", ";p [GeV]", 30, 0, 20);
    TH1D* h_proton_p = new TH1D("h_proton_p", ";p [GeV]", 30, 0, 20);
    TH1D* h_neutron_p = new TH1D("h_neutron_p", ";p [GeV]", 30, 0, 20);
    TH1D* h_electron_p = new TH1D("h_electron_p", ";p [GeV]", 30, 0, 20);
    TH1D* h_pion_p = new TH1D("h_pion_p", ";p [GeV]", 30, 0, 20);
    TH1D* h_other_p = new TH1D("h_other_p", ";p [GeV]", 30, 0, 20);

    TH1D* h_events_theta = new TH1D("h_events_theta", ";#theta [mrad]", 30, 0, 250);
    TH1D* h_events_pt = new TH1D("h_events_pt", ";p_{T} [GeV]", 30, 0, 20);
    TH1D* h_events_p = new TH1D("h_events_p", ";p [GeV]", 30, 0, 20);

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
            particle = getElectron(mcp, particle, pdg_number);
            if (particle.flag == true) {
                h_total_theta->Fill(particle.theta*1000);
                h_total_pt->Fill(particle.pt);
                h_total_p->Fill(particle.p);
                for (const auto& link : mclinks) {
                    Reconstructed_Particle reco_p;
                    reco_p = getReco(mcp, link, reco_p);
                    if (reco_p.flag == true) { 
                        h_events_theta->Fill(particle.theta *1000);
                        h_events_pt->Fill(particle.pt);
                        h_events_p->Fill(particle.p);
                        if(reco_p.reco.getPDG() == mcp.getPDG() && reco_p.trwgt > 0.75) {
                            h_pass_theta->Fill(particle.theta *1000);
                            h_pass_pt->Fill(particle.pt);
                            h_pass_p->Fill(particle.p);
                        }
                        else {
                            if (std::abs(reco_p.reco.getPDG()) == pdg_number) h_electron_theta->Fill(particle.theta*1000);
                            else if (std::abs(reco_p.reco.getPDG()) == 211) h_pion_theta->Fill(particle.theta*1000);
                            else if (std::abs(reco_p.reco.getPDG()) == 22) h_photon_theta->Fill(particle.theta*1000);
                            else if (std::abs(reco_p.reco.getPDG()) == 2212) h_proton_theta->Fill(particle.theta*1000);
                            else if (std::abs(reco_p.reco.getPDG()) == 2112) h_neutron_theta->Fill(particle.theta*1000);
                            else h_other_theta->Fill(particle.theta*1000);
                    
                            if (std::abs(reco_p.reco.getPDG()) == pdg_number) h_electron_pt->Fill(particle.pt);
                            else if (std::abs(reco_p.reco.getPDG()) == 211) h_pion_pt->Fill(particle.pt);
                            else if (std::abs(reco_p.reco.getPDG()) == 22) h_photon_pt->Fill(particle.pt);
                            else if (std::abs(reco_p.reco.getPDG()) == 2212) h_proton_pt->Fill(particle.pt);
                            else if (std::abs(reco_p.reco.getPDG()) == 2112) h_neutron_pt->Fill(particle.pt);
                            else h_other_pt->Fill(particle.pt);

                            if (std::abs(reco_p.reco.getPDG()) == pdg_number) h_electron_p->Fill(particle.p);
                            else if (std::abs(reco_p.reco.getPDG()) == 211) h_pion_p->Fill(particle.p);
                            else if (std::abs(reco_p.reco.getPDG()) == 22) h_photon_p->Fill(particle.p);
                            else if (std::abs(reco_p.reco.getPDG()) == 2212) h_proton_p->Fill(particle.p);
                            else if (std::abs(reco_p.reco.getPDG()) == 2112) h_neutron_p->Fill(particle.p);
                            else h_other_p->Fill(particle.p);
                        }
                    }
                }
            }
        }
    }
    
    // TEfficency
    TEfficiency* eff_theta = new TEfficiency(*h_pass_theta, *h_total_theta);
    eff_theta->CreateGraph();
    TEfficiency* eff_pt = new TEfficiency(*h_pass_pt, *h_total_pt);
    eff_pt->CreateGraph();
    TEfficiency* eff_p = new TEfficiency(*h_pass_p, *h_total_p);
    eff_p->CreateGraph(); 

    // Percentage of reconstructed particles   
    h_pass_theta->Divide(h_events_theta);
    h_photon_theta->Divide(h_events_theta);  
    h_proton_theta->Divide(h_events_theta);    
    h_neutron_theta->Divide(h_events_theta);
    h_electron_theta->Divide(h_events_theta); 
    h_pion_theta->Divide(h_events_theta);
    h_other_theta->Divide(h_events_theta);

    h_pass_pt->Divide(h_events_pt);
    h_photon_pt->Divide(h_events_pt);  
    h_proton_pt->Divide(h_events_pt);    
    h_neutron_pt->Divide(h_events_pt);
    h_electron_pt->Divide(h_events_pt); 
    h_pion_pt->Divide(h_events_pt);
    h_other_pt->Divide(h_events_pt);

    h_pass_p->Divide(h_events_p);
    h_photon_p->Divide(h_events_p);  
    h_proton_p->Divide(h_events_p);    
    h_neutron_p->Divide(h_events_p);
    h_electron_p->Divide(h_events_p); 
    h_pion_p->Divide(h_events_p);
    h_other_p->Divide(h_events_p);

    // Create graphs 
    eff_theta->SetLineColor(kRed);
    eff_theta->SetMarkerColor(kRed);
    eff_theta->SetMarkerSize(2);
    eff_theta->SetLineWidth(2);
    void_theta->GetYaxis()->SetRangeUser(0.0001, 1.);
    void_theta->GetXaxis()->SetRangeUser(0, 250); 
    void_theta->GetYaxis()->SetTitle("True electron rate");
    void_theta->GetYaxis()->SetTitleSize(0.05);
    void_theta->GetXaxis()->SetTitleSize(0.05);
    void_theta->GetXaxis()->SetTitleOffset(0.8);
    void_theta->GetYaxis()->SetTitleOffset(0.8);

    eff_pt->SetLineColor(kRed);
    eff_pt->SetMarkerColor(kRed);
    eff_pt->SetMarkerSize(2);
    eff_pt->SetLineWidth(2);
    void_pt->GetYaxis()->SetRangeUser(0.0001, 1.);
    void_pt->GetXaxis()->SetRangeUser(0, 20); 
    void_pt->GetYaxis()->SetTitle("True electron rate");
    void_pt->GetYaxis()->SetTitleSize(0.05);
    void_pt->GetXaxis()->SetTitleSize(0.05);
    void_pt->GetXaxis()->SetTitleOffset(0.8);
    void_pt->GetYaxis()->SetTitleOffset(0.8);

    eff_p->SetLineColor(kRed);
    eff_p->SetMarkerColor(kRed);
    eff_p->SetMarkerSize(2);
    eff_p->SetLineWidth(2);
    void_p->GetYaxis()->SetRangeUser(0.0001, 1.);
    void_p->GetXaxis()->SetRangeUser(0, 20); 
    void_p->GetYaxis()->SetTitle("True electron rate");
    void_p->GetYaxis()->SetTitleSize(0.05);
    void_p->GetXaxis()->SetTitleSize(0.05);
    void_p->GetXaxis()->SetTitleOffset(0.8);
    void_p->GetYaxis()->SetTitleOffset(0.8);

    h_photon_theta->SetLineColor(kYellow+1);
    h_photon_theta->SetMarkerColor(kYellow+1);
    h_photon_theta->SetMarkerSize(2);
    h_photon_theta->SetLineWidth(2);
    h_photon_theta->GetYaxis()->SetRangeUser(0.0001, 1.05);
    h_photon_theta->GetXaxis()->SetRangeUser(0, 250); 
    h_photon_theta->GetYaxis()->SetTitle("Reconstructed particles (%)");
    h_photon_theta->GetYaxis()->SetTitleSize(0.05);
    h_photon_theta->GetXaxis()->SetTitleSize(0.05);
    h_photon_theta->GetXaxis()->SetTitleOffset(0.8);
    h_photon_theta->GetYaxis()->SetTitleOffset(0.8);
    h_proton_theta->SetLineColor(kGreen+2);
    h_proton_theta->SetMarkerColor(kGreen+2);
    h_proton_theta->SetMarkerSize(2);
    h_proton_theta->SetLineWidth(2);
    h_neutron_theta->SetLineColor(kCyan+1);
    h_neutron_theta->SetMarkerColor(kCyan+1);
    h_neutron_theta->SetMarkerSize(2);
    h_neutron_theta->SetLineWidth(2);
    h_electron_theta->SetLineColor(kRed-7);
    h_electron_theta->SetMarkerColor(kRed-7);
    h_electron_theta->SetMarkerSize(2);
    h_electron_theta->SetLineWidth(2);
    h_pion_theta->SetLineColor(kBlue+1);
    h_pion_theta->SetMarkerColor(kBlue+1);
    h_pion_theta->SetMarkerSize(2);
    h_pion_theta->SetLineWidth(2);
    h_other_theta->SetLineColor(kRed+2);
    h_other_theta->SetMarkerColor(kRed+2);
    h_other_theta->SetMarkerSize(2);
    h_other_theta->SetLineWidth(2);
    h_pass_theta->SetLineColor(kBlack);
    h_pass_theta->SetMarkerColor(kBlack);
    h_pass_theta->SetMarkerSize(2);
    h_pass_theta->SetLineWidth(2);


    h_photon_pt->SetLineColor(kYellow+1);
    h_photon_pt->SetMarkerColor(kYellow+1);
    h_photon_pt->SetMarkerSize(2);
    h_photon_pt->SetLineWidth(2);
    h_photon_pt->GetYaxis()->SetRangeUser(0.0001, 1.05);
    h_photon_pt->GetXaxis()->SetRangeUser(0, 20); 
    h_photon_pt->GetYaxis()->SetTitle("Reconstructed particles (%)");
    h_photon_pt->GetYaxis()->SetTitleSize(0.05);
    h_photon_pt->GetXaxis()->SetTitleSize(0.05);
    h_photon_pt->GetXaxis()->SetTitleOffset(0.8);
    h_photon_pt->GetYaxis()->SetTitleOffset(0.8);
    h_proton_pt->SetLineColor(kGreen+2);
    h_proton_pt->SetMarkerColor(kGreen+2);
    h_proton_pt->SetMarkerSize(2);
    h_proton_pt->SetLineWidth(2);
    h_neutron_pt->SetLineColor(kCyan+1);
    h_neutron_pt->SetMarkerColor(kCyan+1);
    h_neutron_pt->SetMarkerSize(2);
    h_neutron_pt->SetLineWidth(2);
    h_electron_pt->SetLineColor(kRed-7);
    h_electron_pt->SetMarkerColor(kRed-7);
    h_electron_pt->SetMarkerSize(2);
    h_electron_pt->SetLineWidth(2);
    h_pion_pt->SetLineColor(kBlue+1);
    h_pion_pt->SetMarkerColor(kBlue+1);
    h_pion_pt->SetMarkerSize(2);
    h_pion_pt->SetLineWidth(2);
    h_other_pt->SetLineColor(kRed+2);
    h_other_pt->SetMarkerColor(kRed+2);
    h_other_pt->SetMarkerSize(2);
    h_other_pt->SetLineWidth(2);
    h_pass_pt->SetLineColor(kBlack);
    h_pass_pt->SetMarkerColor(kBlack);
    h_pass_pt->SetMarkerSize(2);
    h_pass_pt->SetLineWidth(2);

    h_photon_p->SetLineColor(kYellow+1);
    h_photon_p->SetMarkerColor(kYellow+1);
    h_photon_p->SetMarkerSize(2);
    h_photon_p->SetLineWidth(2);
    h_photon_p->GetYaxis()->SetRangeUser(0.0001, 1.05);
    h_photon_p->GetXaxis()->SetRangeUser(0, 20); 
    h_photon_p->GetYaxis()->SetTitle("Reconstructed particles (%)");
    h_photon_p->GetYaxis()->SetTitleSize(0.05);
    h_photon_p->GetXaxis()->SetTitleSize(0.05);
    h_photon_p->GetXaxis()->SetTitleOffset(0.8);
    h_photon_p->GetYaxis()->SetTitleOffset(0.8);
    h_proton_p->SetLineColor(kGreen+2);
    h_proton_p->SetMarkerColor(kGreen+2);
    h_proton_p->SetMarkerSize(2);
    h_proton_p->SetLineWidth(2);
    h_neutron_p->SetLineColor(kCyan+1);
    h_neutron_p->SetMarkerColor(kCyan+1);
    h_neutron_p->SetMarkerSize(2);
    h_neutron_p->SetLineWidth(2);
    h_electron_p->SetLineColor(kRed-7);
    h_electron_p->SetMarkerColor(kRed-7);
    h_electron_p->SetMarkerSize(2);
    h_electron_p->SetLineWidth(2);
    h_pion_p->SetLineColor(kBlue+1);
    h_pion_p->SetMarkerColor(kBlue+1);
    h_pion_p->SetMarkerSize(2);
    h_pion_p->SetLineWidth(2);
    h_other_p->SetLineColor(kRed+2);
    h_other_p->SetMarkerColor(kRed+2);
    h_other_p->SetMarkerSize(2);
    h_other_p->SetLineWidth(2);
    h_pass_p->SetLineColor(kBlack);
    h_pass_p->SetMarkerColor(kBlack);
    h_pass_p->SetMarkerSize(2);
    h_pass_p->SetLineWidth(2);

    // Create canvas and draw histograms
    TLatex t;
    t.SetTextSize(30);
    t.SetTextFont(63);
    gStyle->Reset("Modern");
    
    TCanvas* c_track_theta = new TCanvas();    
    void_theta->Draw("p");
    eff_theta->Draw("psame");

    t.DrawLatexNDC(0.15, 0.93935, "CLD #font[52]{work in progress}");
    //gPad->SetLogy();
    c_track_theta->Draw();
    c_track_theta->SaveAs("plots/trk_electrons/ElectronRate_theta_lr.pdf"); 

    TCanvas* c_track_pt = new TCanvas();   
    void_pt->Draw("p");
    eff_pt->Draw("psame");

    t.DrawLatexNDC(0.15, 0.93935, "CLD #font[52]{work in progress}");
    //gPad->SetLogy();
    c_track_pt->Draw();
    c_track_pt->SaveAs("plots/trk_electrons/ElectronRate_pt_lr.pdf");

    TCanvas* c_track_p = new TCanvas();   
    void_p->Draw("p");
    eff_p->Draw("psame");


    t.DrawLatexNDC(0.15, 0.93935, "CLD #font[52]{work in progress}");
    //gPad->SetLogy();
    c_track_p->Draw();
    c_track_p->SaveAs("plots/trk_electrons/ElectronRate_p_lr.pdf");

    TCanvas* c_other_theta = new TCanvas();
    h_photon_theta->Draw("hist");
    h_proton_theta->Draw("histsame");
    h_neutron_theta->Draw("histsame");
    h_electron_theta->Draw("histsame");
    h_pion_theta->Draw("histsame");
    h_other_theta->Draw("histsame");
    h_pass_theta->Draw("histsame");

    t.DrawLatexNDC(0.15, 0.93935, "CLD #font[52]{work in progress}");
    TLegend* l = new TLegend(0.1, 0.5, 0.4, 0.8);
    l->AddEntry(h_photon_theta, "#gamma");
    l->AddEntry(h_proton_theta, "p");
    l->AddEntry(h_neutron_theta, "n");
    l->AddEntry(h_electron_theta, "e^{#pm} incorrectly reconstructed");
    l->AddEntry(h_pion_theta, "#pi^{#pm}");
    l->AddEntry(h_other_theta, "Other particles");
    l->AddEntry(h_pass_theta, "e^{#pm} correctly reconstructed");
    l->Draw();
    //gPad->SetLogy();
    c_other_theta->Draw();
    c_other_theta->SaveAs("plots/trk_electrons/NotReconstructed_theta_lr.pdf");

    TCanvas* c_other_pt = new TCanvas();
    h_photon_pt->Draw("hist");
    h_proton_pt->Draw("histsame");
    h_neutron_pt->Draw("histsame");
    h_electron_pt->Draw("histsame");
    h_pion_pt->Draw("histsame");
    h_other_pt->Draw("histsame");
    h_pass_pt->Draw("histsame");

    t.DrawLatexNDC(0.15, 0.93935, "CLD #font[52]{work in progress}");
    TLegend* l2 = new TLegend(0.1, 0.6, 0.4, 0.9);
    l2->AddEntry(h_photon_pt, "#gamma");
    l2->AddEntry(h_proton_pt, "p");
    l2->AddEntry(h_neutron_pt, "n");
    l2->AddEntry(h_electron_pt, "e^{#pm} incorrectly reconstructed");
    l2->AddEntry(h_pion_pt, "#pi^{#pm}");
    l2->AddEntry(h_other_pt, "Other particles");
    l2->AddEntry(h_pass_theta, "e^{#pm} correctly reconstructed");
    l2->Draw(); 
    //gPad->SetLogy();
    c_other_pt->Draw();
    c_other_pt->SaveAs("plots/trk_electrons/NotReconstructed_pt_lr.pdf");

    TCanvas* c_other_p = new TCanvas();
    h_photon_p->Draw("hist");
    h_proton_p->Draw("histsame");
    h_neutron_p->Draw("histsame");
    h_electron_p->Draw("histsame");
    h_pion_p->Draw("histsame");
    h_other_p->Draw("histsame");
    h_pass_p->Draw("histsame");

    t.DrawLatexNDC(0.15, 0.93935, "CLD #font[52]{work in progress}");
    TLegend* l3 = new TLegend(0.2, 0.65, 0.5, 0.9);
    l3->AddEntry(h_photon_p, "#gamma");
    l3->AddEntry(h_proton_p, "p");
    l3->AddEntry(h_neutron_p, "n");
    l3->AddEntry(h_electron_p, "e^{#pm} incorrectly reconstructed");
    l3->AddEntry(h_pion_p, "#pi^{#pm}");
    l3->AddEntry(h_other_p, "Other particles");
    l3->AddEntry(h_pass_theta, "e^{#pm} correctly reconstructed");
    l3->Draw(); 
    //gPad->SetLogy();
    c_other_p->Draw();
    c_other_p->SaveAs("plots/trk_electrons/NotReconstructed_p_lr.pdf");

    outputfile->Write();
    outputfile->Close(); 
    return 0;
}

