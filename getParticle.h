
#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>
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

// Struct to hold MC particle data
struct MC_Particle {
    float theta;
    float pt;
    float p;  
    bool flag; 
};

// Struct to hold reconstructed particle data
struct Reconstructed_Particle {
    std::pair<float, float> trwgt;
    std::pair<float,float> clwgt;
    std::pair<edm4hep::ReconstructedParticle,edm4hep::ReconstructedParticle> reco;
};

// Struct to hold hit data
struct Tracker_Hit {
    float z;
    float radius;
    float time;   
    float path;  
    bool flag;                  
};

// Function to get electron data from MCParticleCollection
MC_Particle getStableMCParticle(const edm4hep::MCParticle& mcp, MC_Particle particle, int pdg_number) {  
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
Reconstructed_Particle getReco(const edm4hep::MCParticle& mcp, const edm4hep::MCRecoParticleAssociationCollection& mclinks, Reconstructed_Particle reco_p) { 
    float trwgt = 0.;
    float clwgt = 0.;
    for (const auto& link : mclinks) {  
        float temporary_trwgt = (int(link.getWeight())%10000)/1000.;
        float temporary_clwgt = (int(link.getWeight())/10000)/1000.;
        if (link.getSim() == mcp && temporary_trwgt > trwgt) {    // in the first one the one with the highest track weight
            reco_p.trwgt.first = temporary_trwgt;
            reco_p.reco.first = link.getRec();
            trwgt = reco_p.trwgt.first;           
        }
        if (link.getSim() == mcp && temporary_clwgt > clwgt) {    // in the second one the one with the highest cluster weight
            reco_p.clwgt.second = temporary_clwgt;
            reco_p.reco.second = link.getRec();
            clwgt = reco_p.clwgt.second;           
        }
    }
    return reco_p;   
}  

// Function to get hits
Tracker_Hit getHit(const edm4hep::MCParticle& mcp, const edm4hep::SimTrackerHit& hit, Tracker_Hit tracker_hit) {    
    tracker_hit.flag = false; 
    if (hit.getParticle() == mcp) { 
        auto pos = hit.getPosition(); 
        TVector3 pos_v(pos.x, pos.y, pos.z);
        tracker_hit.z = pos.z;
        tracker_hit.radius = sqrt(pos.x*pos.x+pos.y*pos.y);
        tracker_hit.time = hit.getTime();   
        tracker_hit.path = hit.getPathLength();  
        tracker_hit.flag = true;                                                                            
    }
    return tracker_hit;
}

