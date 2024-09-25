
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

// Struct to hold min radius and z
struct minimum {
    float radius_min;
    float z_min;
};

// Struct to hold max radius and z
struct maximum {
    float radius_max;
    float z_max;
    float energy; // energy at the hit position
    float pt; // pt at the hit position
    float energy_dep; //percentage
};

// Struct to hold MC particle data
struct MC_Particle {
    edm4hep::MCParticle mc_particle;
    int pdg;
    float theta;
    float pt;
    float p;  
    float vertex; //radius
    float z; // vertex
    float energy; 
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
    float y;
    float x;
    float radius;
    float time;   
    float path; 
    float eDep;
    float pt;
    float p; 
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
        tracker_hit.y = pos.y;
        tracker_hit.x = pos.x;
        tracker_hit.radius = sqrt(pos.x*pos.x+pos.y*pos.y+pos.z*pos.z);
        auto mom = hit.getMomentum();
        TVector3 mom_v(mom.x,mom.y,mom.z);
        tracker_hit.pt = mom_v.Pt();
        tracker_hit.p = sqrt(mom.x*mom.x + mom.y*mom.y + mom.z*mom.z); 
        tracker_hit.time = hit.getTime();   
        tracker_hit.path = hit.getPathLength();  
        tracker_hit.eDep = hit.getEDep();
        tracker_hit.flag = true;                                                                            
    }
    return tracker_hit;
}

// Function to get the number of hits of the track
bool getNumberofHitTrack(const edm4hep::SimTrackerHit& hit, const edm4hep::TrackerHitSimTrackerHitLinkCollection& tracker_link, const edm4hep::TrackCollection& tracks) {
    bool flag = false;   
    for (const auto& link : tracker_link) { 
        if(link.getSim() == hit) {
            for (const auto& track : tracks) {
                auto range = track.getTrackerHits(); 
                for (const auto& hits_track: range) {                                                                 
                    if(link.getRec() == hits_track ){                    
                        flag = true;
                    }
                } 
            }
        }
    } 
    return flag;
} 

// Function to get number of the hits of the reconstructed track
bool getNumberofHitReco( const podio::RelationRange<edm4hep::Track>& tracks, const edm4hep::SimTrackerHit& hit,  const edm4hep::TrackerHitSimTrackerHitLinkCollection& tracker_link ){
    bool flag = false;
    for (const auto& link : tracker_link) { 
        if(link.getSim() == hit) {
            for (const auto& track : tracks) {
                auto range = track.getTrackerHits(); 
                for (const auto& hits_track: range) {                                                                 
                    if(link.getRec() == hits_track ){                    
                        flag = true;
                    }
                }
            } 
        }
    } 
    return flag; 
}

// Function that gets the daughters of the MC particle
MC_Particle getDaughters(const edm4hep::MCParticle& mcp, maximum max, minimum min){
    auto daughters = mcp.getDaughters();
    MC_Particle mc_daughter;
    mc_daughter.energy = 0.;
    for (const auto& daughter : daughters){
        float radius_vertex = sqrt(daughter.getVertex().x*daughter.getVertex().x+daughter.getVertex().y*daughter.getVertex().y+daughter.getVertex().z*daughter.getVertex().z);
        float z_vertex = daughter.getVertex().z;
        float momentum_daughter = sqrt(pow(daughter.getMomentum().x,2)+pow(daughter.getMomentum().y,2)+pow(daughter.getMomentum().z,2));
        float energy_daughter = daughter.getEnergy();
        float energy_frac = energy_daughter/max.energy_dep;
        auto mom = daughter.getMomentum(); 
        TVector3 mom_v(mom.x, mom.y, mom.z);
        float pt_daughter = mom_v.Pt();
        float pt_frac = pt_daughter/max.pt;
        if (radius_vertex > max.radius_max && radius_vertex < min.radius_min && abs(z_vertex) > abs(max.z_max) && abs(z_vertex) < abs(min.z_min) && energy_frac > mc_daughter.energy) {
            mc_daughter.mc_particle = daughter;
            mc_daughter.pdg = daughter.getPDG();
            mc_daughter.vertex = radius_vertex; 
            mc_daughter.energy = energy_daughter;
            mc_daughter.pt = pt_frac;
            mc_daughter.flag = true;                 
        }
    }
    return mc_daughter;
}

// Functions that return the minimum/maximum radius of the missing hits
minimum getMinRadius(Tracker_Hit tracker_hit, minimum min){
    if( tracker_hit.radius < min.radius_min) { 
        min.radius_min = tracker_hit.radius;
        min.z_min = tracker_hit.z;
    }
    return min;
}

maximum getMaxRadius(Tracker_Hit tracker_hit, maximum max, const edm4hep::SimTrackerHit& hit, const edm4hep::MCParticle& mcp){
    if( tracker_hit.radius > max.radius_max) { 
        max.radius_max = tracker_hit.radius;
        max.z_max = tracker_hit.z;
        float momentum_mcp = sqrt(pow(hit.getMomentum().x,2)+pow(hit.getMomentum().y,2)+pow(hit.getMomentum().z,2));
        max.energy = sqrt(pow(mcp.getMass(),2) + pow(momentum_mcp,2))-hit.getEDep();
        auto mom = hit.getMomentum(); 
        TVector3 mom_v(mom.x, mom.y, mom.z);
        max.pt = mom_v.Pt();
        max.energy_dep = hit.getEDep()/sqrt(pow(mcp.getMass(),2) + pow(momentum_mcp,2));
    }
    return max;
}
