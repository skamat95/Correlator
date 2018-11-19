#ifndef ClusterTrackLinker_hh
#define ClusterTrackLinker_hh

#include "ClusterFinder.hh"
#include "ap_int.h"

#define uint10_t ap_uint<10>
#define uint9_t ap_uint<9>
#define uint3_t ap_uint<3>

const uint16_t MaxTracks = 10;

// Track PT is specified in 16-bits
const uint16_t MaxTrackPT = 0xFFFF;
// Barrel calorimeter eta limit is 1.479 and is measured with LSB = 0.0005 - making it up :)
// Temporarily deal with only positive eta - in principle MaxTrackEta should be double this!
const uint16_t MaxTrackEta = int(1.479 / 0.0005);	//295
const uint16_t conv_cluster_eta = MaxTrackEta / NCrystalsInEta;
const ap_ufixed<10,2> conv_track_eta = 0.2881355; //NCrystalsInEta / MaxTrackEta;
// Barrel calorimeter eta phi coverage is measured with LSB = 0.0001 - making it up :)
const uint16_t MaxTrackPhi = int(2.0 * 3.1415927 / 0.0001);
const uint16_t conv_cluster_phi = MaxTrackPhi / NCrystalsInPhi;
const ap_ufixed<16,2> conv_track_phi = 0.0057296; //NCrystalsInPhi / MaxTrackPhi;

const uint16_t MaxNeutralClusters = NCaloLayer1Eta * NCaloLayer1Phi;


bool getClusterTrackLinker(uint16_t clusterET[NCaloLayer1Eta][NCaloLayer1Phi],
			   uint16_t peakEta[NCaloLayer1Eta][NCaloLayer1Phi],
			   uint16_t peakPhi[NCaloLayer1Eta][NCaloLayer1Phi],
			   uint16_t trackPT[MaxTracks],
			   uint16_t trackEta[MaxTracks],
			   uint16_t trackPhi[MaxTracks],
			   uint16_t linkedTrackPT[MaxTracks],
			   uint16_t linkedTrackEta[MaxTracks],
			   uint16_t linkedTrackPhi[MaxTracks],
			   //float linkedTrackQuality[MaxTracks],
			   uint16_t neutralClusterET[MaxNeutralClusters],
			   uint16_t neutralClusterEta[MaxNeutralClusters],
			   uint16_t neutralClusterPhi[MaxNeutralClusters]);

/*
bool getClusterTrackLinker(uint10_t clusterET[NCaloLayer1Eta][NCaloLayer1Phi],
			   uint3_t peakEta[NCaloLayer1Eta][NCaloLayer1Phi],
			   uint3_t peakPhi[NCaloLayer1Eta][NCaloLayer1Phi],
			   uint10_t trackPT[MaxTracks],
			   uint9_t trackEta[MaxTracks],
			   uint10_t trackPhi[MaxTracks],
			   uint10_t linkedTrackPT[MaxTracks],
			   uint9_t linkedTrackEta[MaxTracks],
			   uint10_t linkedTrackPhi[MaxTracks],
			  // float linkedTrackQuality[MaxTracks],
			   uint10_t neutralClusterET[MaxNeutralClusters],
			   uint9_t neutralClusterEta[MaxNeutralClusters],
			   uint10_t neutralClusterPhi[MaxNeutralClusters]);


bool getClusterTrackLinker(uint10_t clusterET[NCaloLayer1Eta][NCaloLayer1Phi],
		   uint3_t peakEta[NCaloLayer1Eta][NCaloLayer1Phi],
		   uint3_t peakPhi[NCaloLayer1Eta][NCaloLayer1Phi],
		   uint10_t trackPT[MaxTracks],
		   uint9_t trackEta[MaxTracks],
		   uint10_t trackPhi[MaxTracks],
		   algo_out &output);
*/

#endif
