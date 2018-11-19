#include "ClusterFinder.hh"
#include "ClusterTrackLinker.hh"
#include <math.h>


#define uint10_t ap_uint<10>
#define uint9_t ap_uint<9>
#define uint3_t ap_uint<3>

using namespace std;

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
			   uint16_t neutralClusterPhi[MaxNeutralClusters])
			  {

#pragma HLS PIPELINE II=6
//#pragma HLS INTERFACE ap_none port=output
#pragma HLS ARRAY_PARTITION variable=clusterET complete dim=0
#pragma HLS ARRAY_PARTITION variable=peakEta complete dim=0
#pragma HLS ARRAY_PARTITION variable=peakPhi complete dim=0
#pragma HLS ARRAY_PARTITION variable=trackPT complete dim=0
#pragma HLS ARRAY_PARTITION variable=trackEta complete dim=0
#pragma HLS ARRAY_PARTITION variable=trackPhi complete dim=0
#pragma HLS ARRAY_PARTITION variable=linkedTrackPT complete dim=0
#pragma HLS ARRAY_PARTITION variable=linkedTrackEta complete dim=0
#pragma HLS ARRAY_PARTITION variable=linkedTrackPhi complete dim=0
//#pragma HLS ARRAY_PARTITION variable=linkedTrackQuality complete dim=0
#pragma HLS ARRAY_PARTITION variable=neutralClusterET complete dim=0
#pragma HLS ARRAY_PARTITION variable=neutralClusterEta complete dim=0
#pragma HLS ARRAY_PARTITION variable=neutralClusterPhi complete dim=0

  uint16_t clusterEta[MaxNeutralClusters];
  uint16_t clusterPhi[MaxNeutralClusters];
#pragma HLS ARRAY_PARTITION variable=clusterEta complete dim=0
#pragma HLS ARRAY_PARTITION variable=clusterPhi complete dim=0
  for(int tEta = 0; tEta < NCaloLayer1Eta; tEta++) {
#pragma HLS UNROLL
    for(int tPhi = 0; tPhi < NCaloLayer1Phi; tPhi++) {
#pragma HLS UNROLL
      int cluster = tEta * NCaloLayer1Phi + tPhi;
      // Convert cruder calorimeter position to track LSB
      // This can be a LUT - perhaps HLS will take care of this efficiently
      clusterEta[cluster] = (tEta * NCrystalsPerEtaPhi) + peakEta[tEta][tPhi];
      clusterPhi[cluster] = (tPhi * NCrystalsPerEtaPhi) + peakPhi[tEta][tPhi];
      //cout << "peaketaphi" << peakEta[tEta][tPhi] << "\t" << peakPhi[tEta][tPhi] << endl;
      //cout <<"I_from_func" << clusterET[tEta][tPhi] << "\t" << clusterEta[cluster] << "\t" << clusterPhi[cluster] << endl;

      // Initialize neutral clusters
      neutralClusterET[cluster] = clusterET[tEta][tPhi];
      neutralClusterEta[cluster] = clusterEta[cluster];
      neutralClusterPhi[cluster] = clusterPhi[cluster];

     // cout <<"II_from_func" << neutralClusterET[cluster] << "\t" <<neutralClusterEta[cluster] << "\t" << neutralClusterPhi[cluster] << endl;
      //cout <<"III_from_func" << clusterET[tEta][tPhi] << "\t" << clusterEta[cluster] << "\t" << clusterPhi[cluster] << endl;
    }
  }

	ap_fixed<16,8> track_peak_eta[MaxTracks];
	ap_fixed<16,8> track_peak_phi[MaxTracks];
//  float track_peak_eta[MaxTracks];
//  float track_peak_phi[MaxTracks];
#pragma HLS ARRAY_PARTITION variable=track_peak_eta dim=0
#pragma HLS ARRAY_PARTITION variable=track_peak_phi dim=0



  // Double loop over tracks and clusters for linking
  for(int track = 0; track < MaxTracks; track++) {
#pragma HLS UNROLL

	  linkedTrackPT[track] = trackPT[track];
	  linkedTrackEta[track] = trackEta[track];
	  linkedTrackPhi[track] = trackPhi[track];
	  //linkedTrackQuality[track] = 0;//ap_fixed<8,6>(0);


	  track_peak_eta[track] = trackEta[track] * conv_track_eta_calc;
	  track_peak_phi[track] = trackPhi[track] * conv_track_phi_calc;

	  cout<<"track_e"<< trackEta[track] << endl;
	  cout<<"track_p"<< trackPhi[track] << endl;

	  cout<<"track_peak_eta"<< track_peak_eta[track] << endl;
	  cout<<"track_peak_phi"<< track_peak_eta[track] << endl;

//	  int comp_track_peak_eta = int(track_peak_eta[track]);
//	  int comp_track_peak_phi = int(track_peak_phi[track]);

	  uint16_t track_tower_eta = int(track_peak_eta[track] / NCrystalsPerEtaPhi); //Base
	  uint16_t track_tower_phi = int(track_peak_phi[track] / NCrystalsPerEtaPhi);

	  cout<<"base: track_tower_eta"<< track_tower_eta << endl;
	  cout<<"base: track_tower_phi"<< track_tower_phi << endl;

//	  uint16_t track_crystal_eta = int(track_peak_eta[track] % NCrystalsPerEtaPhi); //Offset
//	  uint16_t track_crystal_phi = int(track_peak_phi[track] % NCrystalsPerEtaPhi);
//
//	  uint16_t track_crystal_eta = fmod(track_peak_eta[track],  NCrystalsPerEtaPhi); //Offset
//	  uint16_t track_crystal_phi = fmod(track_peak_phi[track], NCrystalsPerEtaPhi);

	  uint16_t track_crystal_eta = track_peak_eta[track] - (track_tower_eta * NCrystalsPerEtaPhi); //Offset
	  uint16_t track_crystal_phi = track_peak_phi[track] - (track_tower_phi * NCrystalsPerEtaPhi);

	  cout<<"offset: track_crystal_eta"<< track_crystal_eta << endl;
	  cout<<"offset: track_crystal_phi"<< track_crystal_phi << endl;



	  uint16_t diff[3][3];
	  //ap_fixed<10,2> diff[3][3];
	  for(int a = 0; a < 3; a++){
		  for (int b = 0; b<3; b++){
			  diff[a][b] = 0;
		  }
	  }
	  //Calculating crystal tEta and tPhi for track like it is done for clusters
	  uint16_t trackEta = (track_tower_eta * NCrystalsPerEtaPhi) + track_crystal_eta;
	  uint16_t trackPhi = (track_tower_phi * NCrystalsPerEtaPhi) + track_crystal_phi;

	  //Do this to find out which clusters we need to search for this track
	  uint16_t track_tower = (track_tower_eta * NCaloLayer1Phi) + track_tower_phi; //This is the tower corresponding to the track

	  cout<< "track tower" << track_tower << endl;
	  cout << "track_location"<< track_tower << endl;
		if(clusterET[track_tower_eta][track_tower_phi] == 0) break; //If no cluster ET, break

		uint16_t diffEta = clusterEta[track_tower] - track_peak_eta[track];
		if(diffEta >= MaxTrackEta) diffEta = track_peak_eta[track] - clusterEta[track_tower];
		uint16_t diffPhi = clusterPhi[track_tower] - track_peak_phi[track];
		if(diffPhi >= MaxTrackPhi) diffPhi = track_peak_phi[track] - clusterPhi[track_tower];

		//If this is the best match,break:
		if(diffEta <= 1 && diffPhi<= 1)
		{
			//This is the best match, so link them
			linkedTrackEta[track] = clusterEta[track_tower];
			linkedTrackPhi[track] = clusterPhi[track_tower];
			//making neutral cluster
			neutralClusterET[track_tower] -= trackPT[track];
			break;
		}

		//If not, we run for all the other 8 towers in a 3*3 matrix around the track-tower level

		ap_ufixed<10,2> least_dist = 100; //random big value initialization
	  for (int i = -1; i < 2; i++){
#pragma HLS UNROLL
		  for (int j = -1; j < 2; j++){
#pragma HLS UNROLL
			  if(i==0 && j==0) break;


				uint16_t eta = track_tower_eta + i;
				uint16_t phi = track_tower_phi + j;
				//to find out the tower coordinate in the 1D cluster array
				uint16_t cluster_location = (eta * NCaloLayer1Phi) + phi;
				if(cluster_location < 0 || cluster_location > 67) break;

				uint16_t diffEta = clusterEta[cluster_location] - track_peak_eta[track];
				if(diffEta >= MaxTrackEta) diffEta = track_peak_eta[track] - clusterEta[cluster_location];
				uint16_t diffPhi = clusterPhi[cluster_location] - track_peak_phi[track];
				if(diffPhi >= MaxTrackPhi) diffPhi = track_peak_phi[track] - clusterPhi[cluster_location];
				uint16_t common, uncommon;
				if(diffEta >= diffPhi) {

					uncommon = diffEta - diffPhi;
					common = diffEta - uncommon;
					}

				else {
					uncommon = diffPhi - diffEta;
					common = diffPhi - uncommon;
				}

				uint16_t temp = (common>>1) + common;
				diff[i][j] = temp + uncommon;
				//diff[i][j] = ap_ufixed<8,6>(1.4*common + uncommon);
				if(diff[i][j] < least_dist) {
					least_dist =  diff[i][j];
					//This is the best match, so link them
					linkedTrackEta[track] = clusterEta[cluster_location];
					linkedTrackPhi[track] = clusterPhi[cluster_location];
					//making neutral cluster
					neutralClusterET[cluster_location] -= trackPT[track];
					//linkedTrackQuality[track] = least_dist;
				}

		  	}

	  }


  }

	  return true;
}
