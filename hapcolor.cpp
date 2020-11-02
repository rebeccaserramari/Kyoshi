#include "hapcolor.h"
#include <cmath>
#include <algorithm>
#include <math.h>
#include <iomanip>

using namespace std;



float emission(int a, int b) {
	if (a != b) return (float)log(exp(-10));
	else return (float)log(1-exp(-10));
}

float transition(int a, int b) {
	if (a != b) return (float)log(exp(-1));
	else return (float)log(1-exp(-1));
}

float logSum(float x, float y) {
	if( !isfinite(x) || !isfinite(y)) {
		if(!isfinite(x))	return y;
		else return x;
	}
	else {
		if(x > y) return  x + log(1 + exp(y - x));
		else return  y + log(1 + exp(x - y));
	}
	return 0;
}

float fwd_bkw(vector<int> H_A, vector<vector<int>> Ref, string outfilepath) {

	vector<vector<int>> states;
	states = Ref;	
	vector<int> observations;
	observations = H_A;
	
	vector<float> initialprobs;
	for (int i =0; i < states.size(); i++) {
		initialprobs.push_back(log((float)1/states.size()));	
	} 
	float prob = (float)1/states.size();
	prob = log(prob);

	vector<vector<float>> forwardProbabilities;
	vector<float> firstVector;
	for (int i = 0; i < states.size(); i++) {
		firstVector.push_back(initialprobs[i] + emission(states[i][0],observations[0])); 		
	}
	forwardProbabilities.push_back(firstVector); 

	for (int t = 1; t < observations.size(); t++) {
		vector<float> forwardVector;
		for (int j = 0; j < states.size(); j++) {
			float logProb = log(0.0); 
			for (int i = 0; i < states.size(); i++ ) {
				logProb = logSum(logProb, forwardProbabilities[t-1][i] + transition(i, j)); 				
			}
			forwardVector.push_back(logProb + emission(states[j][t],observations[t]));		
		}
		forwardProbabilities.push_back(forwardVector);
	}
	
	vector<vector<float>> backwardProbabilities(observations.size(),vector<float>(states.size()));
	backwardProbabilities[observations.size()-1] = vector<float>(states.size());
	for(int t = observations.size()-2; t > 0 ; t --) {
		vector<float> backwardVector;
		for(int i = 0; i < states.size(); i++ ) {
			float logProb = log(0.0);
			for(int j = 0; j < states.size(); j++) {
				logProb = logSum(logProb, transition(i,j) + backwardProbabilities[t+1][j] + emission(states[j][t+1],observations[t+1]));
			}
			backwardVector.push_back(logProb);
		}
		backwardProbabilities[t] = backwardVector;
	}

	vector<vector<float>> posteriorProbabilities(observations.size(),vector<float>(states.size()));
	for(int t = 0; t < observations.size(); t++) {
		float normalizingConst = log(0.0);
		for(int i = 0; i < states.size(); i++) {
			posteriorProbabilities[t][i] = forwardProbabilities[t][i] + backwardProbabilities[t][i];
			normalizingConst = logSum(normalizingConst, posteriorProbabilities[t][i]);
		}
		for (int i = 0; i < states.size(); i++) {
			posteriorProbabilities[t][i] -= normalizingConst;  
		}
	}
	
	for (int var = 0; var < posteriorProbabilities.size(); var++) {
		float currentmax = posteriorProbabilities[var][0];
		int currentmaxindex = 0;
		for (int h = 0; h <  posteriorProbabilities[var].size(); h++) {
			prob = posteriorProbabilities[var][h];
			if (prob > currentmax) {
				currentmax = prob;
				currentmaxindex = h;	
			}
		}
		int maxcounter = 0;
		for (auto h: posteriorProbabilities[var]) {
			if (h == currentmax) {
				maxcounter += 1;		
			}
		}
	}
	
	ofstream outfile;
	outfile.open(outfilepath);

	for (auto var: posteriorProbabilities) {
		for (auto ref: var) {
			outfile << setprecision(3) << ref << ',';		
		}
		outfile << '\n';
	}

	return(0);
}

/*
	Reads the input and calls the function for computation of the scoring matrix
*/
void compute_hapcolor(vector<int>* hap1, vector<int>* hap2, vector<int>* E_whole, vector<string>* panel, float mutation, vector<int>* final_path1, vector<int>* final_path2, vector<float>* recomblist, string path_posterior_h1, string path_posterior_h2){	
	vector<float>& recombs = *recomblist; 
	vector<int>& H_A = *hap1;
	vector<int>& H_B = *hap2;
	vector<int>& E = *E_whole;
	unordered_set<int> e{begin(E), end(E)};
	
	vector<string>& refpanel = *panel;
	vector<vector<int> > Ref;
	
	//iterate over the strings in reference panel
	for(const auto& row: refpanel) {
		//iterate over values in string
		vector<int> row_values;
		for (unsigned int i = 0; i < row.size(); i++) {
			row_values.push_back(row[i] - '0');	
		}
		Ref.push_back(row_values);
	}	
	
	//read in the mutation parameter
	float param_mut = mutation;
	
	//perform the actual computation of the scoring matrix and return the matrix 'Score'
	float** Score;
	vector<int> stored_positions;	

	time_t start_DP, end_DP;
	time(&start_DP);
	float p_bkw_h1;
	float p_bkw_h2;
	p_bkw_h1 = fwd_bkw(H_A, Ref, path_posterior_h1);	

	p_bkw_h2 = fwd_bkw(H_B, Ref, path_posterior_h2);

}