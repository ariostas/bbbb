/*
 * diHiggs.C
 * 
 * This macro analyses signal and background samples for bbbb in the boosted regime.
 * To run use "root bbbb.C"
 * 
 * Code by: Andres Rios
 */

// Include ROOT and C++ libraries
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TH1.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <fstream>

// Include Delphes libraries
#include "modules/Delphes.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "classes/DelphesClasses.h"

// Include fastjet
#include "fastjet/ClusterSequence.hh"

using namespace fastjet;
using namespace std;

// Declare functions

void histogram(TH1D*, TH1D*, TCanvas*, const char*, const char*, const char*);
void histogram(TH1D*, TCanvas*, const char*, const char*, const char*);
void analyzeS();
void analyzeB(TString, Double_t);
void saveResults();
Float_t deltaR( const Float_t, const Float_t, const Float_t, const Float_t);

// Initialize histograms



// Initialyze storage variables

Double_t totalSignal = 0, jetsSignal = 0, selJetsSignal = 0, dijetCutSignal = 0, injetCutSignal = 0;
Double_t totalBackground = 0, jetsBackground = 0, selJetsBackground = 0, dijetCutBackground = 0, injetCutBackground = 0;


/*
 * MAIN FUNCTION
 */

void bbbb(){
	
	analyzeS();
	analyzeB("tt-4p-0-600-v1510_14TEV", 530.89358);
	analyzeB("tt-4p-600-1100-v1510_14TEV", 42.55351);
	analyzeB("tt-4p-1100-1700-v1510_14TEV", 4.48209);
	analyzeB("tt-4p-1700-2500-v1510_14TEV", 0.52795);
	analyzeB("tt-4p-2500-100000-v1510_14TEV", 0.05449);
	saveResults();
	
}


/*
 * SIGNAL ANALYSIS
 */

void analyzeS()
{	
	
	const TString inputFilet = "HHToBBBB.txt";
	ifstream ifs(inputFilet);
	assert(ifs.is_open());

	TString filename;
	TChain chain("Delphes");
	
	cout << "Reading signal samples" << endl;

	while(ifs >> filename){
		
		TString filenamef = "root://eoscms.cern.ch//store/user/arapyan/Delphes_phase2/VBFHHTobbbb_TuneZ2_16TeV-madgraph/files/files/" + filename;
		cout << "Reading " << filenamef << endl;
		chain.Add(filenamef);
		
	}
	
	ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
	Long64_t numberOfEntries = treeReader->GetEntries();
	totalSignal = numberOfEntries;

	// Set up loop variables
	Jet *jet;

	// Set up storage variables
	Jet *bJet1, *bJet2, *bJet3, *bJet4;
	Int_t nbJet1, nbJet2, nbJet3, nbJet4;
	TLorentzVector vbJet1, vbJet2, vbJet3, vbJet4;
	
	TClonesArray *branchJet = treeReader->UseBranch("Jet");


	for (Long64_t iEntry=0; iEntry<numberOfEntries; iEntry++) { // Event loop
		treeReader->ReadEntry(iEntry);

		vector<PseudoJet> particles;
		TLorentzVector vJet;
		
		// Use all jets as input for fat jets and choose the four leading bjets
		for (Int_t iJet=0; iJet<branchJet->GetEntries(); iJet++) {// Jet loop
			
			jet = (Jet*) branchJet->At(iJet);
			vJet.SetPtEtaPhiM(jet->PT, jet->Eta, jet->Phi, jet->Mass);
			particles.push_back( PseudoJet(vJet.Px(),  vJet.Px(),  vJet.Px(), vJet.E()) );
			
			if(jet->BTag){
				if(nbJet1 == -1){
					nbJet1 = iJet;
					bJet1 = (Jet*) branchJet->At(nbJet1);
				}
				else if(jet->PT > bJet1->PT){
						
					nbJet4 = nbJet3;
					bJet4 = (Jet*) branchJet->At(nbJet4);
						
					nbJet3 = nbJet2;
					bJet3 = (Jet*) branchJet->At(nbJet3);
						
					nbJet2 = nbJet1;
					bJet2 = (Jet*) branchJet->At(nbJet2);
						
					nbJet1 = iJet;
					bJet1 = (Jet*) branchJet->At(nbJet1);
						
				}
				else if(nbJet2 == -1){
						
					nbJet2 = iJet;
					bJet2 = (Jet*) branchJet->At(nbJet2);
						
				}
				else if(jet->PT > bJet2->PT){
						
					nbJet4 = nbJet3;
					bJet4 = (Jet*) branchJet->At(nbJet4);
						
					nbJet3 = nbJet2;
					bJet3 = (Jet*) branchJet->At(nbJet3);
						
					nbJet2 = iJet;
					bJet2 = (Jet*) branchJet->At(nbJet2);
						
				}
				else if(nbJet3 == -1){
						
					nbJet3 = iJet;
					bJet3 = (Jet*) branchJet->At(nbJet3);
						
				}
				else if(jet->PT > bJet3->PT){
						
					nbJet4 = nbJet3;
					bJet4 = (Jet*) branchJet->At(nbJet4);
						
					nbJet3 = nbJet2;
					bJet3 = (Jet*) branchJet->At(nbJet3);
						
				}
				else if(nbJet4 == -1){
						
					nbJet4 = iJet;
					bJet4 = (Jet*) branchJet->At(nbJet4);
						
				}
				else if(jet->PT > bJet4->PT){
						
					nbJet4 = iJet;
					bJet4 = (Jet*) branchJet->At(nbJet4);
						
				}
			}
			
		}// End jet loop
		
		// Choose a jet definition
		double R = 1.2;
		JetDefinition jet_def_fat(cambridge_algorithm, R);
		
		// Run the clustering, extract the jets
		ClusterSequence cs_fat(particles, jet_def_fat);
		vector<PseudoJet> fat_jets = sorted_by_pt(cs_fat.inclusive_jets());
		
		
		// Check for two fat jets with Pt > 200
		if(!(fat_jets.size() >= 2)) continue;
		if(!((fat_jets[0].pt() > 200) && (fat_jets[1].pt() > 200))) continue;
		
		
		// Check for no isolated leptons
		
		
		
		
		
		// Check for a difference in rapidity of less than 2.5
		if(!(fabs(fat_jets[0].rapidity() - fat_jets[1].rapidity()) < 2.5)) continue;
		
		// Check for four bjets with pt > 40
		if(!((nbJet1 != -1) && (nbJet2 != -1) && (nbJet3 != -1) && (nbJet4 != -1))) continue;
		if(!((bJet1->PT > 40) && (bJet2->PT > 40) && (bJet3->PT > 40) && (bJet4->PT > 40))) continue;
		
		// Construct micro-jets with Cambridge/Aachen algorithm
		R = 0.2
		JetDefinition jet_def_micro(cambridge_algorithm, R);
		ClusterSequence cs_micro(particles, jet_def_micro);
		vector<PseudoJet> micro_jets = sorted_by_pt(cs_micro.inclusive_jets());
		
		// Use BDRS and/or Shower Deconstruction method
		
		
		
		
	} // End event loop

	
	ifs.close();
}


/*
 * BACKGROUND ANALYSIS
 */

void analyzeB(TString input, Double_t crossSection)
{	
	
	
}

/*
 * FUNCTION FOR SAVING THE RESULTS
 */

void saveResults(){

	TCanvas *c1 = new TCanvas("Histogram", "Histogram", 900, 600);
	
	//gROOT->SetOptStat(kFALSE);

	// Save histograms
	
	
	// Print event yields
	

}





/*
 * FUNCTION FOR SAVING TWO HISTOGRAMS
 */
 
void histogram(TH1D *histoS, TH1D *histoB, TCanvas *can, const char* xTitle, const char* yTitle, const char* name){
	Double_t nS=1, nB=1;
	
	nS/=histoS->Integral();
	nB/=histoB->Integral();
	histoS->Scale(nS);
	histoB->Scale(nB);
	
	Double_t max;
	if((histoS->GetMaximum()) > (histoB->GetMaximum())){max=1.1*(histoS->GetMaximum());}
	else if((histoS->GetMaximum()) <= (histoB->GetMaximum())){max=1.1*(histoB->GetMaximum());}
	
	histoS->SetMaximum(max);
	histoB->SetMaximum(max);
	
	histoS->Draw();
	// add axis labels
	histoS->GetXaxis()->SetTitle(xTitle);
	histoS->GetYaxis()->SetTitle(yTitle);
	histoS->SetTitle(""); // title on top
	
	
	histoB->SetLineColor(kRed);
	histoB->Draw("same");
	

	can->SaveAs(name);
}

/*
 * FUNCTION FOR SAVING ONE HISTOGRAM
 */

void histogram(TH1D *histo, TCanvas *can, const char* xTitle, const char* yTitle, const char* name){
	Double_t norm=1;
	norm/=histo->Integral();
	histo->Scale(norm);
	histo->Draw();
	// add axis labels
	histo->GetXaxis()->SetTitle(xTitle);
	histo->GetYaxis()->SetTitle(yTitle);
	histo->SetTitle(""); // title on top

	can->SaveAs(name);
}

/*
 * FUNCTION TO CALCULATE DELTA R
 */

Float_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 ) {

	const Float_t pi = 3.14159265358979;

	Float_t etaDiff = (eta1-eta2);
	Float_t phiDiff = fabs(phi1-phi2);
	while ( phiDiff>pi ) phiDiff = fabs(phiDiff-2.0*pi);

	Float_t deltaRSquared = etaDiff*etaDiff + phiDiff*phiDiff;

	return TMath::Sqrt(deltaRSquared);

}
