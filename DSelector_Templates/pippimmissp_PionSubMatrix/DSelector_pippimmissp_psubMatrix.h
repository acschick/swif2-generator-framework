#ifndef DSelector_pippimmissp_psubMatrix_h
#define DSelector_pippimmissp_psubMatrix_h

#include <iostream>

#include "DSelector/DSelector.h"
#include "DSelector/DHistogramActions.h"
#include "DSelector/DCutActions.h"

#include "TH1I.h"
#include "TH2I.h"
#include "TMVA/Reader.h"

class DSelector_pippimmissp_psubMatrix : public DSelector
{
	public:

		DSelector_pippimmissp_psubMatrix(TTree* locTree = NULL) : DSelector(locTree){}
		virtual ~DSelector_pippimmissp_psubMatrix(){}

		void Init(TTree *tree);
		Bool_t Process(Long64_t entry);

	private:
		// TMVA Readers for XML-based inference
		TMVA::Reader* dTMVAReader_MLP;
		TMVA::Reader* dTMVAReader_BDT;
		
		// Input variables for TMVA (must be Float_t for TMVA::Reader)
		Float_t dTMVA_EoverPkin_FCAL;
		Float_t dTMVA_TrackFCAL_DOCA;
		Float_t dTMVA_E9E25_FCAL;
		Float_t dTMVA_E1E9_FCAL;
		Float_t dTMVA_SumU_FCAL;
		Float_t dTMVA_SumV_FCAL;
		Float_t dTMVA_UV_Asymmetry_FCAL;
		Float_t dTMVA_Saturation_FCAL;

		Bool_t ApplyMLPClassification = true;
		Bool_t ApplyBDTClassification = true;
		void Get_ComboWrappers(void);
		void Finalize(void);

		// BEAM POLARIZATION INFORMATION
		UInt_t dPreviousRunNumber;
		bool dIsPolarizedFlag; //else is AMO
		bool dIsPARAFlag; //else is PERP or AMO

		bool dIsMC;

		bool foundBestInTimeCombo;

		// ANALYZE CUT ACTIONS
		// // Automatically makes mass histograms where one cut is missing
		DHistogramAction_AnalyzeCutActions* dAnalyzeCutActions;

		//CREATE REACTION-SPECIFIC PARTICLE ARRAYS

		//Step 0
		DParticleComboStep* dStep0Wrapper;
		DBeamParticle* dComboBeamWrapper;
		DChargedTrackHypothesis* dPositiveWrapper;
		DChargedTrackHypothesis* dNegativeWrapper;
		DKinematicData* dRecoilWrapper;

		// DEFINE YOUR HISTOGRAMS HERE
		// EXAMPLES:
		TH1I* dHist_MissingMassSquared;
		TH1I* dHist_BeamEnergy;
		TH1I* dHist_BeamEnergy_BestChiSq;
		TH1D* dHist_InvMass_BestChiSq_total;
		TH1D* dHist_InvMass_BestChiSq_EoPprecut;
		TH2D* dHist_MLP1_vs_MLP2_BestChiSq_total;
		TH2D* dHist_BDT1_vs_BDT2_BestChiSq_total;
		TH2D* dHist_MLP1_vs_MLP2_BestChiSq_EoPprecut;
		TH2D* dHist_BDT1_vs_BDT2_BestChiSq_EoPprecut;
				Double_t dMinKinFitCL;
		Double_t dMinBeamEnergy;
		Double_t dMaxBeamEnergy;
		Int_t dMaxNumUnusedTracks;

		// Cut enable/disable flags
		Bool_t ApplyPreselectionEoPCut;        // 0a: Preselection E/p > 0.4 (measured) for both tracks
		Bool_t ApplyNoUnusedTracksCut;         // 1a: No unused charged tracks
		Bool_t ApplyNoUnusedShowersCut;        // 1b: No unused neutral shower energy
		Bool_t ApplyFCALEnergyNonZeroCut;      // 2a: FCAL E_1 and E_2 > 0
		Bool_t ApplyTOFdEdxNonZeroCut;         // 2b: TOF dE/dx > 0 for both tracks
		Bool_t ApplyBeamEnergyCut;             // 3a: Coherent peak (8.2-8.8 GeV)
		Bool_t ApplyMinThetaCut;               // 3b: Minimum theta > 1.5 deg
		Bool_t ApplyMaxThetaCut;               // 3c: Maximum theta cut (normally off for rho0)
		Bool_t Apply2DThetaAcceptanceMatchCut; // 3d: e+e- Sim/data theta acceptance match (off for rho0)
		Bool_t ApplyMomentumRangeCut;          // 3e: Valid e+e- FCAL correction range (0.45-7.92 GeV)
		Bool_t ApplyInvariantMassCut;          // 3f: Invariant mass window (0.7-0.77 GeV for rho0)
		Bool_t ApplyVertexZCut;                // 3g: Vertex Z position (52-78 cm)
		Bool_t ApplyKinFitCLCut;               // 3h: Kinematic fit confidence level > 1E-6
		Bool_t OnlyBestChiSqComboInFlat;       // Only save combo with best chiSq to flat tree
		Bool_t ThisComboIsBestChiSq;           // Flag to indicate if the current combo is the best chiSq combo

	// Additional test/development cuts (NOT used in official analysis)
	// These are provided for quick testing and can be enabled as needed
	
	Bool_t ApplyMinEoverPCut;              // Test: Minimum E/P cut for both tracks
	Bool_t ApplyMaxEoverPCut;              // Test: Maximum E/P cut for both tracks
	Bool_t ApplyMinFCALElasticityCut;      // Test: Minimum FCAL elasticity (E1+E2)/Ebeam
	Bool_t ApplyMaxFCALElasticityCut;      // Test: Maximum FCAL elasticity (E1+E2)/Ebeam
	Bool_t ApplyMaxFCALDOCACut;            // Test: Maximum FCAL DOCA cut
	Bool_t ApplyMaxTOFdEdxCut;             // Test: Maximum TOF dE/dx cut

	// Cut threshold values for test cuts
	Double_t dMinEoverP;                   // Minimum E/P for tracks
	Double_t dMaxEoverP;                   // Maximum E/P for tracks
	Double_t dMinFCALElasticity;           // Minimum FCAL elasticity
	Double_t dMaxFCALElasticity;           // Maximum FCAL elasticity  
	Double_t dMaxFCALDOCA;                 // Maximum FCAL DOCA (cm)
	Double_t dMaxTOFdEdx;                  // Maximum TOF dE/dx

	ClassDef(DSelector_pippimmissp_psubMatrix, 0);
};

void DSelector_pippimmissp_psubMatrix::Get_ComboWrappers(void)
{
	//Step 0
	dStep0Wrapper = dComboWrapper->Get_ParticleComboStep(0);
	dComboBeamWrapper = static_cast<DBeamParticle*>(dStep0Wrapper->Get_InitialParticle());
	dPositiveWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(0));
	dNegativeWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(1));
	dRecoilWrapper = dStep0Wrapper->Get_FinalParticle(2);
}

#endif // DSelector_pippimmissp_psubMatrix_h
