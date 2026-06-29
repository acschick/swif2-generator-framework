#ifndef DSelector_pippimmissp_h
#define DSelector_pippimmissp_h

#include <iostream>

#include "DSelector/DSelector.h"
#include "DSelector/DHistogramActions.h"
#include "DSelector/DCutActions.h"

#include "TH1I.h"
#include "TH2I.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"

class DSelector_pippimmissp : public DSelector
{
	public:

		DSelector_pippimmissp(TTree* locTree = NULL) : DSelector(locTree){}
		virtual ~DSelector_pippimmissp(){}

		void Init(TTree *tree);
		Bool_t Process(Long64_t entry);

	private:

		void Get_ComboWrappers(void);
		void Finalize(void);

		// BEAM POLARIZATION INFORMATION
		UInt_t dPreviousRunNumber;
		bool dIsPolarizedFlag; //else is AMO
		bool dIsPARAFlag; //else is PERP or AMO
		int dPolarizationAngle;

		bool dIsMC;

		// ANALYZE CUT ACTIONS
		// // Automatically makes mass histograms where one cut is missing
		DHistogramAction_AnalyzeCutActions* dAnalyzeCutActions;

		//CREATE REACTION-SPECIFIC PARTICLE ARRAYS

		//Step 0
		DParticleComboStep* dStep0Wrapper;
		DBeamParticle* dComboBeamWrapper;
		DChargedTrackHypothesis* dPositiveWrapper;  // Generic: π+, K+, e+, etc.
		DChargedTrackHypothesis* dNegativeWrapper;  // Generic: π-, K-, e-, etc.
		DKinematicData* dRecoilWrapper;  // Generic: missing proton, detected proton, etc.		// DEFINE YOUR HISTOGRAMS HERE
		// EXAMPLES:

		// DEFINE YOUR HISTOGRAMS HERE
		// EXAMPLES:
		TH1I* dHist_MissingMassSquared;
		TH1I* dHist_BeamEnergy;
		TH1I* dHist_BeamEnergy_BestChiSq;
		TH1I* dHist_TaggerAccidentals;

		TH1D* dHist_InvMass_TwoTrack;
		TH1D* dHist_InvMass_TwoTrack_BestChiSq;

		TH1D* dHist_FCAL_Energy_pip;
		TH1D* dHist_FCAL_Energy_pim;
		TH1D* dHist_FCAL_EoverP_pip;
		TH1D* dHist_FCAL_EoverP_pim;
		TH1D* dHist_FCAL_EoverPmeas_pip;
		TH1D* dHist_FCAL_EoverPmeas_pim;
		TH2D* dHist_Delta_Efcal_kinfitE_vs_kinPmag_plus;
		TH2D* dHist_Delta_Efcal_kinfitE_vs_kinPmag_minus;
		TH2D* dHist_Delta_Efcal_kinfitE_vs_kinPtheta_plus;
		TH2D* dHist_Delta_Efcal_kinfitE_vs_kinPtheta_minus;
		TH1D* dHist_FCAL_Elasticity;
		TH1D* dHist_FCAL_Asymmetry;
		TH1D* dHist_FCAL_Elasticity_Asym_L0pt2;
		TH1D* dHist_FCAL_Elasticity_Asym_G0pt2_L0pt5;
		TH1D* dHist_FCAL_Elasticity_Asym_G0pt5;
		TH1D* dHist_TrackFCAL_DOCA_plus;
		TH1D* dHist_TrackFCAL_DOCA_minus;
		TH1D* dHist_FCAL_E1E9_plus;
		TH1D* dHist_FCAL_E1E9_minus;
		TH1D* dHist_FCAL_E9E25_plus;
		TH1D* dHist_FCAL_E9E25_minus;
		TH1D* dHist_FCAL_kin_res_plus;
		TH1D* dHist_FCAL_kin_res_minus;
		TH1D* dHist_FCAL_meas_res_plus;
		TH1D* dHist_FCAL_meas_res_minus;
		TH2D* dHist_FCAL_Saturation_vs_Eshower_plus;
		TH2D* dHist_FCAL_Saturation_vs_Eshower_minus;
		TH1D* dHist_FCAL_Saturation_plus;
		TH1D* dHist_FCAL_Saturation_minus;
		TH1D* dHist_FCAL_SumU_plus;
		TH1D* dHist_FCAL_SumU_minus;
		TH1D* dHist_FCAL_SumV_plus;
		TH1D* dHist_FCAL_SumV_minus;
		TH1D* dHist_FCAL_UV_Asymmetry_plus;
		TH1D* dHist_FCAL_UV_Asymmetry_minus;

		// Polarization angles to analyze
		static const int nPolarizations = 4;
		static constexpr const char* Polarizations[nPolarizations] = {"0", "45", "90", "135"};

		// epemPol directory histograms (electron mass hypothesis) - indexed by polarization
		TH1D* dHist_epemPol_InvMass_Kin;
		TH1D* dHist_epemPol_InvMass_Meas;
		TH1D* dHist_epemPol_JTphi_Kin[nPolarizations];
		TH1D* dHist_epemPol_JTphi_Meas[nPolarizations];
		TH1D* dHist_epemPol_qvec2;

		// Rho0Pol directory histograms (pion mass hypothesis) - indexed by polarization
		TH1D* dHist_Rho0Pol_InvMass_Kin;
		TH1D* dHist_Rho0Pol_InvMass_Meas;
		TH1D* dHist_Rho0Pol_CosThetaKin[nPolarizations];
		TH1D* dHist_Rho0Pol_CosThetaMeas[nPolarizations];
		TH1D* dHist_Rho0Pol_phiKin[nPolarizations];
		TH1D* dHist_Rho0Pol_phiMeas[nPolarizations];
		TH1D* dHist_Rho0Pol_PhiKin[nPolarizations];
		TH1D* dHist_Rho0Pol_PhiMeas[nPolarizations];
		TH1D* dHist_Rho0Pol_psiKin[nPolarizations];
		TH1D* dHist_Rho0Pol_psiMeas[nPolarizations];
		TH2D* dHist_Rho0Pol_phikin_vs_Phikin[nPolarizations];


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
		
		


	ClassDef(DSelector_pippimmissp, 0);
};

void DSelector_pippimmissp::Get_ComboWrappers(void)
{
	//Step 0
	dStep0Wrapper = dComboWrapper->Get_ParticleComboStep(0);
	dComboBeamWrapper = static_cast<DBeamParticle*>(dStep0Wrapper->Get_InitialParticle());
	dPositiveWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(0));
	dNegativeWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(1));
	dRecoilWrapper = dStep0Wrapper->Get_FinalParticle(2);
}

#endif // DSelector_pippimmissp_h
