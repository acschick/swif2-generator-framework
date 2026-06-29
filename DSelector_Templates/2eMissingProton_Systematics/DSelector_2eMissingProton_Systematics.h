#ifndef DSelector_2eMissingProton_Systematics_h
#define DSelector_2eMissingProton_Systematics_h

#include <iostream>
#include <string>

#include "DSelector/DSelector.h"
#include "DSelector/DHistogramActions.h"
#include "DSelector/DCutActions.h"
#include "TMVA/Reader.h"

#include "TH1I.h"
#include "TH2I.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TRandom.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"

class DSelector_2eMissingProton_Systematics : public DSelector
{
	public:

		DSelector_2eMissingProton_Systematics(TTree* locTree = NULL) : DSelector(locTree){}
		virtual ~DSelector_2eMissingProton_Systematics(){}

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
		
		struct ComboCutInputs {
			Double_t beamE;
			bool ThisComboIsBestChiSq;
			Double_t theta1;
			Double_t theta2;
			Double_t Wkin;
			Double_t chisqndf;
			Double_t kinFitCL;
			Double_t vertexZ;
			Double_t pMagPlus;
			Double_t pMagMinus;
			Int_t numUnusedTracks;
			Double_t energyUnusedShowers;
			Double_t FCAL_Energy_plus;
			Double_t FCAL_Energy_minus;
			Double_t TOF_dEdx_plus;
			Double_t TOF_dEdx_minus;
			Double_t FCAL_EoverPkin_plus;
			Double_t FCAL_EoverPkin_minus;
			Double_t FCAL_EoverPmeas_plus;
			Double_t FCAL_EoverPmeas_minus;
			Double_t FCAL_Elasticity;
			Double_t TrackFCAL_DOCA_plus;
			Double_t TrackFCAL_DOCA_minus;
			Double_t MLP1;
			Double_t MLP2;
			Double_t BDT1;
			Double_t BDT2;
		};

		struct CutConditions {
			// Section 0: Preselection cut emulation
			bool applyPreselectionEoPCut;         // 0a: Preselection E/p > 0.4 (measured) for both tracks
				double preselectionMinEoverP;

			// Section 1: Exclusivity cuts
			bool applyNoUnusedTracksCut;          // 1a: No unused charged tracks
				int maxNumUnusedTracks;               	// == 0
			bool applyNoUnusedShowersCut;         // 1b: No unused neutral shower energy
				double maxUnusedShowersEnergy;       	 // == 0

			// Section 2: Forward detector cuts
			bool applyFCALEnergyNonZeroCut;       // 2a: FCAL E_1 and E_2 > 0
			bool applyTOFdEdxNonZeroCut;          // 2b: TOF dE/dx > 0 for both tracks

			// Section 3: Fiducial cuts
			bool applyMinBeamEcut;                // 3a: Form Factor: Full Spectrum (7.0-11.0 GeV)
			bool applyMaxBeamEcut;                // 3a: Polarization: Coherent Peak (8.2-8.8 GeV) 
				double minBeamE;					  // 8.2 GeV if Coherent Peak, 7.0 if Full Spectrum
				double maxBeamE;					  // 8.8 GeV if Coherent Peak, 11.0 if Full Spectrum
			bool applyMinThetaCut;                // 3b: Minimum theta > 1.5 deg
			bool applyMaxThetaCut;                // 3c: Maximum theta cut (turn off for rho0 pion sub sample)
				double minTheta;					  // 1.5 deg if GlueX, 1.1 if CPP
				double maxTheta;                      // 7.5 deg if GlueX, 7.5 if CPP 
			                                     
			bool apply2DThetaCut;                 // 3d: Sim/data theta acceptance match (off for rho0 pion sample)
			bool applyMomentumRangeCut;           // 3e: Valid FCAL correction range (0.45-7.92 GeV)
				double minPmag;						  // 0.45 GeV/c 
				double maxPmag;                       // 7.92 GeV/c 

			bool applyMinWkinCut;                 // 3f: Invariant mass window: 0.25  - 0.621 GeV for e+e- analysis
			bool applyMaxWkinCut;                 //                            0.700 - 0.770 GeV for rho0
				double minWkin;						  // 0.25  GeV 
				double maxWkin;				          // 0.621 GeV 
			bool applyVertexZCut;                 // 3g: Vertex Z position (52-78 cm)
				double minVertexZ;					  // 52 cm
				double maxVertexZ;                    // 78 cm
			bool applyKinFitCLCut;                // 3h: Kinematic fit confidence level > 1E-6
				double minKinFitCL;					  // 1E-6	
			bool applyBestChiSqComboCut;         // 3i: Require selected best in-time combo only

			// Section 4: Multivariate Analysis Cuts
			bool applyMVACuts;
		std::string modelChoice;             // "MLP" or "BDT"
		std::string particleChoice;			 // "ee" or "pi" or "none". none automatically sets applyMVACuts = false.
				double BDT_ee;   				 	// > 0.053080951 for e+e-
				double BDT_pi;   				 	// < -0.069788195 for pi+pi-
				double MLP_ee;   					// > 0.8 for e+e-
				double MLP_pi;   					// < 0.4 for pi+pi-

			// Additional test/development cuts (NOT used in official analysis)
			bool applyMinEoverP1Cut;              // T1a: Minimum E/P cut for track 1
				double minEoverP1;
			bool applyMaxEoverP1Cut;              // T1b: Maximum E/P cut for track 1
				double maxEoverP1;
			bool applyMinEoverP2Cut;              // T1c: Minimum E/P cut for track 2
				double minEoverP2;
			bool applyMaxEoverP2Cut;              // T1d: Maximum E/P cut for track 2
				double maxEoverP2;
			bool applyMinFCALElasticityCut;       // T2a: Minimum FCAL elasticity (E1+E2)/Ebeam
				double minFCALElasticity;
			bool applyMaxFCALElasticityCut;       // T2b: Maximum FCAL elasticity (E1+E2)/Ebeam
				double maxFCALElasticity;
			bool applyMinFCALEnergyCut;           // T2c: Minimum FCAL energy for each track
				double minFCALEnergy;
			bool applyMaxFCALDOCACut;             // T3a: Maximum FCAL DOCA cut
				double maxFCALDOCA;
			bool applyMinTOFdEdxCut;			  // T4a: Minimum TOF dE/dx cut
				double minTOFdEdx;
			bool applyMaxTOFdEdxCut;              // T4b: Maximum TOF dE/dx cut
				double maxTOFdEdx;
			bool applyMLPResponseCut;             // T5a: Minimum MLP response cut for e+/e- identification
				double minMLPResponse;
			bool applyChiSqCut;			          // T5b: Max kin fit ChiSq/ndf cut (off in official analysis, 1e-6 cut on CL applied in earlier stage)
				double maxKFChiSq;					
		};


		void Get_ComboWrappers(void);
		void Finalize(void);
		void BookSystematics(void);
		void FillSystematics(const ComboCutInputs& inputs, Double_t q2kin, Double_t JTphi, Double_t Wmeas, Double_t q2kinRes, bool includeQ2ResSystematics, int runPeriodIndex, int polarizationIndex, Double_t fillWeight = 1.0);
		bool ComboPassesCuts(const ComboCutInputs& inputs, const CutConditions& conditions);
		CutConditions BuildActiveFiducialConditions(void) const;
	int GetBeamWindowIndex(double beamE);  // Returns 0=FullSpectrum, 1=CoherentPeak, -1=outside ranges

		// BEAM POLARIZATION INFORMATION
		UInt_t dPreviousRunNumber;
		bool dIsCPPRunPeriod;  // true for CPP run period (100531-101622), false for GlueX-I
		int dRunPeriodIndex;   // 0=1801, 1=1808, 2=CPP, -1=unknown/outside configured ranges
		bool dIsPolarizedFlag; //else is AMO
		bool dIsPARAFlag; //else is PERP or AMO
		int dPolarizationAngle;

		bool dIsMC;  // Set once in Init() - true if analyzing simulated data

		// ANALYZE CUT ACTIONS
		// // Automatically makes mass histograms where one cut is missing
		//DHistogramAction_AnalyzeCutActions* dAnalyzeCutActions;

		// Centralized systematics directory structure (created once in Init)
		TDirectory* dDir_Systematics;
		TDirectory* dDir_Systematics_JTphi;
		TDirectory* dDir_Systematics_q2;
		TDirectory* dDir_Systematics_theta;
		TDirectory* dDir_Systematics_invmass;
		TDirectory* dDir_Systematics_chisq;
		TDirectory* dDir_Systematics_thrown;
		TDirectory* dDir_MVAResponsePlots;

		// Systematics configuration
		static const int kNumMethods = 2;
		static const int kNumPid = 3;
		static const int kNumAngles = 27;
		static const int kNumBeamRegions = 16;
		static const int kNumMassRegions = 10;
		static const int kNumWminCuts = 16;
		static const int kNumChiSqCuts = 16;
		static const int kNumMLPeeThresholds = 9;
		static const int kNumMLPpiThresholds = 9;
		static const int kNumBDTeeThresholds = 14;
		static const int kNumBDTpiThresholds = 18;
		static const int kNumBestChiSqCategories = 6;
		static const int kNumBeamWindows = 2;
		static const int kNumRunPeriods = 3;
		static const int nPolarizations = 5;

		CutConditions dAnalysisCutConditions;
		CutConditions dAnalysisCutConditionsCPP;
		const char* dMethods[kNumMethods];
		const char* dPidChoices[kNumPid];
		double dMinAngles[kNumAngles];
		double dBeamEnergyRegions[kNumBeamRegions][2];
		int dNumActiveBeamRegions;
		bool dBookCPPBeamRegions;
		double dMassRegions[kNumMassRegions][2];
		double dMinInvMass[kNumWminCuts];
		double dMaxChiSq[kNumChiSqCuts];
		float dMLPeeThresholds[kNumMLPeeThresholds];
		float dMLPpiThresholds[kNumMLPpiThresholds];
		float dBDTeeThresholds[kNumBDTeeThresholds];
		float dBDTpiThresholds[kNumBDTpiThresholds];

		// Systematics histograms
		TH2D* dHist_theta2_vs_theta1_2DCutOFF[kNumMethods][kNumPid];
		TH1D* dHist_theta1_noCuts[kNumMethods][kNumPid];
		TH1D* dHist_theta1_fidCuts_noTheta[kNumMethods][kNumPid];
		TH1D* dHist_theta2_noCuts[kNumMethods][kNumPid];
		TH1D* dHist_theta2_fidCuts_noTheta[kNumMethods][kNumPid];
		TH1D* dHist_JTphi_angles[kNumMethods][kNumPid][kNumAngles];
		TH1D* dHist_JTphi_FID_2DCutOFF[kNumMethods][kNumPid];
		TH1D* dHist_q2kin_angles[kNumMethods][kNumPid][kNumAngles];
		TH1D* dHist_q2kin_varWidth_angles[kNumMethods][kNumPid][kNumAngles];
		TH1D* dHist_q2kin_FID_2DCutOFF[kNumMethods][kNumPid];
		TH1D* dHist_q2kin_varWidth_FID_2DCutOFF[kNumMethods][kNumPid];
		TH2D* dHist_qvec2_res_vs_q2kin_FID_2DCutOFF[kNumMethods][kNumPid];
		TH2D* dHist_qvec2_res_vs_q2kin_angles[kNumMethods][kNumPid][kNumAngles];

		TH1D* dHist_JTphi_BeamRegions[kNumMethods][kNumPid][kNumBeamRegions];
		TH1D* dHist_q2kin_BeamRegions[kNumMethods][kNumPid][kNumBeamRegions];
		TH1D* dHist_q2kin_varWidth_BeamRegions[kNumMethods][kNumPid][kNumBeamRegions];
		TH2D* dHist_qvec2_res_vs_q2kin_BeamRegions[kNumMethods][kNumPid][kNumBeamRegions];

		TH1D* dHist_Wepem_kin_noCuts[kNumMethods][kNumPid];
		TH1D* dHist_Wepem_kin_fidCuts_noWcut[kNumMethods][kNumPid];
		TH1D* dHist_Wepem_Measured_noCuts[kNumMethods][kNumPid];
		TH1D* dHist_Wepem_Measured_fidCuts_noWcut[kNumMethods][kNumPid];
		TH1D* dHist_Wkin_MassRegions[kNumMethods][kNumPid][kNumMassRegions];
		TH1D* dHist_JTphi_MassRegions[kNumMethods][kNumPid][kNumMassRegions];
		TH1D* dHist_JTphi_WminCuts[kNumMethods][kNumPid][kNumWminCuts];
		TH1D* dHist_q2kin_MassRegions[kNumMethods][kNumPid][kNumMassRegions];
		TH1D* dHist_q2kin_varWidth_MassRegions[kNumMethods][kNumPid][kNumMassRegions];
		TH1D* dHist_q2kin_WminCuts[kNumMethods][kNumPid][kNumWminCuts];
		TH1D* dHist_q2kin_varWidth_WminCuts[kNumMethods][kNumPid][kNumWminCuts];
		TH2D* dHist_qvec2_res_vs_q2kin_MassRegions[kNumMethods][kNumPid][kNumMassRegions];
		TH2D* dHist_qvec2_res_vs_q2kin_WminCuts[kNumMethods][kNumPid][kNumWminCuts];

		TH1D* dHist_KinFitChiSq_noCuts[kNumMethods][kNumPid];
		TH1D* dHist_KinFitChiSq_fidCuts[kNumMethods][kNumPid];
		TH1D* dHist_JTphi_MaxChiSq[kNumMethods][kNumPid][kNumChiSqCuts];
		TH1D* dHist_q2kin_MaxChiSq[kNumMethods][kNumPid][kNumChiSqCuts];
		TH1D* dHist_q2kin_varWidth_MaxChiSq[kNumMethods][kNumPid][kNumChiSqCuts];
		TH2D* dHist_qvec2_res_vs_q2kin_MaxChiSq[kNumMethods][kNumPid][kNumChiSqCuts];

		TH1D* dHist_JTphi_MLP_ee[kNumMLPeeThresholds];
		TH1D* dHist_JTphi_MLP_pi[kNumMLPpiThresholds];
		TH1D* dHist_q2kin_MLP_ee[kNumMLPeeThresholds];
		TH1D* dHist_q2kin_MLP_pi[kNumMLPpiThresholds];
		TH1D* dHist_q2kin_varWidth_MLP_ee[kNumMLPeeThresholds];
		TH1D* dHist_q2kin_varWidth_MLP_pi[kNumMLPpiThresholds];
		TH2D* dHist_qvec2_res_vs_q2kin_MLP_ee[kNumMLPeeThresholds];
		TH2D* dHist_qvec2_res_vs_q2kin_MLP_pi[kNumMLPpiThresholds];

		TH1D* dHist_JTphi_BDT_ee[kNumBDTeeThresholds];
		TH1D* dHist_JTphi_BDT_pi[kNumBDTpiThresholds];
		TH1D* dHist_q2kin_BDT_ee[kNumBDTeeThresholds];
		TH1D* dHist_q2kin_BDT_pi[kNumBDTpiThresholds];
		TH1D* dHist_q2kin_varWidth_BDT_ee[kNumBDTeeThresholds];
		TH1D* dHist_q2kin_varWidth_BDT_pi[kNumBDTpiThresholds];
		TH2D* dHist_qvec2_res_vs_q2kin_BDT_ee[kNumBDTeeThresholds];
		TH2D* dHist_qvec2_res_vs_q2kin_BDT_pi[kNumBDTpiThresholds];

		// Run-period and polarization labels for polarization-dependent histogram families

		static constexpr const char* RunPeriodTags[kNumRunPeriods] = {"1801", "1808", "CPP"};
		
		static constexpr const char* Polarizations[nPolarizations] = {"0", "45", "90", "135", "AMO"};

		TH1D* dHist_MLP_responsePIP;
		TH1D* dHist_MLP_responsePIM;
		TH2D* dHist_MLP_response_PIP_vs_PIM;
		TH1D* dHist_BDT_responsePIP;
		TH1D* dHist_BDT_responsePIM;
		TH2D* dHist_BDT_response_PIP_vs_PIM;
		// Thrown and Resolution histograms - split by beam window
		TH2D* dHist_qvec2_res_vs_q2kin[kNumBeamWindows];
		TH1D* dHist_qvec2_varWidth_Thrown[kNumBeamWindows];
		TH1D* dHist_qvec2_Thrown[kNumBeamWindows];
		TH2D* dHist_theta_KinRes_vs_theta_Thrown[kNumBeamWindows];
		TH1D* dHist_theta1_Thrown[kNumBeamWindows];
		TH1D* dHist_theta2_Thrown[kNumBeamWindows];
		TH1D* dHist_ep_Pmag_Thrown[kNumBeamWindows];
		TH1D* dHist_em_Pmag_Thrown[kNumBeamWindows];
		TH1D* dHist_ep_Pmag_KinRes[kNumBeamWindows];
		TH1D* dHist_em_Pmag_KinRes[kNumBeamWindows];
		TH2D* dHist_ep_Pmag_KinRes_vs_ep_Pmag_Thrown[kNumBeamWindows];
		TH2D* dHist_em_Pmag_KinRes_vs_em_Pmag_Thrown[kNumBeamWindows];
		TH1D* dHist_Wepem_Thrown[kNumBeamWindows];
		TH1D* dHist_Wepem_KinRes[kNumBeamWindows];
		TH2D* dHist_Wepem_KinRes_vs_Wepem_Thrown[kNumBeamWindows];
		TH1D* dHist_JTphi_Thrown[kNumBeamWindows][kNumRunPeriods][nPolarizations];
		TH1D* dHist_JTphi_kinRes[kNumBeamWindows][kNumRunPeriods][nPolarizations];
		TH1D* dHist_ep_phi_Thrown[kNumBeamWindows][kNumRunPeriods][nPolarizations];
		TH1D* dHist_em_phi_Thrown[kNumBeamWindows][kNumRunPeriods][nPolarizations];
		TH1D* dHist_phi_KinRes[kNumBeamWindows];

		// Recoil proton diagnostics: theta vs momentum
		TH2D* dHist_RecoilThetaVsP_Thrown[kNumBeamWindows];

		// N-1 plots for fiducial cut validation - split by beam window
		TH1I* dHist_FID_Nminus1_BeamEnergy[kNumBeamWindows];
		TH1D* dHist_FID_Nminus1_Wepem[kNumBeamWindows];
		TH1D* dHist_FID_Nminus1_p_ep[kNumBeamWindows];
		TH1D* dHist_FID_Nminus1_p_em[kNumBeamWindows];
		TH1D* dHist_FID_Nminus1_theta1[kNumBeamWindows];
		TH1D* dHist_FID_Nminus1_theta2[kNumBeamWindows];
		TH2D* dHist_FID_Nminus1_theta2_vs_theta1[kNumBeamWindows];
		TH1D* dHist_FID_Nminus1_VertexZ[kNumBeamWindows];
		
		// Theta variants - no theta cuts applied
		TH1D* dHist_FID_Nminus1_theta1_noThetaCuts[kNumBeamWindows];
		TH1D* dHist_FID_Nminus1_theta2_noThetaCuts[kNumBeamWindows];
		TH2D* dHist_FID_Nminus1_theta2_vs_theta1_noThetaCuts[kNumBeamWindows];
		
		// Preselection E/p
		TH1D* dHist_FID_Nminus1_EoverP_ep[kNumBeamWindows];
		TH1D* dHist_FID_Nminus1_EoverP_em[kNumBeamWindows];
		
		// Exclusivity
		TH1I* dHist_FID_Nminus1_NumUnusedTracks[kNumBeamWindows];
		TH1D* dHist_FID_Nminus1_UnusedShowerEnergy[kNumBeamWindows];
		
		// Forward detector
		TH1D* dHist_FID_Nminus1_FCALEnergy_ep[kNumBeamWindows];
		TH1D* dHist_FID_Nminus1_FCALEnergy_em[kNumBeamWindows];
		TH1D* dHist_FID_Nminus1_TOFdEdx_ep[kNumBeamWindows];
		TH1D* dHist_FID_Nminus1_TOFdEdx_em[kNumBeamWindows];
		
		// MVA outputs
		TH1D* dHist_FID_Nminus1_MLP_ep[kNumBeamWindows];
		TH1D* dHist_FID_Nminus1_MLP_em[kNumBeamWindows];
		TH2D* dHist_FID_Nminus1_MLP_ep_vs_em[kNumBeamWindows];
		TH1D* dHist_FID_Nminus1_BDT_ep[kNumBeamWindows];
		TH1D* dHist_FID_Nminus1_BDT_em[kNumBeamWindows];
		TH2D* dHist_FID_Nminus1_BDT_ep_vs_em[kNumBeamWindows];
		
		// Special no-MVA variants
		TH1D* dHist_FID_Nminus1_EoverP_ep_and_noMVA[kNumBeamWindows];
		TH1D* dHist_FID_Nminus1_FCALElasticity_and_noMVA[kNumBeamWindows];

		//CREATE REACTION-SPECIFIC PARTICLE ARRAYS

		//Step 0
		DParticleComboStep* dStep0Wrapper;
		DBeamParticle* dComboBeamWrapper;
		DChargedTrackHypothesis* dPositiveWrapper;  // Generic: e+, π+, K+, etc.
		DChargedTrackHypothesis* dNegativeWrapper;  // Generic: e-, π-, K-, etc.
		DKinematicData* dRecoilWrapper;  // Generic: missing proton, detected proton, etc.
	
		// DEFINE YOUR HISTOGRAMS HERE
		// EXAMPLES:
		TH1I* dHist_MissingMassSquared;
		TH1I* dHist_BeamEnergy;
		TH1I* dHist_BeamEnergy_BestChiSq;
		TH1I* dHist_TaggerAccidentals;
		TH1D* dHist_KinFitChiSq;
		TH1D* dHist_KinFitCL;
		TH1D* dHist_VertexZ_BestChiSq;

		TH1D* dHist_Wepem;
		TH1D* dHist_Wepem_BestChiSq;
		TH1I* dHist_MissingMassSquared_Measured;
		TH1I* dHist_MissingMassSquared_Measured_AllowOneUnused;
		TH1I* dHist_MissingEnergy_Measured;
		TH2D* dHist_preFid_MM2Residual_vs_RecoilP_UpTo1UnusedTrack;
		TH1D* dHist_EoverP_measured_plus;
		TH1D* dHist_EoverP_measured_minus;
		TH1D* dHist_MLPResponse_plus;
		TH1D* dHist_MLPResponse_minus;
		TH2D* dHist_MLPResponsePlus_vs_MLPResponseMinus;

		TH1D* dHist_FCAL_Energy_pip;
		TH1D* dHist_FCAL_Energy_pim;
		TH1D* dHist_FCAL_Energy_pip_PostCorr;
		TH1D* dHist_FCAL_Energy_pim_PostCorr;
		TH1D* dHist_FCAL_EoverP_pip;
		TH1D* dHist_FCAL_EoverP_pim;
		TH1D* dHist_FCAL_EoverP_pip_PostCorr;
		TH1D* dHist_FCAL_EoverP_pim_PostCorr;
		TH1D* dHist_FCAL_EoverPmeas_pip;
		TH1D* dHist_FCAL_EoverPmeas_pim;
		TH1D* dHist_FCAL_EoverPmeas_pip_PostCorr;
		TH1D* dHist_FCAL_EoverPmeas_pim_PostCorr;
		TH2D* dHist_Delta_Efcal_kinfitE_vs_kinPmag_plus;
		TH2D* dHist_Delta_Efcal_kinfitE_vs_kinPmag_minus;
		TH2D* dHist_Delta_Efcal_kinfitE_vs_kinPmag_plus_PostCorr;
		TH2D* dHist_Delta_Efcal_kinfitE_vs_kinPmag_minus_PostCorr;
		TH2D* dHist_Delta_Efcal_kinfitE_vs_kinPtheta_plus;
		TH2D* dHist_Delta_Efcal_kinfitE_vs_kinPtheta_minus;
		TH1D* dHist_FCAL_Elasticity;
		TH1D* dHist_FCAL_Asymmetry;
		TH1D* dHist_FCAL_Elasticity_PostCorr;
		TH1D* dHist_FCAL_Asymmetry_PostCorr;
		TH1D* dHist_FCAL_Elasticity_Asym_L0pt2;
		TH1D* dHist_FCAL_Elasticity_Asym_G0pt2_L0pt5;
		TH1D* dHist_FCAL_Elasticity_Asym_G0pt5;
		TH1D* dHist_FCAL_Elasticity_Asym_L0pt2_PostCorr;
		TH1D* dHist_FCAL_Elasticity_Asym_G0pt2_L0pt5_PostCorr;
		TH1D* dHist_FCAL_Elasticity_Asym_G0pt5_PostCorr;
		TH1D* dHist_TrackFCAL_DOCA_plus;
		TH1D* dHist_TrackFCAL_DOCA_minus;
		TH1D* dHist_FCAL_E1E9_plus;
		TH1D* dHist_FCAL_E1E9_minus;
		TH1D* dHist_FCAL_E9E25_plus;
		TH1D* dHist_FCAL_E9E25_minus;
		TH1D* dHist_FCAL_kin_res_plus;
		TH1D* dHist_FCAL_kin_res_minus;
		TH1D* dHist_FCAL_kin_res_plus_PostCorr;
		TH1D* dHist_FCAL_kin_res_minus_PostCorr;
		TH1D* dHist_FCAL_meas_res_plus;
		TH1D* dHist_FCAL_meas_res_minus;
		TH1D* dHist_FCAL_meas_res_plus_PostCorr;
		TH1D* dHist_FCAL_meas_res_minus_PostCorr;
		TH2D* dHist_FCAL_Saturation_vs_Eshower_plus;
		TH2D* dHist_FCAL_Saturation_vs_Eshower_minus;
		TH2D* dHist_FCAL_Saturation_vs_Eshower_plus_PostCorr;
		TH2D* dHist_FCAL_Saturation_vs_Eshower_minus_PostCorr;
		TH1D* dHist_FCAL_Saturation_plus;
		TH1D* dHist_FCAL_Saturation_minus;
		TH1D* dHist_FCAL_Saturation_plus_PostCorr;
		TH1D* dHist_FCAL_Saturation_minus_PostCorr;
		TH1D* dHist_FCAL_SumU_plus;
		TH1D* dHist_FCAL_SumU_minus;
		TH1D* dHist_FCAL_SumV_plus;
		TH1D* dHist_FCAL_SumV_minus;
		TH1D* dHist_FCAL_UV_Asymmetry_plus;
		TH1D* dHist_FCAL_UV_Asymmetry_minus;

		enum EChargeView {
			kBoth = 0,
			kPlus,
			kMinus,
			kNumChargeViews
		};

		struct FCALThetaDiagnostics {
			TH2D* Energy[kNumChargeViews] = {nullptr};
			TH2D* EoverP[kNumChargeViews] = {nullptr};
			TH2D* EoverPmeas[kNumChargeViews] = {nullptr};
			TH2D* DeltaEfcal_kinfitE[kNumChargeViews] = {nullptr};
			TH2D* TrackDOCA[kNumChargeViews] = {nullptr};
			TH2D* E1E9[kNumChargeViews] = {nullptr};
			TH2D* E9E25[kNumChargeViews] = {nullptr};
			TH2D* KinRes[kNumChargeViews] = {nullptr};
			TH2D* MeasRes[kNumChargeViews] = {nullptr};
			TH2D* Saturation[kNumChargeViews] = {nullptr};
			TH2D* SumU[kNumChargeViews] = {nullptr};
			TH2D* SumV[kNumChargeViews] = {nullptr};
			TH2D* UVAsymmetry[kNumChargeViews] = {nullptr};

			TH2D* Energy_PostCorr[kNumChargeViews] = {nullptr};
			TH2D* EoverP_PostCorr[kNumChargeViews] = {nullptr};
			TH2D* EoverPmeas_PostCorr[kNumChargeViews] = {nullptr};
			TH2D* DeltaEfcal_kinfitE_PostCorr[kNumChargeViews] = {nullptr};
			TH2D* KinRes_PostCorr[kNumChargeViews] = {nullptr};
			TH2D* MeasRes_PostCorr[kNumChargeViews] = {nullptr};
			TH2D* Saturation_PostCorr[kNumChargeViews] = {nullptr};
		};

		struct FCALPkinDiagnostics {
			TH2D* Energy[kNumChargeViews] = {nullptr};
			TH2D* EoverP[kNumChargeViews] = {nullptr};
			TH2D* EoverPmeas[kNumChargeViews] = {nullptr};
			TH2D* DeltaEfcal_kinfitE[kNumChargeViews] = {nullptr};
			TH2D* TrackDOCA[kNumChargeViews] = {nullptr};
			TH2D* E1E9[kNumChargeViews] = {nullptr};
			TH2D* E9E25[kNumChargeViews] = {nullptr};
			TH2D* KinRes[kNumChargeViews] = {nullptr};
			TH2D* MeasRes[kNumChargeViews] = {nullptr};
			TH2D* Saturation[kNumChargeViews] = {nullptr};
			TH2D* SumU[kNumChargeViews] = {nullptr};
			TH2D* SumV[kNumChargeViews] = {nullptr};
			TH2D* UVAsymmetry[kNumChargeViews] = {nullptr};

			TH2D* Energy_PostCorr[kNumChargeViews] = {nullptr};
			TH2D* EoverP_PostCorr[kNumChargeViews] = {nullptr};
			TH2D* EoverPmeas_PostCorr[kNumChargeViews] = {nullptr};
			TH2D* DeltaEfcal_kinfitE_PostCorr[kNumChargeViews] = {nullptr};
			TH2D* KinRes_PostCorr[kNumChargeViews] = {nullptr};
			TH2D* MeasRes_PostCorr[kNumChargeViews] = {nullptr};
			TH2D* Saturation_PostCorr[kNumChargeViews] = {nullptr};
		};

		struct FCALqvec2Diagnostics {
			TH2D* Energy[kNumChargeViews] = {nullptr};
			TH2D* EoverP[kNumChargeViews] = {nullptr};
			TH2D* EoverPmeas[kNumChargeViews] = {nullptr};
			TH2D* DeltaEfcal_kinfitE[kNumChargeViews] = {nullptr};
			TH2D* TrackDOCA[kNumChargeViews] = {nullptr};
			TH2D* E1E9[kNumChargeViews] = {nullptr};
			TH2D* E9E25[kNumChargeViews] = {nullptr};
			TH2D* KinRes[kNumChargeViews] = {nullptr};
			TH2D* MeasRes[kNumChargeViews] = {nullptr};
			TH2D* Saturation[kNumChargeViews] = {nullptr};
			TH2D* SumU[kNumChargeViews] = {nullptr};
			TH2D* SumV[kNumChargeViews] = {nullptr};
			TH2D* UVAsymmetry[kNumChargeViews] = {nullptr};

			// MC-only post-correction mirrors for observables that depend on FCAL energy correction functions.
			TH2D* Energy_PostCorr[kNumChargeViews] = {nullptr};
			TH2D* EoverP_PostCorr[kNumChargeViews] = {nullptr};
			TH2D* EoverPmeas_PostCorr[kNumChargeViews] = {nullptr};
			TH2D* DeltaEfcal_kinfitE_PostCorr[kNumChargeViews] = {nullptr};
			TH2D* KinRes_PostCorr[kNumChargeViews] = {nullptr};
			TH2D* MeasRes_PostCorr[kNumChargeViews] = {nullptr};
			TH2D* Saturation_PostCorr[kNumChargeViews] = {nullptr};
		};

		FCALThetaDiagnostics dFCALvsTheta;
		FCALPkinDiagnostics dFCALvsPkin;
		FCALqvec2Diagnostics dFCALvsqvec2;

		enum EBestChiSqCategory {
			kNoMVA_cutsBased = 0,
			kPreMVASelection,
			kMLP_ee,
			kMLP_pi,
			kBDT_ee,
			kBDT_pi
		};

		enum EBeamWindow {
			kFullSpectrum = 0,
			kCoherentPeak
		};

		struct CategoryHistogramSet {
			TH1I* BeamEnergy = nullptr;
			TH1I* RelBeamBucket = nullptr;
			TH1I* TaggerAccidentals = nullptr;
			TH1I* MissingMassSquared = nullptr;
			TH1I* MissingEnergy = nullptr;
			TH2D* RecoilThetaVsP = nullptr;
			TH1D* Wepem = nullptr;
			TH1D* qvec2_varWidth = nullptr;
			TH1D* qvec2 = nullptr;
			TH1D* theta1 = nullptr;
			TH1D* theta2 = nullptr;
			TH2D* theta2_vs_theta1 = nullptr;
			TH1D* ep_Pmag = nullptr;
			TH1D* em_Pmag = nullptr;
			TH1D* JTphi[kNumRunPeriods][nPolarizations] = {{nullptr}};
			TH1D* ep_phi[kNumRunPeriods][nPolarizations] = {{nullptr}};
			TH1D* em_phi[kNumRunPeriods][nPolarizations] = {{nullptr}};

			TH1D* FCAL_Energy_pip = nullptr;
			TH1D* FCAL_Energy_pim = nullptr;
			TH1D* FCAL_EoverP_pip = nullptr;
			TH1D* FCAL_EoverP_pim = nullptr;
			TH1D* FCAL_EoverPmeas_pip = nullptr;
			TH1D* FCAL_EoverPmeas_pim = nullptr;
			TH2D* Delta_Efcal_kinfitE_vs_kinPmag_plus = nullptr;
			TH2D* Delta_Efcal_kinfitE_vs_kinPmag_minus = nullptr;
			TH1D* FCAL_Elasticity = nullptr;
			TH1D* FCAL_Asymmetry = nullptr;
			TH1D* FCAL_Elasticity_Asym_L0pt2 = nullptr;
			TH1D* FCAL_Elasticity_Asym_G0pt2_L0pt5 = nullptr;
			TH1D* FCAL_Elasticity_Asym_G0pt5 = nullptr;
			TH1D* TrackFCAL_DOCA_plus = nullptr;
			TH1D* TrackFCAL_DOCA_minus = nullptr;
			TH1D* FCAL_E1E9_plus = nullptr;
			TH1D* FCAL_E1E9_minus = nullptr;
			TH1D* FCAL_E9E25_plus = nullptr;
			TH1D* FCAL_E9E25_minus = nullptr;
			TH1D* FCAL_kin_res_plus = nullptr;
			TH1D* FCAL_kin_res_minus = nullptr;
			TH1D* FCAL_meas_res_plus = nullptr;
			TH1D* FCAL_meas_res_minus = nullptr;
			TH2D* FCAL_Saturation_vs_Eshower_plus = nullptr;
			TH2D* FCAL_Saturation_vs_Eshower_minus = nullptr;
			TH1D* FCAL_Saturation_plus = nullptr;
			TH1D* FCAL_Saturation_minus = nullptr;
			TH1D* FCAL_SumU_plus = nullptr;
			TH1D* FCAL_SumU_minus = nullptr;
			TH1D* FCAL_SumV_plus = nullptr;
			TH1D* FCAL_SumV_minus = nullptr;
			TH1D* FCAL_UV_Asymmetry_plus = nullptr;
			TH1D* FCAL_UV_Asymmetry_minus = nullptr;

			TH1D* qvec2_Meas = nullptr;
			TH1D* theta1_Meas = nullptr;
			TH1D* theta2_Meas = nullptr;
			TH2D* theta2_vs_theta1_Meas = nullptr;
			TH1D* ep_Pmag_Meas = nullptr;
			TH1D* em_Pmag_Meas = nullptr;
			TH1D* JTphi_Meas[kNumRunPeriods][nPolarizations] = {{nullptr}};
			TH1D* ep_phi_Meas[kNumRunPeriods][nPolarizations] = {{nullptr}};
			TH1D* em_phi_Meas[kNumRunPeriods][nPolarizations] = {{nullptr}};
		};

		CategoryHistogramSet dBestChiSqHistograms[kNumBeamWindows][kNumBestChiSqCategories];
		CategoryHistogramSet dRawAccSubdHistograms[kNumBeamWindows][kNumBestChiSqCategories];


		Bool_t OnlyBestChiSqComboInFlat;       // Only save combo with best chiSq to flat tree
		Bool_t ThisComboIsBestChiSq;           // Flag to indicate if the current combo is the best chiSq combo

		// FCAL-specific flags for 2eMissingProton analysis
		Bool_t ApplyFCALEnergyCorrections;      // Apply momentum-dependent FCAL energy corrections
	Bool_t ApplyMLPClassification;		    // Apply MLP-based e+/e- identification
	Bool_t ApplyBDTClassification;		    // Apply BDT-based e+/e- identification

		// Runtime histogram/fill switches for staged debugging (set in Init())
		struct RuntimeFillSwitches {
			bool fillPreFidDiagnostics;     // Pre-fiducial diagnostics block
			bool runPostPreFidBlocks;       // Master switch for blocks after pre-fid fills
			bool fillSystematics;           // FillSystematics() and systematics directories
			bool fillThrownResolution;      // MC thrown/resolution histograms
			bool fillCategoryDirectories;   // CategoryHistogramSet fills (BestChiSq dirs)
			bool fillFCALStudy;             // FCAL study block (standalone FCAL / FCALvsTheta / FCALvsPkin / FCALvsqvec2)
			bool fillNminus1;              // N-1 fiducial validation directories
		};
		RuntimeFillSwitches dRuntimeFillSwitches;

		// FCAL correction functions for e+ and e-
		TF1* ep_FCAL_cor1;
		TF1* ep_FCAL_cor2;
		TF1* em_FCAL_cor1;
		TF1* em_FCAL_cor2;
		TF1* ep_FCAL_cor2data1;
		TF1* ep_FCAL_cor2data2;
		TF1* em_FCAL_cor2data1;
		TF1* em_FCAL_cor2data2;

		// Constants
		Double_t PI = 3.14159;
		Double_t MProton = 0.93827231;  //(GeV/c^2)
		Double_t ElectronMass = 0.000511;

	ClassDef(DSelector_2eMissingProton_Systematics, 0);
};

void DSelector_2eMissingProton_Systematics::Get_ComboWrappers(void)
{
	//Step 0
	dStep0Wrapper = dComboWrapper->Get_ParticleComboStep(0);
	dComboBeamWrapper = static_cast<DBeamParticle*>(dStep0Wrapper->Get_InitialParticle());
	dPositiveWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(0));
	dNegativeWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(1));
	dRecoilWrapper = dStep0Wrapper->Get_FinalParticle(2);
}

#endif // DSelector_2eMissingProton_Systematics_h
