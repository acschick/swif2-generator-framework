#include "DSelector_2piMissingProton_Systematics.h"
#include "TMath.h"

#include <cctype>

// ========== PARTICLE TYPE CONFIGURATION ==========
// Define particle labels for easy porting between channels
// These constants control histogram names and ROOT LaTeX formatting
// 
// TO PORT TO A DIFFERENT CHANNEL:
// Simply change these three definitions and all histograms will update automatically
//
// Current configuration: π+π- channel
static const char* PARTICLE_PLUS = "pip";      // Histogram name suffix for positive track
static const char* PARTICLE_MINUS = "pim";     // Histogram name suffix for negative track  
static const char* PARTICLE_LATEX = "#pi";     // ROOT LaTeX symbol for plot titles

//
// Examples for other channels:
//   Electrons: PARTICLE_PLUS="ep",  PARTICLE_MINUS="em",  PARTICLE_LATEX="e"
//   Kaons:     PARTICLE_PLUS="kp",  PARTICLE_MINUS="km",  PARTICLE_LATEX="K"
//   Muons:     PARTICLE_PLUS="mup", PARTICLE_MINUS="mum", PARTICLE_LATEX="#mu"
//
// This gives you:
//   - Histogram names like "FCAL_Energy_pip" or "FCAL_Energy_ep"
//   - Invariant mass like "InvMass_pipim" or "InvMass_epem"
//   - Axis labels like "E/p π+" or "E/p e+"
// =================================================

// Sanitize function to protect against NaN/Inf values in flat tree
// Replaces non-finite values with sentinel value (default -999)
static inline Float_t sanitize_f(Float_t v, Float_t sentinel = -999.f) {
    return std::isfinite((double)v) ? v : sentinel;
}

static bool EqualsIgnoreCase(const TString& a, const char* b);
static bool IsAutoOrBlank(const TString& value);
static bool IsTokenMode(const TString& value);
static TString ExtractBRunFileTag(TObjArray* tokens);
static TString BuildStructuredTagFromTokens(TObjArray* tokens);

// | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ============================= INIT INIT INIT INIT INIT INIT INIT ==============================
// -----------------------------------------------------------------------------------------------
// | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |
void DSelector_2piMissingProton_Systematics::Init(TTree *locTree)
{
	//SET OUTPUT FILE NAME //can be overriden by user in PROOF
	dOutputFileName = "auto"; // Use auto naming - includes run/file numbers in filename
	dOutputTreeFileName = "";//"ee_MVA_7-20-2020_ComboTrees.root"; //"" for none
	dFlatTreeFileName = ""; //output flat tree (one combo per tree entry), "" for none
	dFlatTreeName = ""; //if blank, default name will be chosen

	// AUTO-GENERATE OUTPUT FILE NAMES FROM INPUT FILE NAME WHEN REQUESTED.
	// "auto" means: use input file stem and replace "tree" -> "hists" (case-insensitive).
	// "tokens" means: derive a compact tokenized tag from full path tokens.
	if(locTree->GetCurrentFile() != NULL) {
		TString inputPath = locTree->GetCurrentFile()->GetName();
		TString baseName = gSystem->BaseName(inputPath);
		TString baseNoExt = baseName;
		if(baseNoExt.EndsWith(".root"))
			baseNoExt.ReplaceAll(".root", "");

		TString autoStem = baseNoExt;
		autoStem.ReplaceAll("tree", "hists");
		autoStem.ReplaceAll("Tree", "hists");
		autoStem.ReplaceAll("TREE", "hists");

		TString trimmedTreeStem = baseNoExt;
		trimmedTreeStem.ReplaceAll("tree", "trimmedtree");
		trimmedTreeStem.ReplaceAll("Tree", "trimmedtree");
		trimmedTreeStem.ReplaceAll("TREE", "trimmedtree");

		TString tokenStem = autoStem;
		const bool wantTokenMode = IsTokenMode(dOutputFileName) || IsTokenMode(dOutputTreeFileName) || IsTokenMode(dFlatTreeFileName);
		if(wantTokenMode) {
			TObjArray* tokens = inputPath.Tokenize("_/.");
			TString structuredTag = BuildStructuredTagFromTokens(tokens);
			TString bRunFileTag = ExtractBRunFileTag(tokens);
			delete tokens;

			if(structuredTag != "")
				tokenStem = structuredTag;
			else if(bRunFileTag != "")
				tokenStem = bRunFileTag;
		}

		if(IsAutoOrBlank(dOutputFileName))
			dOutputFileName = autoStem + ".root";
		else if(IsTokenMode(dOutputFileName))
			dOutputFileName = "Hists_" + tokenStem + ".root";
		if(EqualsIgnoreCase(dOutputTreeFileName, "auto"))
			dOutputTreeFileName = trimmedTreeStem + ".root";
		else if(IsTokenMode(dOutputTreeFileName))
			dOutputTreeFileName = "Tree_" + tokenStem + ".root";
		if(EqualsIgnoreCase(dFlatTreeFileName, "auto"))
			dFlatTreeFileName = autoStem + "_flat.root";
		else if(IsTokenMode(dFlatTreeFileName))
			dFlatTreeFileName = "Flat_" + tokenStem + ".root";
		if(EqualsIgnoreCase(dFlatTreeName, "auto") || dFlatTreeName == "")
			dFlatTreeName = "";

		cout << "Output naming source path: " << inputPath << endl;
		cout << "  Auto stem: " << autoStem << endl;
		if(wantTokenMode)
			cout << "  Token stem: " << tokenStem << endl;
		cout << "  Histogram file: " << dOutputFileName << endl;
		cout << "  Combo tree file: " << dOutputTreeFileName << endl;
		cout << "  Flat tree file: " << dFlatTreeFileName << endl;
		cout << "  Flat tree name: " << dFlatTreeName << endl;
	}
	  
	//Because this function gets called for each TTree in the TChain, we must be careful:
		//We need to re-initialize the tree interface & branch wrappers, but don't want to recreate histograms
	bool locInitializedPriorFlag = dInitializedFlag; //save whether have been initialized previously
	DSelector::Init(locTree); //This must be called to initialize wrappers for each new TTree
	//gDirectory now points to the output file with name dOutputFileName (if any)
	if(locInitializedPriorFlag)
		return; //have already created histograms, etc. below: exit

	Get_ComboWrappers();
	dPreviousRunNumber = 0;
	dIsCPPRunPeriod = false;
	dRunPeriodIndex = -1;

	/*********************************** EXAMPLE USER INITIALIZATION: ANALYSIS ACTIONS **********************************/

	// EXAMPLE: Create deque for histogramming particle masses:
	// // For histogramming the phi mass in phi -> K+ K-
	// // Be sure to change this and dAnalyzeCutActions to match reaction
	//std::deque<Particle_t> MyPhi;
	//MyPhi.push_back(KPlus); MyPhi.push_back(KMinus);

	//ANALYSIS ACTIONS: //Executed in order if added to dAnalysisActions
	//false/true below: use measured/kinfit data

	//PID
	//	dAnalysisActions.push_back(new DHistogramAction_ParticleID(dComboWrapper, false)); //This was DEFAULT setting. I copied from pi file whats below
//	dAnalysisActions.push_back(new DHistogramAction_ParticleID(dComboWrapper, false)); //false: use measured data                                                  
//	dAnalysisActions.push_back(new DHistogramAction_ParticleID(dComboWrapper, true, "KinFit")); //true: use kinfit data 
	//below: value: +/- N ns, Unknown: All PIDs, SYS_NULL: all timing systems
	//dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 0.5, KPlus, SYS_BCAL));

	//MASSES
	//dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, false, Lambda, 1000, 1.0, 1.2, "Lambda"));
	//dAnalysisActions.push_back(new DHistogramAction_MissingMassSquared(dComboWrapper, false, 1000, -0.1, 0.1));

	//KINFIT RESULTS
//	dAnalysisActions.push_back(new DHistogramAction_KinFitResults(dComboWrapper));

	//CUT MISSING MASS
	//dAnalysisActions.push_back(new DCutAction_MissingMassSquared(dComboWrapper, false, -0.03, 0.02));

	//BEAM ENERGY
//	dAnalysisActions.push_back(new DHistogramAction_BeamEnergy(dComboWrapper, false));
	//dAnalysisActions.push_back(new DCutAction_BeamEnergy(dComboWrapper, false, 8.4, 9.05));

	//KINEMATICS
//	dAnalysisActions.push_back(new DHistogramAction_ParticleComboKinematics(dComboWrapper, false));

	// ANALYZE CUT ACTIONS
	 // Change MyPhi to match reaction
	//dAnalyzeCutActions = new DHistogramAction_AnalyzeCutActions( dAnalysisActions, dComboWrapper, false, 0, MyPhi, 1000, 0.9, 2.4, "CutActionEffect" );

	//INITIALIZE ACTIONS
	//If you create any actions that you want to run manually (i.e. don't add to dAnalysisActions), be sure to initialize them here as well
	Initialize_Actions();
	//	dAnalyzeCutActions->Initialize(); // manual action, must call Initialize()

	/************************************** DETERMINE IF ANALYZING SIMULATED DATA *************************************/
	dIsMC = (dTreeInterface->Get_Branch("MCWeight") != NULL);
	if(dIsMC) {
		cout << "Simulated data detected - will book thrown/resolution histograms" << endl;
	} else {
		cout << "Real data - skipping thrown/resolution histograms" << endl;
	}

	/******************** USER INITIALIZATION: STAND-ALONE HISTOGRAMS ********************/





	// Create directory structure for organizing histograms
	TDirectory *dir_Main = gDirectory->mkdir("LeptonPairStudies");
		TDirectory *dir_PreFiducialCutDiagnostics = dir_Main->mkdir("PreFiducialCutDiagnostics");
		TDirectory *dir_FiducialCuts = dir_Main->mkdir("FiducialCuts");
			TDirectory *dir_BestChiSq = dir_FiducialCuts->mkdir("BestChiSq");
				TDirectory *dir_BestChiSq_FullSpectrum = dir_BestChiSq->mkdir("FullSpectrum");
				TDirectory *dir_BestChiSq_CoherentPeak = dir_BestChiSq->mkdir("CoherentPeak");

				const char* bestChiSqCategoryNames[kNumBestChiSqCategories] = {
					"NoMVA_cutsBased", "PreMVASelection", "MLP_ee", "MLP_pi", "BDT_ee", "BDT_pi"
				};

				TDirectory* bestChiSqFidDir[kNumBeamWindows][kNumBestChiSqCategories] = {{nullptr}};
				TDirectory* bestChiSqMeasuredDir[kNumBeamWindows][kNumBestChiSqCategories] = {{nullptr}};
				TDirectory* bestChiSqFCALDir[kNumBeamWindows][kNumBestChiSqCategories] = {{nullptr}};

				TDirectory* bestChiSqBeamDir[kNumBeamWindows] = {dir_BestChiSq_FullSpectrum, dir_BestChiSq_CoherentPeak};
				for(int beamWindowIndex = 0; beamWindowIndex < kNumBeamWindows; ++beamWindowIndex) {
					for(int categoryIndex = 0; categoryIndex < kNumBestChiSqCategories; ++categoryIndex) {
						bestChiSqFidDir[beamWindowIndex][categoryIndex] = bestChiSqBeamDir[beamWindowIndex]->mkdir(bestChiSqCategoryNames[categoryIndex]);
						bestChiSqMeasuredDir[beamWindowIndex][categoryIndex] = bestChiSqFidDir[beamWindowIndex][categoryIndex]->mkdir("Measured");
						bestChiSqFCALDir[beamWindowIndex][categoryIndex] = bestChiSqFidDir[beamWindowIndex][categoryIndex]->mkdir("FCAL");
					}
				}
				
			TDirectory *dir_Thrown = nullptr;
			TDirectory *dir_Thrown_FullSpectrum = nullptr;
			TDirectory *dir_Thrown_CoherentPeak = nullptr;
			if(dIsMC) {
				dir_Thrown = dir_BestChiSq->mkdir("Thrown_and_Resolutions");
				dir_Thrown_FullSpectrum = dir_Thrown->mkdir("FullSpectrum");
				dir_Thrown_CoherentPeak = dir_Thrown->mkdir("CoherentPeak");
			}
			
			TDirectory *dir_FID_Nminus1_Plots = dir_BestChiSq->mkdir("Fid_Nminus1_Plots");
				TDirectory *dir_FID_Nminus1_FullSpectrum = dir_FID_Nminus1_Plots->mkdir("FullSpectrum");
				TDirectory *dir_FID_Nminus1_CoherentPeak = dir_FID_Nminus1_Plots->mkdir("CoherentPeak");

		TDirectory *dir_RawAccSubd = dir_Main->mkdir("RawAccSubd");
			TDirectory *dir_RawAccSubd_FullSpectrum = dir_RawAccSubd->mkdir("FullSpectrum");
				TDirectory *dir_RawAccSubd_CoherentPeak = dir_RawAccSubd->mkdir("CoherentPeak");
				TDirectory* rawAccSubdFidDir[kNumBeamWindows][kNumBestChiSqCategories] = {{nullptr}};
				TDirectory* rawAccBeamDir[kNumBeamWindows] = {dir_RawAccSubd_FullSpectrum, dir_RawAccSubd_CoherentPeak};
				for(int beamWindowIndex = 0; beamWindowIndex < kNumBeamWindows; ++beamWindowIndex) {
					for(int categoryIndex = 0; categoryIndex < kNumBestChiSqCategories; ++categoryIndex) {
						rawAccSubdFidDir[beamWindowIndex][categoryIndex] = rawAccBeamDir[beamWindowIndex]->mkdir(bestChiSqCategoryNames[categoryIndex]);
					}
				}
			
		TDirectory *dir_FCAL = dir_Main->mkdir("FCAL");
				TDirectory *dir_FCALvsTheta = dir_FCAL->mkdir("FCALvsTheta");
					TDirectory *dir_FCALvsTheta_both = dir_FCALvsTheta->mkdir("FCALvsTheta_Both");
					TDirectory *dir_FCALvsTheta_plus = dir_FCALvsTheta->mkdir(Form("FCALvsTheta_%s", PARTICLE_PLUS));
					TDirectory *dir_FCALvsTheta_minus = dir_FCALvsTheta->mkdir(Form("FCALvsTheta_%s", PARTICLE_MINUS));
				TDirectory *dir_FCALvsMomentum = dir_FCAL->mkdir("FCALvsMomentum");
					TDirectory *dir_FCALvsMomentum_both = dir_FCALvsMomentum->mkdir("FCALvsMomentum_Both");
					TDirectory *dir_FCALvsMomentum_plus = dir_FCALvsMomentum->mkdir(Form("FCALvsMomentum_%s", PARTICLE_PLUS));
					TDirectory *dir_FCALvsMomentum_minus = dir_FCALvsMomentum->mkdir(Form("FCALvsMomentum_%s", PARTICLE_MINUS));
				TDirectory *dir_FCALvsqvec2 = dir_FCAL->mkdir("FCALvsqvec2");

		// Centralized systematics directory structure
		dDir_Systematics = dir_Main->mkdir("systematics");
			dDir_Systematics_JTphi = dDir_Systematics->mkdir("JTphi");
			dDir_Systematics_q2 = dDir_Systematics->mkdir("q2");
			dDir_Systematics_theta = dDir_Systematics->mkdir("theta_sys");
			dDir_Systematics_invmass = dDir_Systematics->mkdir("invmass_sys");
			dDir_Systematics_chisq = dDir_Systematics->mkdir("chisq_sys");
			dDir_Systematics_thrown = dDir_Systematics->mkdir("thrown");
		dDir_MVAResponsePlots = dir_Main->mkdir("MVA_response_plots");

	// Variable width bin scheme to match q^2 resolution  
	const Int_t qvec2_NBINS = 200;
  	Double_t qvec2_edges[qvec2_NBINS + 1] = {0.0000, 0.0002, 0.0006, 0.0013, 0.002, 0.003, 0.0042, 0.0055, 0.007, 0.009, 0.011,
    	0.013, 0.015, 0.017, 0.019, 0.021, 0.023, 0.025, 0.027, 0.029, 0.031, 0.033, 0.035, 0.037, 0.039, 0.041, 0.043,
    	0.045, 0.047, 0.049, 0.051, 0.053, 0.055, 0.057, 0.059, 0.061, 0.063, 0.065, 0.067, 0.069, 0.071, 0.073, 0.075,
		0.077, 0.079, 0.081, 0.083, 0.085, 0.087, 0.089, 0.091, 0.093, 0.095, 0.097, 0.099, 0.101, 0.103, 0.105, 0.107,
		0.109, 0.111, 0.113, 0.115, 0.117, 0.119, 0.121, 0.123, 0.125, 0.127, 0.129, 0.131, 0.133, 0.135, 0.137, 0.139,
		0.141, 0.143, 0.145, 0.147, 0.149, 0.151, 0.153, 0.155, 0.157, 0.159, 0.161, 0.163, 0.165, 0.167, 0.169, 0.171,
		0.173, 0.175, 0.177, 0.179, 0.181, 0.183, 0.185, 0.187, 0.189, 0.191, 0.193, 0.195, 0.197, 0.199, 0.201, 0.203,
		0.205, 0.207, 0.209, 0.211, 0.213, 0.215, 0.217, 0.219, 0.221, 0.223, 0.225, 0.227, 0.229, 0.231, 0.233, 0.235,
		0.237, 0.239, 0.241, 0.243, 0.245, 0.247, 0.249, 0.251, 0.253, 0.255, 0.257, 0.259, 0.261, 0.263, 0.265, 0.267,
		0.269, 0.271, 0.273, 0.275, 0.277, 0.279, 0.281, 0.283, 0.285, 0.287, 0.289, 0.291, 0.293, 0.295, 0.297, 0.299,
		0.301, 0.303, 0.305, 0.307, 0.309, 0.311, 0.313, 0.315, 0.317, 0.319, 0.321, 0.323, 0.325, 0.327, 0.329, 0.331,
		0.333, 0.335, 0.337, 0.339, 0.341, 0.343, 0.345, 0.347, 0.349, 0.351, 0.353, 0.355, 0.357, 0.359, 0.361, 0.363,
		0.365, 0.367, 0.369, 0.371, 0.373, 0.375, 0.377, 0.379, 0.381, 0.383, 0.385, 0.387, 0.389, 0.391};

	//================= Begin Diagnostic Histograms ==================
	// These histograms are for checking the distributions of key variables before applying fiducial cuts, MVA selections, or chi^2 selection.
	// Good for showing bi-modal peaks for electrons and pions, or making before and after plots, e.g. splitting the original W_epem histogram into
	// pions and electrons with the neural net. 
	dir_PreFiducialCutDiagnostics->cd();
		dHist_BeamEnergy = new TH1I("BeamEnergy", ";Beam Energy (GeV)", 600, 0.0, 11.4);
		dHist_TaggerAccidentals = new TH1I("TaggerPicketFence", "Vertex time - RF (ns)", 400,-20,20);
		// best chi-sq for the ones directly below
		dHist_BeamEnergy_BestChiSq = new TH1I("BeamEnergy_BestChiSq", ";Beam Energy (GeV)", 300, 0.0, 11.4);
		dHist_Wepem_BestChiSq = new TH1D("Wepem_BestChiSq", Form(";Invariant Mass %s^{+}%s^{-} (GeV/c^{2})", PARTICLE_LATEX, PARTICLE_LATEX), 300, 0.00, 1.4);
		dHist_MissingMassSquared_Measured = new TH1I("pre-MissingMassSquared_Measured_BestChiSq", ";M_{miss}^{2} - M_{target}^{2} ((GeV/c^{2})^{2})", 200, -2.0, 5.0);
		dHist_MissingEnergy_Measured = new TH1I("pre-MissingEnergy_Measured_BestChiSq", ";E_{miss} - E_{target} (GeV)", 200, -2.0, 2.0);
		dHist_preFid_MM2Residual_vs_RecoilP_UpTo1UnusedTrack = new TH2D(
			"pre-MM2Residual_vs_RecoilP_UpTo1UnusedTrack_BestChiSq",
			";|#vec{p}_{p,miss}^{kinfit}| (GeV/c);M_{miss}^{2} - M_{target}^{2} ((GeV/c^{2})^{2})",
			200, 0.0, 2.0, 240, -3.0, 3.0);
		dHist_MissingMassSquared_Measured_AllowOneUnused = new TH1I("pre-MissingMassSquared_Measured_AllowOneUnused_BestChiSq", ";M_{miss}^{2} - M_{target}^{2} ((GeV/c^{2})^{2})", 200, -2.0, 5.0);
		dHist_Wepem = new TH1D(Form("W_%s%s_BestChiSq", PARTICLE_PLUS, PARTICLE_MINUS), 
			Form(";Invariant Mass %s^{+}%s^{-} (GeV/c^{2})", PARTICLE_LATEX, PARTICLE_LATEX), 300, 0.00, 1.4);		
		dHist_KinFitChiSq = new TH1D("pre-KinFitChiSq_BestChiSq", ";Kinematic Fit #chi^{2}", 250, 0.0, 25.0);
		dHist_KinFitCL = new TH1D("pre-KinFitCL_BestChiSq", ";Kinematic Fit Confidence Level", 200, 0.0, 1.0);
		dHist_VertexZ_BestChiSq = new TH1D("pre-VertexZ_BestChiSq", ";Vertex Z (cm)", 100, 40.0, 90.0);

		dHist_EoverP_measured_plus = new TH1D(Form("EoverP_measured_%s_BestChiSq", PARTICLE_PLUS), 
			Form(";E/p %s^{+}", PARTICLE_LATEX), 100, 0.0, 2.0);
		dHist_EoverP_measured_minus = new TH1D(Form("EoverP_measured_%s_BestChiSq", PARTICLE_MINUS), 
			Form(";E/p %s^{-}", PARTICLE_LATEX), 100, 0.0, 2.0);
		dHist_MLPResponse_plus = new TH1D(Form("MLPResponse_%s_BestChiSq", PARTICLE_PLUS), 
			Form(";MLP Response %s^{+}", PARTICLE_LATEX), 100, 0.0, 1.0);
		dHist_MLPResponse_minus = new TH1D(Form("MLPResponse_%s_BestChiSq", PARTICLE_MINUS), 
			Form(";MLP Response %s^{-}", PARTICLE_LATEX), 100, 0.0, 1.0);

		dHist_MLPResponsePlus_vs_MLPResponseMinus = new TH2D("MLPResponsePlus_vs_MLPResponseMinus_BestChiSq", 
			";MLP Response e^{+};MLP Response e^{-}", 100, 0.0, 1.0, 100, 0.0, 1.0);
	//================= End Diagnostic Histograms ==================




	//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>			
    //<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
	//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>




	// ============== Begin Fiducial Analysis Histograms ====================
	// These histograms are for checking the distributions of key variables in different categories of events, e.g. different chi^2 selections, MVA selections, or beam energy windows.

	auto BookCategoryHistogramSet = [&](CategoryHistogramSet& histSet, 
		TDirectory* fiducialDirectory, 
			TDirectory* measuredDirectory, 
			TDirectory* FCALDirectory, 
		bool bookMeasured, bool bookFCAL) {
		fiducialDirectory->cd();
		histSet.BeamEnergy = new TH1I("BeamEnergy", ";Beam Energy (GeV)", 600, 0.0, 11.4);
		histSet.RelBeamBucket = new TH1I("RelBeamBucket", ";Relative Beam Bucket", 21, -10.5, 10.5);
		histSet.TaggerAccidentals = new TH1I("TaggerAccidentals", "Vertex time - RF (ns)", 400,-20,20);
		histSet.MissingMassSquared = new TH1I("MissingMassSquared", ";M_{miss}^{2} - M_{target}^{2} ((GeV/c^{2})^{2})", 200, -2.00, 5.00);
		histSet.MissingEnergy = new TH1I("MissingEnergy", ";E_{miss} - E_{target} (GeV)", 200, -2.00, 2.00);
		histSet.RecoilThetaVsP = new TH2D("RecoilThetaVsP", ";|#vec{p}_{p,miss}^{kinfit}| (GeV/c);#theta_{p,miss}^{kinfit} (deg)", 120, 0.0, 0.30, 120, 0.0, 100.0);
		histSet.Wepem = new TH1D(Form("W_%s%s", PARTICLE_PLUS, PARTICLE_MINUS), Form(";Invariant Mass %s^{+}%s^{-} (GeV/c^{2})", PARTICLE_LATEX, PARTICLE_LATEX), 300, 0.0, 1.2);
		histSet.qvec2_varWidth = new TH1D("qvec2_varWidth", ";#vec{q}^{2} (GeV/c)^{2}", qvec2_NBINS, qvec2_edges);
		histSet.qvec2 = new TH1D("qvec2", ";#vec{q}^{2} (GeV/c)^{2}", 30, 0.000000, 0.015);
		histSet.theta1 = new TH1D("theta1", ";lab #theta_{1} (deg)", 100, 0, 15);
		histSet.theta2 = new TH1D("theta2", ";lab #theta_{2} (deg)", 100, 0, 15);
		histSet.theta2_vs_theta1 = new TH2D("theta2_vs_theta1", ";lab #theta_{1} (deg);lab #theta_{2} (deg)", 100, 0, 15, 100, 0, 15);
		histSet.ep_Pmag = new TH1D("ep_Pmag", ";|#vec{p}_{e^{+}}| (GeV/c)", 100, 0, 10);
		histSet.em_Pmag = new TH1D("em_Pmag", ";|#vec{p}_{e^{-}}| (GeV/c)", 100, 0, 10);
		for(int runPeriodIndex = 0; runPeriodIndex < kNumRunPeriods; ++runPeriodIndex) {
			for(int polIndex = 0; polIndex < nPolarizations; ++polIndex) {
				histSet.JTphi[runPeriodIndex][polIndex] = new TH1D(
					Form("JTphi_%s_%s", RunPeriodTags[runPeriodIndex], Polarizations[polIndex]),
					";J_{T}.#phi (deg)", 90, -180, 180);
				histSet.ep_phi[runPeriodIndex][polIndex] = new TH1D(
					Form("ep_phi_%s_%s", RunPeriodTags[runPeriodIndex], Polarizations[polIndex]),
					Form(";#phi_{e^{+}} (deg) run=%s #theta_{pol}=%s#circ", RunPeriodTags[runPeriodIndex], Polarizations[polIndex]),
					180, -180, 180);
				histSet.em_phi[runPeriodIndex][polIndex] = new TH1D(
					Form("em_phi_%s_%s", RunPeriodTags[runPeriodIndex], Polarizations[polIndex]),
					Form(";#phi\t_{e^{-}} (deg) run=%s #theta_{pol}=%s#circ", RunPeriodTags[runPeriodIndex], Polarizations[polIndex]),
					180, -180, 180);
			}
		}

		if(bookFCAL && FCALDirectory != nullptr) {
			FCALDirectory->cd();
			histSet.FCAL_Energy_pip = new TH1D(Form("FCAL_Energy_%s", PARTICLE_PLUS), ";GeV", 200, 0, 10);
			histSet.FCAL_Energy_pim = new TH1D(Form("FCAL_Energy_%s", PARTICLE_MINUS), ";GeV", 200, 0, 10);
			histSet.FCAL_EoverP_pip = new TH1D(Form("FCAL_EoverP_%s", PARTICLE_PLUS), Form(";E/p %s^{+}", PARTICLE_LATEX), 200, 0.0, 1.5);
			histSet.FCAL_EoverP_pim = new TH1D(Form("FCAL_EoverP_%s", PARTICLE_MINUS), Form(";E/p %s^{-}", PARTICLE_LATEX), 200, 0.0, 1.5);
			histSet.FCAL_EoverPmeas_pip = new TH1D(Form("FCAL_EoverPmeas_%s", PARTICLE_PLUS), Form(";E_{FCAL}/P_{meas} %s^{+}", PARTICLE_LATEX), 200, 0.0, 1.5);
			histSet.FCAL_EoverPmeas_pim = new TH1D(Form("FCAL_EoverPmeas_%s", PARTICLE_MINUS), Form(";E_{FCAL}/P_{meas} %s^{-}", PARTICLE_LATEX), 200, 0.0, 1.5);
			histSet.Delta_Efcal_kinfitE_vs_kinPmag_plus = new TH2D(Form("Delta_Efcal-kinfitE_vs_kinPmag_%s", PARTICLE_PLUS), Form(";p_{kinfit} %s^{+};E_{FCAL} - E_{kinfit} %s^{+}", PARTICLE_LATEX, PARTICLE_LATEX), 200, 0.0, 11.0, 200, -2.0, 2.0);
			histSet.Delta_Efcal_kinfitE_vs_kinPmag_minus = new TH2D(Form("Delta_Efcal-kinfitE_vs_kinPmag_%s", PARTICLE_MINUS), Form(";p_{kinfit} %s^{-};E_{FCAL} - E_{kinfit} %s^{-}", PARTICLE_LATEX, PARTICLE_LATEX), 200, 0.0, 11.0, 200, -2.0, 2.0);
			histSet.FCAL_Elasticity = new TH1D("FCAL_Elasticity", ";(E1 + E2)/Ebeam", 200, 0.0, 1.4);
			histSet.FCAL_Asymmetry = new TH1D("FCAL_Asymmetry", ";|(E1 - E2)/(E1 + E2)|", 200, 0.0, 1.0);
			histSet.FCAL_Elasticity_Asym_L0pt2 = new TH1D("FCAL_Elasticity_Asym_L0pt2", ";(E1 + E2)/Ebeam", 200, 0.0, 1.4);
			histSet.FCAL_Elasticity_Asym_G0pt2_L0pt5 = new TH1D("FCAL_Elasticity_Asym_G0pt2_L0pt5", ";(E1 + E2)/Ebeam", 200, 0.0, 1.4);
			histSet.FCAL_Elasticity_Asym_G0pt5 = new TH1D("FCAL_Elasticity_Asym_G0pt5", ";(E1 + E2)/Ebeam", 200, 0.0, 1.4);
			histSet.TrackFCAL_DOCA_plus = new TH1D(Form("FCAL_DOCA_%s", PARTICLE_PLUS), Form(";FCAL_DOCA_%s", PARTICLE_PLUS), 200, 0, 4);
			histSet.TrackFCAL_DOCA_minus = new TH1D(Form("FCAL_DOCA_%s", PARTICLE_MINUS), Form(";FCAL_DOCA_%s", PARTICLE_MINUS), 200, 0, 4);
			histSet.FCAL_E1E9_plus = new TH1D(Form("FCAL_E1E9_%s", PARTICLE_PLUS), Form(";(Energy in central cell)/(Energy in 3x3) %s^{+}", PARTICLE_LATEX), 200, 0.0, 1.1);
			histSet.FCAL_E1E9_minus = new TH1D(Form("FCAL_E1E9_%s", PARTICLE_MINUS), Form(";(Energy in central cell)/(Energy in 3x3) %s^{-}", PARTICLE_LATEX), 200, 0.0, 1.1);
			histSet.FCAL_E9E25_plus = new TH1D(Form("FCAL_E9E25_%s", PARTICLE_PLUS), Form(";(Energy in 3x3)/(Energy in 5x5) %s^{+}", PARTICLE_LATEX), 200, 0.0, 1.1);
			histSet.FCAL_E9E25_minus = new TH1D(Form("FCAL_E9E25_%s", PARTICLE_MINUS), Form(";(Energy in 3x3)/(Energy in 5x5) %s^{-}", PARTICLE_LATEX), 200, 0.0, 1.1);
			histSet.FCAL_kin_res_plus = new TH1D(Form("FCAL_kin_res_%s", PARTICLE_PLUS), Form(";|(E_{FCAL} - E_{kinfit})|/E_{kinfit} %s^{+}", PARTICLE_LATEX), 200, 0.0, 1.0);
			histSet.FCAL_kin_res_minus = new TH1D(Form("FCAL_kin_res_%s", PARTICLE_MINUS), Form(";|(E_{FCAL} - E_{kinfit})|/E_{kinfit} %s^{-}", PARTICLE_LATEX), 200, 0.0, 1.0);
			histSet.FCAL_meas_res_plus = new TH1D(Form("FCAL_meas_res_%s", PARTICLE_PLUS), Form(";|(E_{FCAL} - E_{meas})|/E_{meas} %s^{+}", PARTICLE_LATEX), 200, 0.0, 1.0);
			histSet.FCAL_meas_res_minus = new TH1D(Form("FCAL_meas_res_%s", PARTICLE_MINUS), Form(";|(E_{FCAL} - E_{meas})|/E_{meas} %s^{-}", PARTICLE_LATEX), 200, 0.0, 1.0);
			histSet.FCAL_Saturation_vs_Eshower_plus = new TH2D(Form("FCAL_saturation_vs_Eshower_%s", PARTICLE_PLUS), ";E shower ;Eshower*E1E9*E9E25", 100, 0, 10, 200, 0, 14);
			histSet.FCAL_Saturation_vs_Eshower_minus = new TH2D(Form("FCAL_saturation_vs_Eshower_%s", PARTICLE_MINUS), ";E shower ;Eshower*E1E9*E9E25", 100, 0, 10, 200, 0, 14);
			histSet.FCAL_Saturation_plus = new TH1D(Form("FCAL_saturation_%s", PARTICLE_PLUS), ";Eshower*E1E9*E9E25", 200, 0, 14);
			histSet.FCAL_Saturation_minus = new TH1D(Form("FCAL_saturation_%s", PARTICLE_MINUS), ";Eshower*E1E9*E9E25", 200, 0, 14);
			histSet.FCAL_SumU_plus = new TH1D(Form("FCAL_SumU_%s", PARTICLE_PLUS), Form(";SumU %s^{+}", PARTICLE_LATEX), 200, 0., 14);
			histSet.FCAL_SumU_minus = new TH1D(Form("FCAL_SumU_%s", PARTICLE_MINUS), Form(";SumU %s^{-}", PARTICLE_LATEX), 200, 0., 14);
			histSet.FCAL_SumV_plus = new TH1D(Form("FCAL_SumV_%s", PARTICLE_PLUS), Form(";SumV %s^{+}", PARTICLE_LATEX), 200, 0., 20);
			histSet.FCAL_SumV_minus = new TH1D(Form("FCAL_SumV_%s", PARTICLE_MINUS), Form(";SumV %s^{-}", PARTICLE_LATEX), 200, 0., 20);
			histSet.FCAL_UV_Asymmetry_plus = new TH1D(Form("FCAL_UV_Asymmetry_%s", PARTICLE_PLUS), Form(";A_{uv} = |#sigma^{2}_{u} - #sigma^{2}_{v}|/|#sigma^{2}_{u} + #sigma^{2}_{v}| %s^{+}", PARTICLE_LATEX), 200, 0., 1.0);
			histSet.FCAL_UV_Asymmetry_minus = new TH1D(Form("FCAL_UV_Asymmetry_%s", PARTICLE_MINUS), Form(";A_{uv} = |#sigma^{2}_{u} - #sigma^{2}_{v}|/|#sigma^{2}_{u} + #sigma^{2}_{v}| %s^{-}", PARTICLE_LATEX), 200, 0., 1.0);
		}

		if(bookMeasured && measuredDirectory != nullptr) {
			measuredDirectory->cd();
			histSet.qvec2_Meas = new TH1D("qvec2", ";#vec{q}^{2} (GeV/c)^{2}", qvec2_NBINS, qvec2_edges);
			histSet.theta1_Meas = new TH1D("theta1", ";lab #theta_{1} (deg)", 100, 0, 15);
			histSet.theta2_Meas = new TH1D("theta2", ";lab #theta_{2} (deg)", 100, 0, 15);
			histSet.theta2_vs_theta1_Meas = new TH2D("theta2_vs_theta1", ";lab #theta_{1} (deg);lab #theta_{2} (deg)", 100, 0, 15, 100, 0, 15);
			histSet.ep_Pmag_Meas = new TH1D("ep_Pmag", ";|#vec{p}_{e^{+}}| (GeV/c)", 100, 0, 10);
			histSet.em_Pmag_Meas = new TH1D("em_Pmag", ";|#vec{p}_{e^{-}}| (GeV/c)", 100, 0, 10);
			for(int runPeriodIndex = 0; runPeriodIndex < kNumRunPeriods; ++runPeriodIndex) {
				for(int polIndex = 0; polIndex < nPolarizations; ++polIndex) {
					histSet.JTphi_Meas[runPeriodIndex][polIndex] = new TH1D(
						Form("JTphi_%s_%s", RunPeriodTags[runPeriodIndex], Polarizations[polIndex]),
						";J_{T}.#phi (deg)", 90, -180, 180);
					histSet.ep_phi_Meas[runPeriodIndex][polIndex] = new TH1D(
						Form("ep_phi_%s_%s", RunPeriodTags[runPeriodIndex], Polarizations[polIndex]),
						Form(";#phi_{e^{+}} (deg) run=%s #theta_{pol}=%s#circ", RunPeriodTags[runPeriodIndex], Polarizations[polIndex]),
						180, -180, 180);
					histSet.em_phi_Meas[runPeriodIndex][polIndex] = new TH1D(
						Form("em_phi_%s_%s", RunPeriodTags[runPeriodIndex], Polarizations[polIndex]),
						Form(";#phi\t_{e^{-}} (deg) run=%s #theta_{pol}=%s#circ", RunPeriodTags[runPeriodIndex], Polarizations[polIndex]),
						180, -180, 180);
				}
			}
		}
	};

	for(int beamWindowIndex = 0; beamWindowIndex < kNumBeamWindows; ++beamWindowIndex) {
		for(int categoryIndex = 0; categoryIndex < kNumBestChiSqCategories; ++categoryIndex) {
			BookCategoryHistogramSet(
				dBestChiSqHistograms[beamWindowIndex][categoryIndex],
				bestChiSqFidDir[beamWindowIndex][categoryIndex],
				bestChiSqMeasuredDir[beamWindowIndex][categoryIndex],
				bestChiSqFCALDir[beamWindowIndex][categoryIndex],
				true,
				true
			);

			BookCategoryHistogramSet(
				dRawAccSubdHistograms[beamWindowIndex][categoryIndex],
				rawAccSubdFidDir[beamWindowIndex][categoryIndex],
				nullptr,
				nullptr,
				false,
				false
			);
		}
	}

	// ============== End Fiducial Analysis Histograms ====================





	//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>			
    //<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
	//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>



	
	//                       THROWN AND RESOLUTION HISTOGRAMS
	//                  
	//            ------->           ------->          ---->
	//                o
	//            O /                    o                  q
	//           /|                                      o >-|.
	//           / \'                                   ___  m>
	/*========================================================================================*/
	if(dIsMC) {
		const char* beamWindowNames[kNumBeamWindows] = {"FullSpectrum", "CoherentPeak"};
		TDirectory* thrownBeamDirs[kNumBeamWindows] = {dir_Thrown_FullSpectrum, dir_Thrown_CoherentPeak};
		
		for(int beamWindowIndex = 0; beamWindowIndex < kNumBeamWindows; ++beamWindowIndex) {
			thrownBeamDirs[beamWindowIndex]->cd();
			const char* beamSuffix = beamWindowNames[beamWindowIndex];

			dHist_RecoilThetaVsP_Thrown[beamWindowIndex] = new TH2D(
				Form("RecoilThetaVsP_Thrown_%s", beamSuffix),
				";|#vec{p}_{p}^{thrown}| (GeV/c);#theta_{p}^{thrown} (deg)",
				120, 0.0, 0.25, 120, 0.0, 100.0);

			dHist_qvec2_varWidth_Thrown[beamWindowIndex] = new TH1D(Form("qvec2_varWidth_Thrown_%s", beamSuffix), ";#vec{q}^{2} (GeV/c)^{2}", qvec2_NBINS, qvec2_edges);
			dHist_qvec2_Thrown[beamWindowIndex] = new TH1D(Form("qvec2_Thrown_%s", beamSuffix), ";#vec{q}^{2} (GeV/c)^{2}", 30, 0.000000, 0.015);
			/*<o><o><o><o><o><o><o><o> FOR FORM FACTOR STUDY: DETERMINE VARIABLE WIDTH BIN SCHEME ACCORDING TO RESOLUTION <o><o><o><o><o><o><o><o>*/
			/*<o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o>*/
			dHist_qvec2_res_vs_q2kin[beamWindowIndex] = new TH2D(Form("q2kinRes_vs_q2kin_%s", beamSuffix), ";#vec{q}^{2}_{kin} (GeV/c)^{2};#vec{q}^{2}_{kin} - #vec{q}^{2}_{thrown} (GeV/c)^{2}", 200, 0.00, 0.04, 200, -0.0015, 0.0015);
			/*<o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o><o>*/
			/******************************************************************************************************************************************/

			dHist_theta_KinRes_vs_theta_Thrown[beamWindowIndex] = new TH2D(Form("thetaKinRes_vs_thetaThrown_%s", beamSuffix), ";e #theta_{thrown} (deg);#theta_{kin} - #theta_{thrown} (deg)", 200, 0.0, 15.0, 200, -1.0, 1.0);
			dHist_theta1_Thrown[beamWindowIndex] = new TH1D(Form("theta1_Thrown_%s", beamSuffix), ";lab #theta_{1} (deg)", 100, 0, 15);
			dHist_theta2_Thrown[beamWindowIndex] = new TH1D(Form("theta2_Thrown_%s", beamSuffix), ";lab #theta_{2} (deg)", 100, 0, 15);

			dHist_ep_Pmag_Thrown[beamWindowIndex] = new TH1D(Form("ep_Pmag_Thrown_%s", beamSuffix), ";|#vec{p}_{e^{+}}| (GeV/c)", 100, 0, 10);
			dHist_em_Pmag_Thrown[beamWindowIndex] = new TH1D(Form("em_Pmag_Thrown_%s", beamSuffix), ";|#vec{p}_{e^{-}}| (GeV/c)", 100, 0, 10);
			dHist_ep_Pmag_KinRes[beamWindowIndex] = new TH1D(Form("ep_Pmag_KinRes_%s", beamSuffix), ";|#vec{p}_{e^{+}}| (GeV/c)", 100, -0.05, 0.05);
			dHist_em_Pmag_KinRes[beamWindowIndex] = new TH1D(Form("em_Pmag_KinRes_%s", beamSuffix), ";|#vec{p}_{e^{-}}| (GeV/c)", 100, -0.05, 0.05);
			dHist_ep_Pmag_KinRes_vs_ep_Pmag_Thrown[beamWindowIndex] = new TH2D(Form("ep_Pmag_KinRes_vs_ep_Pmag_Thrown_%s", beamSuffix), ";|#vec{p}_{e^{+}}|_{thrown} (GeV/c);|#vec{p}_{e^{+}}|_{kin} - |#vec{p}_{e^{+}}|_{thrown} (GeV/c)", 100, 0, 10, 100, -0.05, 0.05);
			dHist_em_Pmag_KinRes_vs_em_Pmag_Thrown[beamWindowIndex] = new TH2D(Form("em_Pmag_KinRes_vs_em_Pmag_Thrown_%s", beamSuffix), ";|#vec{p}_{e^{-}}|_{thrown} (GeV/c);|#vec{p}_{e^{-}}|_{kin} - |#vec{p}_{e^{-}}|_{thrown} (GeV/c)", 100, 0, 10, 100, -0.05, 0.05);

			dHist_Wepem_Thrown[beamWindowIndex] = new TH1D(Form("Wepem_Thrown_%s", beamSuffix), ";Invariant Mass %s^{+}%s^{-} (GeV/c^{2})", 300, 0.0, 1.0);
			dHist_Wepem_KinRes[beamWindowIndex] = new TH1D(Form("Wepem_KinRes_%s", beamSuffix), ";Invariant Mass %s^{+}%s^{-} (GeV/c^{2})", 300, 0.0, 1.0);
			dHist_Wepem_KinRes_vs_Wepem_Thrown[beamWindowIndex] = new TH2D(Form("Wepem_KinRes_vs_Wepem_Thrown_%s", beamSuffix), ";Invariant Mass %s^{+}%s^{-} thrown (GeV/c^{2});Invariant Mass %s^{+}%s^{-} kin - thrown (GeV/c^{2})", 300, 0.00, 1.0, 100, -0.05, 0.05);

			for(int runPeriodIndex = 0; runPeriodIndex < kNumRunPeriods; ++runPeriodIndex) {
				for(int i = 0; i < nPolarizations; i++){
					dHist_JTphi_Thrown[beamWindowIndex][runPeriodIndex][i] = new TH1D(
						Form("JTphi_Thrown_%s_%s_%s", beamSuffix, RunPeriodTags[runPeriodIndex], Polarizations[i]),
						";#vec{J}_{T}.#phi (deg)", 90, -180, 180);
					dHist_JTphi_kinRes[beamWindowIndex][runPeriodIndex][i] = new TH1D(
						Form("JTphi_KinRes_%s_%s_%s", beamSuffix, RunPeriodTags[runPeriodIndex], Polarizations[i]),
						";#vec{J}_{T}.#phi (deg)", 90, -180, 180);
		
					dHist_ep_phi_Thrown[beamWindowIndex][runPeriodIndex][i] = new TH1D(
						Form("ep_phi_Thrown_%s_%s_%s", beamSuffix, RunPeriodTags[runPeriodIndex], Polarizations[i]),
						Form(";e^{+} #phi (deg), run=%s #theta_{pol}=%s#circ", RunPeriodTags[runPeriodIndex], Polarizations[i]),
						180, -180, 180);
					dHist_em_phi_Thrown[beamWindowIndex][runPeriodIndex][i] = new TH1D(
						Form("em_phi_Thrown_%s_%s_%s", beamSuffix, RunPeriodTags[runPeriodIndex], Polarizations[i]),
						Form(";e^{-} #phi (deg), run=%s #theta_{pol}=%s#circ", RunPeriodTags[runPeriodIndex], Polarizations[i]),
						180, -180, 180);
				}
			}
			dHist_phi_KinRes[beamWindowIndex] = new TH1D(Form("phi_KinRes_%s", beamSuffix), ";e^{#pm} #phi (deg)", 180, -180, 180);
		} // end beam window loop

	} // end if(dIsMC)

	// --------------------------> END THROWN HISTOGRAMS --------------------------->
	



	//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>			
    //<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
	//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>



	//N-1 N-1 N-1 N-1 N-1 N-1 N-1 N-1 N-1 N-1 N-1 N-1 N-1 N-1 N-1 N-1 N-1 N-1 
	// 					FIDUCIAL CUT N-1 HISTOGRAMS
	//N-1 N-1 N-1 N-1 N-1 N-1 N-1 N-1 N-1 N-1 N-1 N-1 N-1 N-1 N-1 N-1 N-1 N-1 
	
	const char* beamWindowNames[kNumBeamWindows] = {"FullSpectrum", "CoherentPeak"};
	TDirectory* n1BeamDirs[kNumBeamWindows] = {dir_FID_Nminus1_FullSpectrum, dir_FID_Nminus1_CoherentPeak};
	
	for(int beamWindowIndex = 0; beamWindowIndex < kNumBeamWindows; ++beamWindowIndex) {
		n1BeamDirs[beamWindowIndex]->cd();
		const char* beamSuffix = beamWindowNames[beamWindowIndex];
		
		// N-1 plots. Nminus1 followed by the variable name implies all fiducial cuts are applied except the one in the variable name.
		dHist_FID_Nminus1_BeamEnergy[beamWindowIndex] = new TH1I(Form("FID_Nminus1_BeamEnergy_%s", beamSuffix), ";Beam Energy (GeV)", 600, 0.0, 11.4);
		// No cut on Invariant Mass
		dHist_FID_Nminus1_Wepem[beamWindowIndex] = new TH1D(Form("FID_Nminus1_Wepem_noWCuts_%s", beamSuffix), Form(";Invariant Mass %s^{+}%s^{-} (GeV/c^{2})", PARTICLE_LATEX, PARTICLE_LATEX), 300, 0.0, 1.2);
		// Momentum distributions with no momentum cuts
		dHist_FID_Nminus1_p_ep[beamWindowIndex] = new TH1D(Form("FID_Nminus1_p_ep_%s", beamSuffix), ";|#vec{p}_{e^{+}}| (GeV/c)", 100, 0, 10);
		dHist_FID_Nminus1_p_em[beamWindowIndex] = new TH1D(Form("FID_Nminus1_p_em_%s", beamSuffix), ";|#vec{p}_{e^{-}}| (GeV/c)", 100, 0, 10);
		
		// theta plots: standard N - 1 
		dHist_FID_Nminus1_theta1[beamWindowIndex] = new TH1D(Form("FID_Nminus1_theta1_%s", beamSuffix), ";lab #theta_{e^{+}} (deg)", 100, 0, 15);
		dHist_FID_Nminus1_theta2[beamWindowIndex] = new TH1D(Form("FID_Nminus1_theta2_%s", beamSuffix), ";lab #theta_{e^{-}} (deg)", 100, 0, 15);
		dHist_FID_Nminus1_theta2_vs_theta1[beamWindowIndex] = new TH2D(Form("FID_Nminus1_theta2_vs_theta1_%s", beamSuffix), ";lab #theta_{1} (deg);lab #theta_{2} (deg)", 100, 0, 15, 100, 0, 15);
		dHist_FID_Nminus1_VertexZ[beamWindowIndex] = new TH1D(Form("FID_Nminus1_VertexZ_noVertexZCut_allOthers_%s", beamSuffix), ";Vertex Z (cm)", 100, 40.0, 90.0);
		// theta plots: no theta cuts whatsoever
		dHist_FID_Nminus1_theta1_noThetaCuts[beamWindowIndex] = new TH1D(Form("FID_Nminus1_theta1_noThetaCuts_%s", beamSuffix), ";lab #theta_{e^{+}} (deg)", 100, 0, 15);
		dHist_FID_Nminus1_theta2_noThetaCuts[beamWindowIndex] = new TH1D(Form("FID_Nminus1_theta2_noThetaCuts_%s", beamSuffix), ";lab #theta_{e^{-}} (deg)", 100, 0, 15);
		dHist_FID_Nminus1_theta2_vs_theta1_noThetaCuts[beamWindowIndex] = new TH2D(Form("FID_Nminus1_theta2_vs_theta1_noThetaCuts_%s", beamSuffix), ";lab #theta_{1} (deg);lab #theta_{2} (deg)", 100, 0,	15, 100, 0, 15);

		// Measured E/p Preselection Cut 
		dHist_FID_Nminus1_EoverP_ep[beamWindowIndex] = new TH1D(Form("FID_Nminus1_EoverP_ep_noEoPCut_allOthers_%s", beamSuffix), ";E/p_{e^{+}}", 100, 0, 1.3);
		dHist_FID_Nminus1_EoverP_em[beamWindowIndex] = new TH1D(Form("FID_Nminus1_EoverP_em_noEoPCut_allOthers_%s", beamSuffix), ";E/p_{e^{-}}", 100, 0, 1.3);

		// Exclusivity cuts
		dHist_FID_Nminus1_NumUnusedTracks[beamWindowIndex] = new TH1I(Form("FID_Nminus1_NumUnusedTracks_noUnusedTracksCut_allOthers_%s", beamSuffix), ";Number of Unused Tracks", 10, 0, 10);
		dHist_FID_Nminus1_UnusedShowerEnergy[beamWindowIndex] = new TH1D(Form("FID_Nminus1_UnusedShowerEnergy_noUnusedShowersCut_allOthers_%s", beamSuffix), ";Unused Shower Energy (GeV)", 100, 0, 10);

		// Forward detector
		dHist_FID_Nminus1_FCALEnergy_ep[beamWindowIndex] = new TH1D(Form("FID_Nminus1_FCALEnergy_ep_noFCALEnergyCut_allOthers_%s", beamSuffix), ";FCAL Energy e^{+} (GeV)", 100, 0, 10);
		dHist_FID_Nminus1_FCALEnergy_em[beamWindowIndex] = new TH1D(Form("FID_Nminus1_FCALEnergy_em_noFCALEnergyCut_allOthers_%s", beamSuffix), ";FCAL Energy e^{-} (GeV)", 100, 0, 10);
		dHist_FID_Nminus1_TOFdEdx_ep[beamWindowIndex] = new TH1D(Form("FID_Nminus1_TOFdEdx_ep_noTOFdEdxCut_allOthers_%s", beamSuffix), ";TOF dE/dx e^{+} (keV/cm)", 100, 0, 10);
		dHist_FID_Nminus1_TOFdEdx_em[beamWindowIndex] = new TH1D(Form("FID_Nminus1_TOFdEdx_em_noTOFdEdxCut_allOthers_%s", beamSuffix), ";TOF dE/dx e^{-} (keV/cm)", 100, 0, 10);
		
		// MLP output
		dHist_FID_Nminus1_MLP_ep[beamWindowIndex] = new TH1D(Form("FID_Nminus1_MLP_ep_noMLPCut_allOthers_%s", beamSuffix), ";MLP Output e^{+}", 100, 0, 1);
		dHist_FID_Nminus1_MLP_em[beamWindowIndex] = new TH1D(Form("FID_Nminus1_MLP_em_noMLPCut_allOthers_%s", beamSuffix), ";MLP Output e^{-}", 100, 0, 1);
		dHist_FID_Nminus1_MLP_ep_vs_em[beamWindowIndex] = new TH2D(Form("FID_Nminus1_MLP_ep_vs_em_noMLPCut_allOthers_%s", beamSuffix), ";MLP Output e^{+};MLP Output e^{-}", 100, 0, 1, 100, 0, 1);	

		//BDT output
		dHist_FID_Nminus1_BDT_ep[beamWindowIndex] = new TH1D(Form("FID_Nminus1_BDT_ep_noBDTCut_allOthers_%s", beamSuffix), ";BDT Output e^{+}", 100, -0.5, 0.5);
		dHist_FID_Nminus1_BDT_em[beamWindowIndex] = new TH1D(Form("FID_Nminus1_BDT_em_noBDTCut_allOthers_%s", beamSuffix), ";BDT Output e^{-}", 100, -0.5, 0.5);	
		dHist_FID_Nminus1_BDT_ep_vs_em[beamWindowIndex] = new TH2D(Form("FID_Nminus1_BDT_ep_vs_em_noBDTCut_allOthers_%s", beamSuffix), ";BDT Output e^{+};BDT Output e^{-}", 100, -0.5, 0.5, 100, -0.5, 0.5);

		// Special No MVA, cuts based analysis E/p and Elasticity N-1 plots
		dHist_FID_Nminus1_EoverP_ep_and_noMVA[beamWindowIndex] = new TH1D(Form("FID_Nminus1_EoverP_ep_noEoP1Cut_noMVA_allOthers_%s", beamSuffix), ";E/p_{e^{+}}", 100, 0, 1.3);
		dHist_FID_Nminus1_FCALElasticity_and_noMVA[beamWindowIndex] = new TH1D(Form("FID_Nminus1_FCALElasticity_noElasticityCut_noMVA_allOthers_%s", beamSuffix), ";(E1 + E2)/Ebeam", 100, 0, 1.4);
	} // end beam window loop

	//N-1 N-1 N-1 N-1 N-1 N-1 N-1 N-1 N-1 N-1 N-1 N-1 N-1 N-1 N-1 N-1 N-1 N-1 
	// 					END FIDUCIAL CUT N-1 HISTOGRAMS
	//N-1 N-1 N-1 N-1 N-1 N-1 N-1 N-1 N-1 N-1 N-1 N-1 N-1 N-1 N-1 N-1 N-1 N-1



	// ^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|
	// <><><><><><><><><><><><><><>><><><><><><><><><><><><><><
	// <><><><><><><><><><><><><><>><><><><><><><><><><><><><><
	/*==========================================================================*/
			// +----+----+----+
			// |    |    |    |
			// +----+----+----+
			// |    |  O |    | <~~~~~~~~~~~~~ FCAL STUDY ~~~~~~~~~~~~~~~~
			// +----+----+----+
			// |    |    |    |
			// +----+----+----+
	// --------------------------------------------------------------------------		
	// Need to determine the range in theta, momentum and q^2 that 
	// FCAL shower variables match between simulation and data.
	/*==========================================================================*/
	dir_FCAL->cd();
		dHist_FCAL_Energy_pip = new TH1D(Form("FCAL_Energy_%s", PARTICLE_PLUS), ";GeV", 200, 0, 10);
		dHist_FCAL_Energy_pim = new TH1D(Form("FCAL_Energy_%s", PARTICLE_MINUS), ";GeV", 200, 0, 10);
		dHist_FCAL_EoverP_pip = new TH1D(Form("FCAL_EoverP_%s", PARTICLE_PLUS), 
			Form(";E/p %s^{+}", PARTICLE_LATEX), 200, 0.0, 1.5);
		dHist_FCAL_EoverP_pim = new TH1D(Form("FCAL_EoverP_%s", PARTICLE_MINUS), 
			Form(";E/p %s^{-}", PARTICLE_LATEX), 200, 0.0, 1.5);
		dHist_FCAL_EoverPmeas_pip = new TH1D(Form("FCAL_EoverPmeas_%s", PARTICLE_PLUS), 
			Form(";E_{FCAL}/P_{meas} %s^{+}", PARTICLE_LATEX), 200, 0.0, 1.5);
		dHist_FCAL_EoverPmeas_pim = new TH1D(Form("FCAL_EoverPmeas_%s", PARTICLE_MINUS), 
			Form(";E_{FCAL}/P_{meas} %s^{-}", PARTICLE_LATEX), 200, 0.0, 1.5);
		dHist_Delta_Efcal_kinfitE_vs_kinPmag_plus = new TH2D(
			Form("Delta_Efcal_kinfitE_vs_kinPmag_%s", PARTICLE_PLUS),
			Form(";p_{kinfit} %s^{+};E_{FCAL} - E_{kinfit} %s^{+}", PARTICLE_LATEX, PARTICLE_LATEX),
			200, 0.0, 11.0, 200, -2.0, 2.0);
		dHist_Delta_Efcal_kinfitE_vs_kinPmag_minus = new TH2D(
			Form("Delta_Efcal_kinfitE_vs_kinPmag_%s", PARTICLE_MINUS),
			Form(";p_{kinfit} %s^{-};E_{FCAL} - E_{kinfit} %s^{-}", PARTICLE_LATEX, PARTICLE_LATEX),
			200, 0.0, 11.0, 200, -2.0, 2.0);
	
		dHist_FCAL_Elasticity = new TH1D("FCAL_Elasticity", ";(E1 + E2)/Ebeam", 200, 0.0, 1.4);
		dHist_FCAL_Asymmetry = new TH1D("FCAL_Asymmetry", ";|(E1 - E2)/(E1 + E2)|", 200, 0.0, 1.0);
		//3 FCAL Asymmetry regions: FCAL_Asymmetry < 0.2, 0.2 < FCAL_Asymmetry < 0.5, FCAL_Asymmetry > 0.5. Create these histograms to study each region separately.
		dHist_FCAL_Elasticity_Asym_L0pt2 = new TH1D("FCAL_Elasticity_Asym_L0pt2", ";(E1 + E2)/Ebeam", 200, 0.0, 1.4);
		dHist_FCAL_Elasticity_Asym_G0pt2_L0pt5 = new TH1D("FCAL_Elasticity_Asym_G0pt2_L0pt5", ";(E1 + E2)/Ebeam", 200, 0.0, 1.4);
		dHist_FCAL_Elasticity_Asym_G0pt5 = new TH1D("FCAL_Elasticity_Asym_G0pt5", ";(E1 + E2)/Ebeam", 200, 0.0, 1.4);

		dHist_TrackFCAL_DOCA_plus = new TH1D(Form("FCAL_DOCA_%s", PARTICLE_PLUS), 
			Form(";FCAL_DOCA_%s", PARTICLE_PLUS), 200, 0, 4);
		dHist_TrackFCAL_DOCA_minus = new TH1D(Form("FCAL_DOCA_%s", PARTICLE_MINUS), 
			Form(";FCAL_DOCA_%s", PARTICLE_MINUS), 200, 0, 4);
		
		dHist_FCAL_E1E9_plus = new TH1D(Form("FCAL_E1E9_%s", PARTICLE_PLUS), 
			Form(";(Energy in central cell)/(Energy in 3x3) %s^{+}", PARTICLE_LATEX), 200, 0.0, 1.1);
		dHist_FCAL_E1E9_minus = new TH1D(Form("FCAL_E1E9_%s", PARTICLE_MINUS), 
			Form(";(Energy in central cell)/(Energy in 3x3) %s^{-}", PARTICLE_LATEX), 200, 0.0, 1.1);
		dHist_FCAL_E9E25_plus = new TH1D(Form("FCAL_E9E25_%s", PARTICLE_PLUS), 
			Form(";(Energy in 3x3)/(Energy in 5x5) %s^{+}", PARTICLE_LATEX), 200, 0.0, 1.1);
		dHist_FCAL_E9E25_minus = new TH1D(Form("FCAL_E9E25_%s", PARTICLE_MINUS), 
			Form(";(Energy in 3x3)/(Energy in 5x5) %s^{-}", PARTICLE_LATEX), 200, 0.0, 1.1);

		dHist_FCAL_kin_res_plus = new TH1D(Form("FCAL_kin_res_%s", PARTICLE_PLUS), 
			Form(";|(E_{FCAL} - E_{kinfit})|/E_{kinfit} %s^{+}", PARTICLE_LATEX), 200, 0.0, 1.0);
		dHist_FCAL_kin_res_minus = new TH1D(Form("FCAL_kin_res_%s", PARTICLE_MINUS), 
			Form(";|(E_{FCAL} - E_{kinfit})|/E_{kinfit} %s^{-}", PARTICLE_LATEX), 200, 0.0, 1.0);

		dHist_FCAL_meas_res_plus = new TH1D(Form("FCAL_meas_res_%s", PARTICLE_PLUS), 
			Form(";|(E_{FCAL} - E_{meas})|/E_{meas} %s^{+}", PARTICLE_LATEX), 200, 0.0, 1.0);
		dHist_FCAL_meas_res_minus = new TH1D(Form("FCAL_meas_res_%s", PARTICLE_MINUS), 
			Form(";|(E_{FCAL} - E_{meas})|/E_{meas} %s^{-}", PARTICLE_LATEX), 200, 0.0, 1.0);

		dHist_FCAL_Saturation_vs_Eshower_plus = new TH2D(Form("FCAL_saturation_vs_Eshower_%s", PARTICLE_PLUS), 
			";E shower ;Eshower*E1E9*E9E25", 100, 0, 10, 200, 0, 14);

		dHist_FCAL_Saturation_vs_Eshower_minus = new TH2D(Form("FCAL_saturation_vs_Eshower_%s", PARTICLE_MINUS), 
			";E shower ;Eshower*E1E9*E9E25", 100, 0, 10, 200, 0, 14);

		dHist_FCAL_Saturation_plus = new TH1D(Form("FCAL_saturation_%s", PARTICLE_PLUS), 
			";Eshower*E1E9*E9E25", 200, 0, 14);

		dHist_FCAL_Saturation_minus = new TH1D(Form("FCAL_saturation_%s", PARTICLE_MINUS), 
			";Eshower*E1E9*E9E25", 200, 0, 14);

		dHist_FCAL_SumU_plus = new TH1D(Form("FCAL_SumU_%s", PARTICLE_PLUS), 
			Form(";SumU %s^{+}", PARTICLE_LATEX), 200, 0., 14);

		dHist_FCAL_SumU_minus = new TH1D(Form("FCAL_SumU_%s", PARTICLE_MINUS), 
			Form(";SumU %s^{-}", PARTICLE_LATEX), 200, 0., 14);

		dHist_FCAL_SumV_plus = new TH1D(Form("FCAL_SumV_%s", PARTICLE_PLUS), 
			Form(";SumV %s^{+}", PARTICLE_LATEX), 200, 0., 20);

		dHist_FCAL_SumV_minus = new TH1D(Form("FCAL_SumV_%s", PARTICLE_MINUS), 
			Form(";SumV %s^{-}", PARTICLE_LATEX), 200, 0., 20);

		dHist_FCAL_UV_Asymmetry_plus = new TH1D(Form("FCAL_UV_Asymmetry_%s", PARTICLE_PLUS), 
			Form(";A_{uv} = |#sigma^{2}_{u} - #sigma^{2}_{v}|/|#sigma^{2}_{u} + #sigma^{2}_{v}| %s^{+}", PARTICLE_LATEX), 200, 0., 1.0);

		dHist_FCAL_UV_Asymmetry_minus = new TH1D(Form("FCAL_UV_Asymmetry_%s", PARTICLE_MINUS), 
			Form(";A_{uv} = |#sigma^{2}_{u} - #sigma^{2}_{v}|/|#sigma^{2}_{u} + #sigma^{2}_{v}| %s^{-}", PARTICLE_LATEX), 200, 0., 1.0);

		if(dIsMC) {
			dHist_FCAL_Energy_pip_PostCorr = new TH1D(Form("FCAL_Energy_%s_PostCorr", PARTICLE_PLUS), ";GeV", 200, 0, 10);
			dHist_FCAL_Energy_pim_PostCorr = new TH1D(Form("FCAL_Energy_%s_PostCorr", PARTICLE_MINUS), ";GeV", 200, 0, 10);
			dHist_FCAL_EoverP_pip_PostCorr = new TH1D(Form("FCAL_EoverP_%s_PostCorr", PARTICLE_PLUS), Form(";E/p %s^{+} post-corr", PARTICLE_LATEX), 200, 0.0, 1.5);
			dHist_FCAL_EoverP_pim_PostCorr = new TH1D(Form("FCAL_EoverP_%s_PostCorr", PARTICLE_MINUS), Form(";E/p %s^{-} post-corr", PARTICLE_LATEX), 200, 0.0, 1.5);
			dHist_FCAL_EoverPmeas_pip_PostCorr = new TH1D(Form("FCAL_EoverPmeas_%s_PostCorr", PARTICLE_PLUS), Form(";E_{FCAL}/P_{meas} %s^{+} post-corr", PARTICLE_LATEX), 200, 0.0, 1.5);
			dHist_FCAL_EoverPmeas_pim_PostCorr = new TH1D(Form("FCAL_EoverPmeas_%s_PostCorr", PARTICLE_MINUS), Form(";E_{FCAL}/P_{meas} %s^{-} post-corr", PARTICLE_LATEX), 200, 0.0, 1.5);
			dHist_Delta_Efcal_kinfitE_vs_kinPmag_plus_PostCorr = new TH2D(Form("Delta_Efcal_kinfitE_vs_kinPmag_%s_PostCorr", PARTICLE_PLUS), Form(";p_{kinfit} %s^{+};E_{FCAL}^{post} - E_{kinfit} %s^{+}", PARTICLE_LATEX, PARTICLE_LATEX), 200, 0.0, 11.0, 200, -2.0, 2.0);
			dHist_Delta_Efcal_kinfitE_vs_kinPmag_minus_PostCorr = new TH2D(Form("Delta_Efcal_kinfitE_vs_kinPmag_%s_PostCorr", PARTICLE_MINUS), Form(";p_{kinfit} %s^{-};E_{FCAL}^{post} - E_{kinfit} %s^{-}", PARTICLE_LATEX, PARTICLE_LATEX), 200, 0.0, 11.0, 200, -2.0, 2.0);
			dHist_FCAL_Elasticity_PostCorr = new TH1D("FCAL_Elasticity_PostCorr", ";(E1 + E2)/Ebeam post-corr", 200, 0.0, 1.4);
			dHist_FCAL_Asymmetry_PostCorr = new TH1D("FCAL_Asymmetry_PostCorr", ";|(E1 - E2)/(E1 + E2)| post-corr", 200, 0.0, 1.0);
			dHist_FCAL_Elasticity_Asym_L0pt2_PostCorr = new TH1D("FCAL_Elasticity_Asym_L0pt2_PostCorr", ";(E1 + E2)/Ebeam post-corr", 200, 0.0, 1.4);
			dHist_FCAL_Elasticity_Asym_G0pt2_L0pt5_PostCorr = new TH1D("FCAL_Elasticity_Asym_G0pt2_L0pt5_PostCorr", ";(E1 + E2)/Ebeam post-corr", 200, 0.0, 1.4);
			dHist_FCAL_Elasticity_Asym_G0pt5_PostCorr = new TH1D("FCAL_Elasticity_Asym_G0pt5_PostCorr", ";(E1 + E2)/Ebeam post-corr", 200, 0.0, 1.4);
			dHist_FCAL_kin_res_plus_PostCorr = new TH1D(Form("FCAL_kin_res_%s_PostCorr", PARTICLE_PLUS), Form(";|(E_{FCAL}^{post} - E_{kinfit})|/E_{kinfit} %s^{+}", PARTICLE_LATEX), 200, 0.0, 1.0);
			dHist_FCAL_kin_res_minus_PostCorr = new TH1D(Form("FCAL_kin_res_%s_PostCorr", PARTICLE_MINUS), Form(";|(E_{FCAL}^{post} - E_{kinfit})|/E_{kinfit} %s^{-}", PARTICLE_LATEX), 200, 0.0, 1.0);
			dHist_FCAL_meas_res_plus_PostCorr = new TH1D(Form("FCAL_meas_res_%s_PostCorr", PARTICLE_PLUS), Form(";|(E_{FCAL}^{post} - E_{meas})|/E_{meas} %s^{+}", PARTICLE_LATEX), 200, 0.0, 1.0);
			dHist_FCAL_meas_res_minus_PostCorr = new TH1D(Form("FCAL_meas_res_%s_PostCorr", PARTICLE_MINUS), Form(";|(E_{FCAL}^{post} - E_{meas})|/E_{meas} %s^{-}", PARTICLE_LATEX), 200, 0.0, 1.0);
			dHist_FCAL_Saturation_vs_Eshower_plus_PostCorr = new TH2D(Form("FCAL_saturation_vs_Eshower_%s_PostCorr", PARTICLE_PLUS), ";E shower ;Eshower*E1E9*E9E25 post-corr", 100, 0, 10, 200, 0, 14);
			dHist_FCAL_Saturation_vs_Eshower_minus_PostCorr = new TH2D(Form("FCAL_saturation_vs_Eshower_%s_PostCorr", PARTICLE_MINUS), ";E shower ;Eshower*E1E9*E9E25 post-corr", 100, 0, 10, 200, 0, 14);
			dHist_FCAL_Saturation_plus_PostCorr = new TH1D(Form("FCAL_saturation_%s_PostCorr", PARTICLE_PLUS), ";Eshower*E1E9*E9E25 post-corr", 200, 0, 14);
			dHist_FCAL_Saturation_minus_PostCorr = new TH1D(Form("FCAL_saturation_%s_PostCorr", PARTICLE_MINUS), ";Eshower*E1E9*E9E25 post-corr", 200, 0, 14);
		}


		const int kThetaBins = 28;
		const double kThetaMin = 0.0;
		const double kThetaMax = 8.0;
		const int kPkinBins = 32;
		const double kPkinMin = 0.0;
		const double kPkinMax = 11.0;
		const int kqvec2Bins = 10;
		const double kqvec2Min = 0.0;
		const double kqvec2Max = 0.001;

		TDirectory* fcalThetaDirByView[kNumChargeViews] = {dir_FCALvsTheta_both, dir_FCALvsTheta_plus, dir_FCALvsTheta_minus};
		TDirectory* fcalPkinDirByView[kNumChargeViews] = {dir_FCALvsMomentum_both, dir_FCALvsMomentum_plus, dir_FCALvsMomentum_minus};

		for(int chargeView = 0; chargeView < kNumChargeViews; ++chargeView) {
			std::string chargeLabel;
			if(chargeView == kBoth)
				chargeLabel = Form("%s^{#pm}", PARTICLE_LATEX);
			else if(chargeView == kPlus)
				chargeLabel = Form("%s^{+}", PARTICLE_LATEX);
			else
				chargeLabel = Form("%s^{-}", PARTICLE_LATEX);

			fcalThetaDirByView[chargeView]->cd();
			dFCALvsTheta.Energy[chargeView] = new TH2D("Energy_vs_theta", Form(";#theta %s (deg);E_{FCAL} (GeV)", chargeLabel.c_str()), kThetaBins, kThetaMin, kThetaMax, 200, 0, 10);
			dFCALvsTheta.EoverP[chargeView] = new TH2D("EoverP_vs_theta", Form(";#theta %s (deg);E/p", chargeLabel.c_str()), kThetaBins, kThetaMin, kThetaMax, 200, 0.0, 1.5);
			dFCALvsTheta.EoverPmeas[chargeView] = new TH2D("EoverPmeas_vs_theta", Form(";#theta %s (deg);E_{FCAL}/P_{meas}", chargeLabel.c_str()), kThetaBins, kThetaMin, kThetaMax, 200, 0.0, 1.5);
			dFCALvsTheta.DeltaEfcal_kinfitE[chargeView] = new TH2D("DeltaEfcal_kinfitE_vs_theta", Form(";#theta_{kinfit} %s (deg);E_{FCAL} - E_{kinfit} (GeV)", chargeLabel.c_str()), kThetaBins, kThetaMin, kThetaMax, 200, -2.0, 2.0);
			dFCALvsTheta.TrackDOCA[chargeView] = new TH2D("TrackFCAL_DOCA_vs_theta", Form(";#theta %s (deg);FCAL DOCA", chargeLabel.c_str()), kThetaBins, kThetaMin, kThetaMax, 200, 0, 4);
			dFCALvsTheta.E1E9[chargeView] = new TH2D("E1E9_vs_theta", Form(";#theta %s (deg);E1E9", chargeLabel.c_str()), kThetaBins, kThetaMin, kThetaMax, 200, 0.0, 1.1);
			dFCALvsTheta.E9E25[chargeView] = new TH2D("E9E25_vs_theta", Form(";#theta %s (deg);E9E25", chargeLabel.c_str()), kThetaBins, kThetaMin, kThetaMax, 200, 0.0, 1.1);
			dFCALvsTheta.KinRes[chargeView] = new TH2D("KinRes_vs_theta", Form(";#theta %s (deg);|(E_{FCAL} - E_{kinfit})|/E_{kinfit}", chargeLabel.c_str()), kThetaBins, kThetaMin, kThetaMax, 200, 0.0, 1.0);
			dFCALvsTheta.MeasRes[chargeView] = new TH2D("MeasRes_vs_theta", Form(";#theta %s (deg);|(E_{FCAL} - E_{meas})|/E_{meas}", chargeLabel.c_str()), kThetaBins, kThetaMin, kThetaMax, 200, 0.0, 1.0);
			dFCALvsTheta.Saturation[chargeView] = new TH2D("Saturation_vs_theta", Form(";#theta %s (deg);Eshower*E1E9*E9E25", chargeLabel.c_str()), kThetaBins, kThetaMin, kThetaMax, 200, 0, 14);
			dFCALvsTheta.SumU[chargeView] = new TH2D("SumU_vs_theta", Form(";#theta %s (deg);SumU", chargeLabel.c_str()), kThetaBins, kThetaMin, kThetaMax, 200, 0., 14);
			dFCALvsTheta.SumV[chargeView] = new TH2D("SumV_vs_theta", Form(";#theta %s (deg);SumV", chargeLabel.c_str()), kThetaBins, kThetaMin, kThetaMax, 200, 0., 20);
			dFCALvsTheta.UVAsymmetry[chargeView] = new TH2D("UVAsymmetry_vs_theta", Form(";#theta %s (deg);A_{uv}", chargeLabel.c_str()), kThetaBins, kThetaMin, kThetaMax, 200, 0., 1.0);
			if(dIsMC) {
				dFCALvsTheta.Energy_PostCorr[chargeView] = new TH2D("Energy_PostCorr_vs_theta", Form(";#theta %s (deg);E_{FCAL}^{post} (GeV)", chargeLabel.c_str()), kThetaBins, kThetaMin, kThetaMax, 200, 0, 10);
				dFCALvsTheta.EoverP_PostCorr[chargeView] = new TH2D("EoverP_PostCorr_vs_theta", Form(";#theta %s (deg);E/p post-corr", chargeLabel.c_str()), kThetaBins, kThetaMin, kThetaMax, 200, 0.0, 1.5);
				dFCALvsTheta.EoverPmeas_PostCorr[chargeView] = new TH2D("EoverPmeas_PostCorr_vs_theta", Form(";#theta %s (deg);E_{FCAL}/P_{meas} post-corr", chargeLabel.c_str()), kThetaBins, kThetaMin, kThetaMax, 200, 0.0, 1.5);
				dFCALvsTheta.DeltaEfcal_kinfitE_PostCorr[chargeView] = new TH2D("DeltaEfcal_kinfitE_PostCorr_vs_theta", Form(";#theta_{kinfit} %s (deg);E_{FCAL}^{post} - E_{kinfit} (GeV)", chargeLabel.c_str()), kThetaBins, kThetaMin, kThetaMax, 200, -2.0, 2.0);
				dFCALvsTheta.KinRes_PostCorr[chargeView] = new TH2D("KinRes_PostCorr_vs_theta", Form(";#theta %s (deg);|(E_{FCAL}^{post} - E_{kinfit})|/E_{kinfit}", chargeLabel.c_str()), kThetaBins, kThetaMin, kThetaMax, 200, 0.0, 1.0);
				dFCALvsTheta.MeasRes_PostCorr[chargeView] = new TH2D("MeasRes_PostCorr_vs_theta", Form(";#theta %s (deg);|(E_{FCAL}^{post} - E_{meas})|/E_{meas}", chargeLabel.c_str()), kThetaBins, kThetaMin, kThetaMax, 200, 0.0, 1.0);
				dFCALvsTheta.Saturation_PostCorr[chargeView] = new TH2D("Saturation_PostCorr_vs_theta", Form(";#theta %s (deg);Eshower*E1E9*E9E25 post-corr", chargeLabel.c_str()), kThetaBins, kThetaMin, kThetaMax, 200, 0, 14);
			}

			fcalPkinDirByView[chargeView]->cd();
			dFCALvsPkin.Energy[chargeView] = new TH2D("Energy_vs_kinPmag", Form(";|#vec{p}_{kinfit}| %s (GeV/c);E_{FCAL} (GeV)", chargeLabel.c_str()), kPkinBins, kPkinMin, kPkinMax, 200, 0, 10);
			dFCALvsPkin.EoverP[chargeView] = new TH2D("EoverP_vs_kinPmag", Form(";|#vec{p}_{kinfit}| %s (GeV/c);E/p", chargeLabel.c_str()), kPkinBins, kPkinMin, kPkinMax, 200, 0.0, 1.5);
			dFCALvsPkin.EoverPmeas[chargeView] = new TH2D("EoverPmeas_vs_kinPmag", Form(";|#vec{p}_{kinfit}| %s (GeV/c);E_{FCAL}/P_{meas}", chargeLabel.c_str()), kPkinBins, kPkinMin, kPkinMax, 200, 0.0, 1.5);
			dFCALvsPkin.DeltaEfcal_kinfitE[chargeView] = new TH2D("DeltaEfcal_kinfitE_vs_kinPmag", Form(";|#vec{p}_{kinfit}| %s (GeV/c);E_{FCAL} - E_{kinfit} (GeV)", chargeLabel.c_str()), kPkinBins, kPkinMin, kPkinMax, 200, -2.0, 2.0);
			dFCALvsPkin.TrackDOCA[chargeView] = new TH2D("TrackFCAL_DOCA_vs_kinPmag", Form(";|#vec{p}_{kinfit}| %s (GeV/c);FCAL DOCA", chargeLabel.c_str()), kPkinBins, kPkinMin, kPkinMax, 200, 0, 4);
			dFCALvsPkin.E1E9[chargeView] = new TH2D("E1E9_vs_kinPmag", Form(";|#vec{p}_{kinfit}| %s (GeV/c);E1E9", chargeLabel.c_str()), kPkinBins, kPkinMin, kPkinMax, 200, 0.0, 1.1);
			dFCALvsPkin.E9E25[chargeView] = new TH2D("E9E25_vs_kinPmag", Form(";|#vec{p}_{kinfit}| %s (GeV/c);E9E25", chargeLabel.c_str()), kPkinBins, kPkinMin, kPkinMax, 200, 0.0, 1.1);
			dFCALvsPkin.KinRes[chargeView] = new TH2D("KinRes_vs_kinPmag", Form(";|#vec{p}_{kinfit}| %s (GeV/c);|(E_{FCAL} - E_{kinfit})|/E_{kinfit}", chargeLabel.c_str()), kPkinBins, kPkinMin, kPkinMax, 200, 0.0, 1.0);
			dFCALvsPkin.MeasRes[chargeView] = new TH2D("MeasRes_vs_kinPmag", Form(";|#vec{p}_{kinfit}| %s (GeV/c);|(E_{FCAL} - E_{meas})|/E_{meas}", chargeLabel.c_str()), kPkinBins, kPkinMin, kPkinMax, 200, 0.0, 1.0);
			dFCALvsPkin.Saturation[chargeView] = new TH2D("Saturation_vs_kinPmag", Form(";|#vec{p}_{kinfit}| %s (GeV/c);Eshower*E1E9*E9E25", chargeLabel.c_str()), kPkinBins, kPkinMin, kPkinMax, 200, 0, 14);
			dFCALvsPkin.SumU[chargeView] = new TH2D("SumU_vs_kinPmag", Form(";|#vec{p}_{kinfit}| %s (GeV/c);SumU", chargeLabel.c_str()), kPkinBins, kPkinMin, kPkinMax, 200, 0., 14);
			dFCALvsPkin.SumV[chargeView] = new TH2D("SumV_vs_kinPmag", Form(";|#vec{p}_{kinfit}| %s (GeV/c);SumV", chargeLabel.c_str()), kPkinBins, kPkinMin, kPkinMax, 200, 0., 20);
			dFCALvsPkin.UVAsymmetry[chargeView] = new TH2D("UVAsymmetry_vs_kinPmag", Form(";|#vec{p}_{kinfit}| %s (GeV/c);A_{uv}", chargeLabel.c_str()), kPkinBins, kPkinMin, kPkinMax, 200, 0., 1.0);
			if(dIsMC) {
				dFCALvsPkin.Energy_PostCorr[chargeView] = new TH2D("Energy_PostCorr_vs_kinPmag", Form(";|#vec{p}_{kinfit}| %s (GeV/c);E_{FCAL}^{post} (GeV)", chargeLabel.c_str()), kPkinBins, kPkinMin, kPkinMax, 200, 0, 10);
				dFCALvsPkin.EoverP_PostCorr[chargeView] = new TH2D("EoverP_PostCorr_vs_kinPmag", Form(";|#vec{p}_{kinfit}| %s (GeV/c);E/p post-corr", chargeLabel.c_str()), kPkinBins, kPkinMin, kPkinMax, 200, 0.0, 1.5);
				dFCALvsPkin.EoverPmeas_PostCorr[chargeView] = new TH2D("EoverPmeas_PostCorr_vs_kinPmag", Form(";|#vec{p}_{kinfit}| %s (GeV/c);E_{FCAL}/P_{meas} post-corr", chargeLabel.c_str()), kPkinBins, kPkinMin, kPkinMax, 200, 0.0, 1.5);
				dFCALvsPkin.DeltaEfcal_kinfitE_PostCorr[chargeView] = new TH2D("DeltaEfcal_kinfitE_PostCorr_vs_kinPmag", Form(";|#vec{p}_{kinfit}| %s (GeV/c);E_{FCAL}^{post} - E_{kinfit} (GeV)", chargeLabel.c_str()), kPkinBins, kPkinMin, kPkinMax, 200, -2.0, 2.0);
				dFCALvsPkin.KinRes_PostCorr[chargeView] = new TH2D("KinRes_PostCorr_vs_kinPmag", Form(";|#vec{p}_{kinfit}| %s (GeV/c);|(E_{FCAL}^{post} - E_{kinfit})|/E_{kinfit}", chargeLabel.c_str()), kPkinBins, kPkinMin, kPkinMax, 200, 0.0, 1.0);
				dFCALvsPkin.MeasRes_PostCorr[chargeView] = new TH2D("MeasRes_PostCorr_vs_kinPmag", Form(";|#vec{p}_{kinfit}| %s (GeV/c);|(E_{FCAL}^{post} - E_{meas})|/E_{meas}", chargeLabel.c_str()), kPkinBins, kPkinMin, kPkinMax, 200, 0.0, 1.0);
				dFCALvsPkin.Saturation_PostCorr[chargeView] = new TH2D("Saturation_PostCorr_vs_kinPmag", Form(";|#vec{p}_{kinfit}| %s (GeV/c);Eshower*E1E9*E9E25 post-corr", chargeLabel.c_str()), kPkinBins, kPkinMin, kPkinMax, 200, 0, 14);
			}
		}

		TDirectory* dir_FCALvsqvec2_both = dir_FCALvsqvec2->mkdir("FCALvsqvec2_Both");
		TDirectory* dir_FCALvsqvec2_plus = dir_FCALvsqvec2->mkdir(Form("FCALvsqvec2_%s", PARTICLE_PLUS));
		TDirectory* dir_FCALvsqvec2_minus = dir_FCALvsqvec2->mkdir(Form("FCALvsqvec2_%s", PARTICLE_MINUS));
		TDirectory* fcalQvec2DirByView[kNumChargeViews] = {dir_FCALvsqvec2_both, dir_FCALvsqvec2_plus, dir_FCALvsqvec2_minus};

		for(int chargeView = 0; chargeView < kNumChargeViews; ++chargeView) {
			std::string chargeLabel;
			std::string chargeSuffix;
			if(chargeView == kBoth) {
				chargeLabel = Form("%s^{#pm}", PARTICLE_LATEX);
				chargeSuffix = "avg";
			} else if(chargeView == kPlus) {
				chargeLabel = Form("%s^{+}", PARTICLE_LATEX);
				chargeSuffix = Form("%s^{+}", PARTICLE_LATEX);
			} else {
				chargeLabel = Form("%s^{-}", PARTICLE_LATEX);
				chargeSuffix = Form("%s^{-}", PARTICLE_LATEX);
			}

			fcalQvec2DirByView[chargeView]->cd();
			dFCALvsqvec2.Energy[chargeView] = new TH2D("Energy_vs_qvec2", Form(";#vec{q}^{2} (GeV/c)^{2};E_{FCAL}^{%s} (GeV)", chargeSuffix.c_str()), 200, 0.0, 0.04, 200, 0, 10);
			dFCALvsqvec2.EoverP[chargeView] = new TH2D("EoverP_vs_qvec2", Form(";#vec{q}^{2} (GeV/c)^{2};E/p %s", chargeLabel.c_str()), 200, 0.0, 0.04, 200, 0.0, 1.5);
			dFCALvsqvec2.EoverPmeas[chargeView] = new TH2D("EoverPmeas_vs_qvec2", Form(";#vec{q}^{2} (GeV/c)^{2};E_{FCAL}/P_{meas} %s", chargeLabel.c_str()), 200, 0.0, 0.04, 200, 0.0, 1.5);
			dFCALvsqvec2.DeltaEfcal_kinfitE[chargeView] = new TH2D("DeltaEfcal_kinfitE_vs_qvec2", Form(";#vec{q}^{2} (GeV/c)^{2};E_{FCAL}^{%s} - E_{kinfit}^{%s} (GeV)", chargeSuffix.c_str(), chargeSuffix.c_str()), 200, 0.0, 0.04, 200, -2.0, 2.0);
			dFCALvsqvec2.TrackDOCA[chargeView] = new TH2D("TrackFCAL_DOCA_vs_qvec2", Form(";#vec{q}^{2} (GeV/c)^{2};FCAL DOCA %s", chargeLabel.c_str()), 200, 0.0, 0.04, 200, 0, 4);
			dFCALvsqvec2.E1E9[chargeView] = new TH2D("E1E9_vs_qvec2", Form(";#vec{q}^{2} (GeV/c)^{2};E1E9 %s", chargeLabel.c_str()), 200, 0.0, 0.04, 200, 0.0, 1.1);
			dFCALvsqvec2.E9E25[chargeView] = new TH2D("E9E25_vs_qvec2", Form(";#vec{q}^{2} (GeV/c)^{2};E9E25 %s", chargeLabel.c_str()), 200, 0.0, 0.04, 200, 0.0, 1.1);
			dFCALvsqvec2.KinRes[chargeView] = new TH2D("KinRes_vs_qvec2", Form(";#vec{q}^{2} (GeV/c)^{2};|(E_{FCAL} - E_{kinfit})|/E_{kinfit} %s", chargeLabel.c_str()), 200, 0.0, 0.04, 200, 0.0, 1.0);
			dFCALvsqvec2.MeasRes[chargeView] = new TH2D("MeasRes_vs_qvec2", Form(";#vec{q}^{2} (GeV/c)^{2};|(E_{FCAL} - E_{meas})|/E_{meas} %s", chargeLabel.c_str()), 200, 0.0, 0.04, 200, 0.0, 1.0);
			dFCALvsqvec2.Saturation[chargeView] = new TH2D("Saturation_vs_qvec2", Form(";#vec{q}^{2} (GeV/c)^{2};Eshower*E1E9*E9E25 %s", chargeLabel.c_str()), 200, 0.0, 0.04, 200, 0, 14);
			dFCALvsqvec2.SumU[chargeView] = new TH2D("SumU_vs_qvec2", Form(";#vec{q}^{2} (GeV/c)^{2};SumU %s", chargeLabel.c_str()), 200, 0.0, 0.04, 200, 0., 14);
			dFCALvsqvec2.SumV[chargeView] = new TH2D("SumV_vs_qvec2", Form(";#vec{q}^{2} (GeV/c)^{2};SumV %s", chargeLabel.c_str()), 200, 0.0, 0.04, 200, 0., 20);
			dFCALvsqvec2.UVAsymmetry[chargeView] = new TH2D("UVAsymmetry_vs_qvec2", Form(";#vec{q}^{2} (GeV/c)^{2};A_{uv} %s", chargeLabel.c_str()), 200, 0.0, 0.04, 200, 0., 1.0);

			if(dIsMC) {
				dFCALvsqvec2.Energy_PostCorr[chargeView] = new TH2D("Energy_PostCorr_vs_qvec2", Form(";#vec{q}^{2} (GeV/c)^{2};E_{FCAL}^{%s} post-corr (GeV)", chargeSuffix.c_str()), 200, 0.0, 0.04, 200, 0, 10);
				dFCALvsqvec2.EoverP_PostCorr[chargeView] = new TH2D("EoverP_PostCorr_vs_qvec2", Form(";#vec{q}^{2} (GeV/c)^{2};E/p %s post-corr", chargeLabel.c_str()), 200, 0.0, 0.04, 200, 0.0, 1.5);
				dFCALvsqvec2.EoverPmeas_PostCorr[chargeView] = new TH2D("EoverPmeas_PostCorr_vs_qvec2", Form(";#vec{q}^{2} (GeV/c)^{2};E_{FCAL}/P_{meas} %s post-corr", chargeLabel.c_str()), 200, 0.0, 0.04, 200, 0.0, 1.5);
				dFCALvsqvec2.DeltaEfcal_kinfitE_PostCorr[chargeView] = new TH2D("DeltaEfcal_kinfitE_PostCorr_vs_qvec2", Form(";#vec{q}^{2} (GeV/c)^{2};E_{FCAL}^{%s} post-corr - E_{kinfit}^{%s} (GeV)", chargeSuffix.c_str(), chargeSuffix.c_str()), 200, 0.0, 0.04, 200, -2.0, 2.0);
				dFCALvsqvec2.KinRes_PostCorr[chargeView] = new TH2D("KinRes_PostCorr_vs_qvec2", Form(";#vec{q}^{2} (GeV/c)^{2};|(E_{FCAL} - E_{kinfit})|/E_{kinfit} %s post-corr", chargeLabel.c_str()), 200, 0.0, 0.04, 200, 0.0, 1.0);
				dFCALvsqvec2.MeasRes_PostCorr[chargeView] = new TH2D("MeasRes_PostCorr_vs_qvec2", Form(";#vec{q}^{2} (GeV/c)^{2};|(E_{FCAL} - E_{meas})|/E_{meas} %s post-corr", chargeLabel.c_str()), 200, 0.0, 0.04, 200, 0.0, 1.0);
				dFCALvsqvec2.Saturation_PostCorr[chargeView] = new TH2D("Saturation_PostCorr_vs_qvec2", Form(";#vec{q}^{2} (GeV/c)^{2};Eshower*E1E9*E9E25 %s post-corr", chargeLabel.c_str()), 200, 0.0, 0.04, 200, 0, 14);
			}
		}

	



	
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//      ~~~~         ~~~~
	//    ~~~   ~~~   ~~~    ~~~
	//  ~~~       ~~~~         ~~~
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		

	// Return to the analysis top-level directory so BookSystematics is created in
	// LeptonPairStudies/systematics (not under FCAL/systematics).
	dir_Main->cd();
	BookSystematics();
	

	/************************************** USER INITIALIZATION: CUT THRESHOLDS *************************************/

	// Runtime histogram/fill switches for staged debugging.
	// Toggle these booleans to enable directories/features incrementally.
	dRuntimeFillSwitches.fillPreFidDiagnostics = true;
	dRuntimeFillSwitches.runPostPreFidBlocks = true;
	dRuntimeFillSwitches.fillSystematics = true;
	dRuntimeFillSwitches.fillThrownResolution = true;
	dRuntimeFillSwitches.fillCategoryDirectories = true;
	dRuntimeFillSwitches.fillFCALStudy = true;
	dRuntimeFillSwitches.fillNminus1 = true;
	
	ApplyFCALEnergyCorrections = false;       // Correct simulated FCAL energies based on track momentum
	// TEMPORARY: Turned off to speed up compilation (no DEPIClassifiers.h)
	ApplyMLPClassification = true;           // Apply MLP classifier to get particle ID scores (turn OFF when producing training files)
	// fMLPModelType = EPI_MLP;               // EPI_MLP = comprehensive 8-input (default), EPI_MLP_OLD = legacy 3-input
	ApplyBDTClassification = true;           // Apply BDT classifier to get particle ID scores (turn OFF when producing training files)

	// Consolidated process-level analysis cut settings.
	// This is the single source of truth for fiducial/analysis thresholds.
	dAnalysisCutConditions = {};

	// Section 0: Preselection cut emulation
	dAnalysisCutConditions.applyPreselectionEoPCut = false; // 0a: E/p > 0.4 preselection (off by default)
	dAnalysisCutConditions.preselectionMinEoverP = 0.4;

	// Section 1: Exclusivity cuts
	dAnalysisCutConditions.applyNoUnusedTracksCut = true; // 1a: No unused charged tracks
	dAnalysisCutConditions.maxNumUnusedTracks = 0;
	dAnalysisCutConditions.applyNoUnusedShowersCut = true; // 1b: No unused neutral shower energy
	dAnalysisCutConditions.maxUnusedShowersEnergy = 0.0;

	// Section 2: Forward detector cuts
	dAnalysisCutConditions.applyFCALEnergyNonZeroCut = true; // 2a: FCAL E > 0 for both tracks
	dAnalysisCutConditions.applyTOFdEdxNonZeroCut = true; // 2b: TOF dE/dx > 0 for both tracks

	// Section 3: Fiducial cuts
	dAnalysisCutConditions.applyMinBeamEcut = true; // 3a: Coherent peak is 8.2-8.8 GeV. FFStudy: 7-11.4 GeV
	dAnalysisCutConditions.applyMaxBeamEcut = true; // in practice, we have two beam regions automatically overwritten with CoherentPeak and FullSpectrum, so for most histograms, these cuts are overwritten. They are set to the FFStudy settings now
	// but CoherentPeak will still be filled automatically.
	dAnalysisCutConditions.minBeamE = 7.0;
	dAnalysisCutConditions.maxBeamE = 11.4;
	dAnalysisCutConditions.applyMinThetaCut = true; // 3b: Minimum theta > 1.5 deg
	dAnalysisCutConditions.applyMaxThetaCut = false; // 3c: Maximum theta < 7.5 deg
	dAnalysisCutConditions.minTheta = 1.5;
	dAnalysisCutConditions.maxTheta = 7.5;
	dAnalysisCutConditions.apply2DThetaCut = false; // 3d: 2D theta acceptance match
	dAnalysisCutConditions.applyMomentumRangeCut = true; // 3e: Valid FCAL correction range (0.45-7.92 GeV)
	dAnalysisCutConditions.minPmag = 0.45;
	dAnalysisCutConditions.maxPmag = 7.92;
	dAnalysisCutConditions.applyMinWkinCut = true; // 3f: Invariant mass window
	dAnalysisCutConditions.applyMaxWkinCut = true;
	dAnalysisCutConditions.minWkin = 0.700;
	dAnalysisCutConditions.maxWkin = 0.770;
	dAnalysisCutConditions.applyVertexZCut = true; // 3g: Vertex Z position (52-78 cm)
	dAnalysisCutConditions.minVertexZ = 52.0;
	dAnalysisCutConditions.maxVertexZ = 78.0;
	dAnalysisCutConditions.applyKinFitCLCut = true; // 3h: Kinematic fit CL > 1e-6
	dAnalysisCutConditions.minKinFitCL = 1e-6;
	dAnalysisCutConditions.applyBestChiSqComboCut = true; // 3i: in-time photon with best chi^2 from kinematic fit only
	dAnalysisCutConditions.applyMVACuts = true;
	dAnalysisCutConditions.modelChoice = "MLP";
	dAnalysisCutConditions.particleChoice = "ee";
	dAnalysisCutConditions.BDT_ee = 0.053080951;
	dAnalysisCutConditions.BDT_pi = -0.069788195;
	dAnalysisCutConditions.MLP_ee = 0.8;
	dAnalysisCutConditions.MLP_pi = 0.4;

	// Additional cuts for testing--off in any official analyses
	dAnalysisCutConditions.applyMinEoverP1Cut = false;
	dAnalysisCutConditions.minEoverP1 = 0.7;
	dAnalysisCutConditions.applyMaxEoverP1Cut = false;
	dAnalysisCutConditions.maxEoverP1 = 1.2;
	dAnalysisCutConditions.applyMinEoverP2Cut = false;
	dAnalysisCutConditions.minEoverP2 = 0.7;
	dAnalysisCutConditions.applyMaxEoverP2Cut = false;
	dAnalysisCutConditions.maxEoverP2 = 1.2;
	dAnalysisCutConditions.applyMinFCALElasticityCut = false;
	dAnalysisCutConditions.minFCALElasticity = 0.8;
	dAnalysisCutConditions.applyMaxFCALElasticityCut = false;
	dAnalysisCutConditions.maxFCALElasticity = 1.2;
	dAnalysisCutConditions.applyMaxFCALDOCACut = false;
	dAnalysisCutConditions.maxFCALDOCA = 2.0;
	dAnalysisCutConditions.applyMaxTOFdEdxCut = false;
	dAnalysisCutConditions.maxTOFdEdx = 2.0 * 0.00263;
	dAnalysisCutConditions.applyMLPResponseCut = false;
	dAnalysisCutConditions.minMLPResponse = 0.8;

	// CPP-specific analysis/fiducial cut settings.
	// Keep this directly below the GlueX block so both experiment definitions
	// are edited in one obvious place.
	dAnalysisCutConditionsCPP = dAnalysisCutConditions;
	dAnalysisCutConditionsCPP.minBeamE = 5.35;
	dAnalysisCutConditionsCPP.maxBeamE = 5.75;
	dAnalysisCutConditionsCPP.minTheta = 1.1;
	dAnalysisCutConditionsCPP.apply2DThetaCut = false;
	dAnalysisCutConditionsCPP.minWkin = 0.0;
	dAnalysisCutConditionsCPP.minVertexZ = -10.0;
	dAnalysisCutConditionsCPP.maxVertexZ = 10.0;
	// ======================================================================
	// ======================================================================
	
	//	======================================================================
	// FCAL non-linear energy corrections

	ep_FCAL_cor1 = new TF1("ep_FCAL_0pt45to3pt6", 
		"-0.0815097 + 0.113626*x -0.112878*x*x + 0.0455309*x*x*x -0.00615324*x*x*x*x",
		0.0,4.0);
	ep_FCAL_cor2 = new TF1("ep_FCAL_3pt6to8pt01",
		"140.568 - 189.004*x + 106.853*x*x -32.8906*x*x*x + 5.94911*x*x*x*x -0.63234*x*x*x*x*x +0.0365855*x*x*x*x*x*x -0.000889301*x*x*x*x*x*x*x",
		0.0,10.0);

	ep_FCAL_cor2data1 = new TF1("ep_FCAL_Correct0pt36to3pt6", 
		"-0.102153 + 0.167009*x -0.168846*x*x + 0.059494*x*x*x -0.00709742*x*x*x*x",
		0.0, 4.0);
	ep_FCAL_cor2data2 = new TF1("ep_FCAL_Correct3pt6andup", 
		"-15.5234 + 14.3966*x - 5.18514*x*x + 0.904002*x*x*x -0.0770312*x*x*x*x + 0.00257271*x*x*x*x*x"
		,0.0,10.0);

	em_FCAL_cor1 = new TF1("em_FCAL_0pt45to3pt6", 
		"-0.052032 + 0.00782447*x",
		0.0, 4.0);
	em_FCAL_cor2 = new TF1("em_FCAL_3pt6to7pt92",
		"18.2227 -22.007*x + 11.0161*x*x -2.90704*x*x*x + 0.423587*x*x*x*x -0.0322319*x*x*x*x*x + 0.000999988*x*x*x*x*x*x",
		0.0, 10.0);

	em_FCAL_cor2data1 = new TF1("em_FCAL_Correct0pt36to3pt6",
		"-0.0788966 + 0.191044*x - 0.283841*x*x + 0.153231*x*x*x -0.0371609*x*x*x*x + 0.00331943*x*x*x*x*x",
		0.0, 4.0);
	em_FCAL_cor2data2 = new TF1("em_FCAL_Correct3pt6andup",
		"-2.53871 + 2.08889*x -0.630479*x*x + 0.0775441*x*x*x -0.00343892*x*x*x*x", 
		0.0, 10.0);

	/******************************** TMVA READER INITIALIZATION *******************************/
	
	std::cout << "Initializing TMVA Readers for e+/e- PID..." << std::endl;
	
	// Initialize TMVA Reader for MLP classifier
	dTMVAReader_MLP = new TMVA::Reader("!Color:!Silent");
	std::cout << "TMVA MLP Reader created." << std::endl;
	
	// Add input variables to the reader (must match training variable names and order)
	dTMVAReader_MLP->AddVariable("EoverPkin_FCAL", &dTMVA_EoverPkin_FCAL);
	dTMVAReader_MLP->AddVariable("TrackFCAL_DOCA", &dTMVA_TrackFCAL_DOCA);
	dTMVAReader_MLP->AddVariable("E9E25_FCAL", &dTMVA_E9E25_FCAL);
	dTMVAReader_MLP->AddVariable("E1E9_FCAL", &dTMVA_E1E9_FCAL);
	dTMVAReader_MLP->AddVariable("SumU_FCAL", &dTMVA_SumU_FCAL);
	dTMVAReader_MLP->AddVariable("SumV_FCAL", &dTMVA_SumV_FCAL);
	dTMVAReader_MLP->AddVariable("UV_Asymmetry_FCAL", &dTMVA_UV_Asymmetry_FCAL);
	dTMVAReader_MLP->AddVariable("Saturation_FCAL", &dTMVA_Saturation_FCAL);
	std::cout << "TMVA MLP variables added." << std::endl;
	
	// Book the MVA method from XML weights file
	// Method name must match what's in the XML: "MLP::MLP"
	std::cout << "Booking MLP method from TMVAClassification_MLP.weights.xml..." << std::endl;
	dTMVAReader_MLP->BookMVA("MLP::MLP", "/work/halld/home/acschick/tmva-hpo/jobs/job_10005/dataset/weights/TMVAClassification_MLP.weights.xml");
	std::cout << "TMVA MLP Reader successfully initialized!" << std::endl;
	
	// Initialize TMVA Reader for BDT classifier
	dTMVAReader_BDT = new TMVA::Reader("!Color:!Silent");
	std::cout << "TMVA BDT Reader created." << std::endl;
	
	// Add input variables to the BDT reader (must use SAME variables as MLP)
	dTMVAReader_BDT->AddVariable("EoverPkin_FCAL", &dTMVA_EoverPkin_FCAL);
	dTMVAReader_BDT->AddVariable("TrackFCAL_DOCA", &dTMVA_TrackFCAL_DOCA);
	dTMVAReader_BDT->AddVariable("E9E25_FCAL", &dTMVA_E9E25_FCAL);
	dTMVAReader_BDT->AddVariable("E1E9_FCAL", &dTMVA_E1E9_FCAL);
	dTMVAReader_BDT->AddVariable("SumU_FCAL", &dTMVA_SumU_FCAL);
	dTMVAReader_BDT->AddVariable("SumV_FCAL", &dTMVA_SumV_FCAL);
	dTMVAReader_BDT->AddVariable("UV_Asymmetry_FCAL", &dTMVA_UV_Asymmetry_FCAL);
	dTMVAReader_BDT->AddVariable("Saturation_FCAL", &dTMVA_Saturation_FCAL);
	std::cout << "TMVA BDT variables added." << std::endl;
	
	// Book the BDT MVA method from XML weights file
	std::cout << "Booking BDT method from TMVAClassification_BDT.weights.xml..." << std::endl;
	dTMVAReader_BDT->BookMVA("BDT", "/work/halld/home/acschick/tmva-hpo/jobs/job_10001/dataset/weights/TMVAClassification_BDT.weights.xml");
	std::cout << "TMVA BDT Reader successfully initialized!" << std::endl;


	/************************************* ADVANCED EXAMPLE: CHOOSE BRANCHES TO READ ************************************/

	//TO SAVE PROCESSING TIME
		//If you know you don't need all of the branches/data, but just a subset of it, you can speed things up
		//By default, for each event, the data is retrieved for all branches
		//If you know you only need data for some branches, you can skip grabbing data from the branches you don't need
		//Do this by doing something similar to the commented code below

	//dTreeInterface->Clear_GetEntryBranches(); //now get none
	//dTreeInterface->Register_GetEntryBranch("Proton__P4"); //manually set the branches you want
}

DSelector_2eMissingProton_Systematics::CutConditions DSelector_2eMissingProton_Systematics::BuildActiveFiducialConditions(void) const
{
	return dIsCPPRunPeriod ? dAnalysisCutConditionsCPP : dAnalysisCutConditions;
}

Bool_t DSelector_2eMissingProton_Systematics::Process(Long64_t locEntry)
{
	// The Process() function is called for each entry in the tree. The entry argument
	// specifies which entry in the currently loaded tree is to be processed.
	//
	// This function should contain the "body" of the analysis. It can contain
	// simple or elaborate selection criteria, run algorithms on the data
	// of the event and typically fill histograms.
	//
	// The processing can be stopped by calling Abort().
	// Use fStatus to set the return value of TTree::Process().
	// The return value is currently not used.

	//CALL THIS FIRST
	DSelector::Process(locEntry); //Gets the data from the tree for the entry

	const bool kRunPostPreFidBlocks = dRuntimeFillSwitches.runPostPreFidBlocks;
	//const bool kMinimalStableMode = true I think is now the same as just turning off all the fill switches (except maybe fillCategoryDirectories, which is needed for some of the systematics)

	// DEBUG: Print progress every 1000 events
	if(locEntry % 1000 == 0) {
		// cout << "Processing entry: " << locEntry << endl;
	}


	/******************************************** GET POLARIZATION ORIENTATION ******************************************/

	//Only if the run number changes
	//RCDB environment must be setup in order for this to work! (Will return false otherwise)
	UInt_t locRunNumber = Get_RunNumber();
	if(locRunNumber != dPreviousRunNumber)
	{
		dIsCPPRunPeriod = (locRunNumber >= 100531 && locRunNumber <= 101622);
		if(locRunNumber >= 40856 && locRunNumber <= 42550)
			dRunPeriodIndex = 0; // 1801 official range
		else if(locRunNumber >= 50685 && locRunNumber <= 51768)
			dRunPeriodIndex = 1; // 1808 range
		else if(dIsCPPRunPeriod)
			dRunPeriodIndex = 2; // CPP
		else
			dRunPeriodIndex = -1;
		dIsPolarizedFlag = dAnalysisUtilities.Get_IsPolarizedBeam(locRunNumber, dIsPARAFlag);
		dAnalysisUtilities.Get_PolarizationAngle(locRunNumber, dPolarizationAngle);
		dPreviousRunNumber = locRunNumber;
		//cout << "Run " << locRunNumber << ": IsPolarized = " << dIsPolarizedFlag << ", IsPARA = " << dIsPARAFlag
		//     << ", PolarizationAngle = " << dPolarizationAngle << " deg" << endl;
	}

	const int activeRunPeriodIndex = dRunPeriodIndex;
	const CutConditions activeFiducialConditions = BuildActiveFiducialConditions();

	// Map run polarization angle to histogram slot: 0,45,90,135,AMO.
	const int activePolarizationIndex = [&]() -> int {
		if(!dIsPolarizedFlag)
			return 4; // AMO
		const int normalizedAngle = ((dPolarizationAngle % 180) + 180) % 180;
		switch(normalizedAngle) {
			case 0: return 0;
			case 45: return 1;
			case 90: return 2;
			case 135: return 3;
			default: return 4; // unknown polarized angle: treat as AMO bin
		}
	}();

	/********************************************* SETUP UNIQUENESS TRACKING ********************************************/

	//ANALYSIS ACTIONS: Reset uniqueness tracking for each action
	//For any actions that you are executing manually, be sure to call Reset_NewEvent() on them here
	Reset_Actions_NewEvent();
	//dAnalyzeCutActions->Reset_NewEvent(); // manual action, must call Reset_NewEvent()

	//PREVENT-DOUBLE COUNTING WHEN HISTOGRAMMING
	//Sometimes, some content is the exact same between one combo and the next
	//e.g. maybe two combos have different beam particles, but the same data for the final-state
	//When histogramming, you don't want to double-count when this happens: artificially inflates your signal (or background)
	//So, for each quantity you histogram, keep track of what particles you used (for a given combo)
	//Then for each combo, just compare to what you used before, and make sure it's unique

	//EXAMPLE 1: Particle-specific info:
	set<Int_t> locUsedSoFar_BeamEnergy; //Int_t: Unique ID for beam particles. set: easy to use, fast to search
	set<Int_t> locUsedSoFar_Positron, locUsedSoFar_Electron;

	//EXAMPLE 2: Combo-specific info:
	//In general: Could have multiple particles with the same PID: Use a set of Int_t's
	//In general: Multiple PIDs, so multiple sets: Contain within a map
	//Multiple combos: Contain maps within a set (easier, faster to search)
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_2e, locUsedSoFar_Angles;

	//INSERT USER ANALYSIS UNIQUENESS TRACKING HERE

	//EXAMPLE 0: Event-specific info:


	/******************************************* LOOP OVER THROWN DATA (OPTIONAL) ***************************************/
	
	// Declare thrown MC variables outside the dThrownBeam block so they're in scope for entire Process() function
	Bool_t foundThrownParticles = false;
	Double_t qvec2_thrown = -999.0;
	Double_t theta_plus_thrown = -999.0;
	Double_t theta_minus_thrown = -999.0;
	Double_t Pmag_plus_thrown = -999.0;
	Double_t Pmag_minus_thrown = -999.0;
	Double_t Wepem_thrown = -999.0;
	Double_t phi_plus_thrown = -999.0;
	Double_t phi_minus_thrown = -999.0;
	Double_t JTphi_thrown = -999.0;
	Double_t theta_recoil_thrown = -999.0;
	Double_t pmag_recoil_thrown = -999.0;
	Bool_t foundThrownProton = false;
	                                                                
        double locEbeam_Thrown = 0;
        if(dThrownBeam != NULL) {
	  		locEbeam_Thrown = dThrownBeam->Get_P4().T();
        	TLorentzVector locBeamP4_Thrown = dThrownBeam->Get_P4();
        	
        	// Access thrown particles directly by index
        	// For reaction gamma p -> e+ e- (p missing), thrown final state particles are at fixed indices
        	TLorentzVector locProton_Thrown;
        	TLorentzVector locPositronP4_Thrown;
        	TLorentzVector locElectronP4_Thrown;
		
			// maybe need to revisit this.
		if (Get_NumThrown() >= 2) {
			for (UInt_t locThrownIndex = 0; locThrownIndex < Get_NumThrown(); ++locThrownIndex) {
				dThrownWrapper->Set_ArrayIndex(locThrownIndex);
				if (dThrownWrapper->Get_PID() == Proton) {
					locProton_Thrown = dThrownWrapper->Get_P4();
					foundThrownProton = true;
					break;
				}
			}

			dThrownWrapper->Set_ArrayIndex(0);
			TLorentzVector locP4_Thrown_0 = dThrownWrapper->Get_P4();
			Particle_t pid_0 = dThrownWrapper->Get_PID();
			
			dThrownWrapper->Set_ArrayIndex(1);
			TLorentzVector locP4_Thrown_1 = dThrownWrapper->Get_P4();
			Particle_t pid_1 = dThrownWrapper->Get_PID();
			
			// Assign based on PID
			if (pid_0 == Positron && pid_1 == Electron) {
				locPositronP4_Thrown = locP4_Thrown_0;
				locElectronP4_Thrown = locP4_Thrown_1;
				foundThrownParticles = true;
			} else if (pid_0 == Electron && pid_1 == Positron) {
				locElectronP4_Thrown = locP4_Thrown_0;
				locPositronP4_Thrown = locP4_Thrown_1;
				foundThrownParticles = true;
			}
		}

		if (foundThrownParticles) {
			qvec2_thrown = (locBeamP4_Thrown.Vect() - locPositronP4_Thrown.Vect() - locElectronP4_Thrown.Vect()).Mag2();
			theta_plus_thrown = locPositronP4_Thrown.Theta() * TMath::RadToDeg();
			theta_minus_thrown = locElectronP4_Thrown.Theta() * TMath::RadToDeg();
			Pmag_plus_thrown = locPositronP4_Thrown.Vect().Mag();
			Pmag_minus_thrown = locElectronP4_Thrown.Vect().Mag();
			
			// Invariant mass
			TLorentzVector locDileptonP4_Thrown = locPositronP4_Thrown + locElectronP4_Thrown;
			Wepem_thrown = locDileptonP4_Thrown.M();
			
			// JTphi calculation for thrown (using same formula as kinfit variables)
			Double_t p1_theta_thrown = locPositronP4_Thrown.Theta();
			Double_t p2_theta_thrown = locElectronP4_Thrown.Theta();
			Double_t JTx_thrown = locPositronP4_Thrown.X() * 2*locElectronP4_Thrown.E()/(locPositronP4_Thrown.E() - Pmag_plus_thrown * cos(p1_theta_thrown))  +
								  locElectronP4_Thrown.X() * 2*locPositronP4_Thrown.E()/(locElectronP4_Thrown.E() - Pmag_minus_thrown * cos(p2_theta_thrown));
			Double_t JTy_thrown = locPositronP4_Thrown.Y() * 2*locElectronP4_Thrown.E()/(locPositronP4_Thrown.E() - Pmag_plus_thrown * cos(p1_theta_thrown))  +
								  locElectronP4_Thrown.Y() * 2*locPositronP4_Thrown.E()/(locElectronP4_Thrown.E() - Pmag_minus_thrown * cos(p2_theta_thrown));
			JTphi_thrown = atan2(JTy_thrown, JTx_thrown) * TMath::RadToDeg();
			
			phi_plus_thrown = locPositronP4_Thrown.Phi() * TMath::RadToDeg();
			phi_minus_thrown = locElectronP4_Thrown.Phi() * TMath::RadToDeg();

			if (foundThrownProton) {
				theta_recoil_thrown = locProton_Thrown.Theta() * TMath::RadToDeg();
				pmag_recoil_thrown = locProton_Thrown.Vect().Mag();
			}
		}
	}

	/************************************************* LOOP OVER COMBOS *************************************************/

	// Vector to store combo information
	std::vector<std::pair<UInt_t, Double_t>> loc_combos;

	// Pre-loop to gather kinfit ComboIndex-chiSq pairing and sort by chiSq value ascendingly
	for(UInt_t loc_i = 0; loc_i < Get_NumCombos(); ++loc_i)
	{
		dComboWrapper->Set_ComboIndex(loc_i);
		Double_t locChiSq = dComboWrapper->Get_ChiSq_KinFit("");
		loc_combos.push_back(std::make_pair(loc_i, locChiSq));
	}
	
	// DEBUG: Print combo count for first few events
	if(locEntry < 10 || locEntry % 1000 == 0) {
		// cout << "  Entry " << locEntry << " has " << Get_NumCombos() << " combos" << endl;
	}
	// Sort the combos by ChiSq
	std::sort(loc_combos.begin(), loc_combos.end(), [](const std::pair<UInt_t, Double_t>& a, const std::pair<UInt_t, Double_t>& b) {
		return a.second < b.second;
	});

	//Loop over combos
	// Track whether we've found the best in-time combo yet (for BestChiSq fills)
	bool foundBestInTimeCombo = false;

	for(const auto& loc_combo : loc_combos){
		UInt_t loc_i = loc_combo.first; // .first 
		//Set branch array indices for combo and all combo particles
		dComboWrapper->Set_ComboIndex(loc_i);
		ThisComboIsBestChiSq = false;
		
		// Is used to indicate when combos have been cut
		if(dComboWrapper->Get_IsComboCut()) // Is false when tree originally created
		continue; // Combo has been cut previously

		/********************************************** GET PARTICLE INDICES *********************************************/
		//Used for tracking uniqueness when filling histograms, and for determining unused particles
		
		//Step 0
		Int_t locBeamID = dComboBeamWrapper->Get_BeamID();

		/*********************************************** GET FOUR-MOMENTUM **********************************************/

		// Get P4's: //is kinfit if kinfit performed, else is measured
		//dTargetP4 is target p4
		//Step 0
		TLorentzVector locBeamP4 = dComboBeamWrapper->Get_P4();
		TLorentzVector locPositiveP4 = dPositiveWrapper->Get_P4();
		TLorentzVector locNegativeP4 = dNegativeWrapper->Get_P4();
		TLorentzVector locRecoilP4 = dRecoilWrapper->Get_P4();

		// Get Measured P4's:
		//Step 0
		TLorentzVector locBeamP4_Measured = dComboBeamWrapper->Get_P4_Measured();
		TLorentzVector locPositiveP4_Measured = dPositiveWrapper->Get_P4_Measured();
		TLorentzVector locNegativeP4_Measured = dNegativeWrapper->Get_P4_Measured();
		// Note: Recoil proton is missing - no measured P4 available

	  
	  	/********************************************* GET COMBO RF TIMING INFO *****************************************/
		TLorentzVector locBeam_X4 = dComboBeamWrapper->Get_X4();
		TLorentzVector locBeamX4_Measured = dComboBeamWrapper->Get_X4_Measured();
		Double_t locBunchPeriod = dAnalysisUtilities.Get_BeamBunchPeriod(Get_RunNumber());
		Double_t locDeltaT_RF = dAnalysisUtilities.Get_DeltaT_RF(Get_RunNumber(), locBeamX4_Measured, dComboWrapper);
		Int_t locRelBeamBucket = dAnalysisUtilities.Get_RelativeBeamBucket(Get_RunNumber(), locBeamX4_Measured, dComboWrapper); // 0 for in-time events, non-zero integer for out-of-time photons
		Int_t locNumOutOfTimeBunchesInTree = 2; //YOU need to specify this number
			//Number of out-of-time beam bunches in tree (on a single side, so that total number out-of-time bunches accepted is 2 times this number for left + right bunches)

		// Mark the first in-time combo encountered (chiSq-sorted) for this event.
		if((locRelBeamBucket == 0) && !foundBestInTimeCombo) {
			foundBestInTimeCombo = true;
			ThisComboIsBestChiSq = true;
		}

		// Bool_t locSkipNearestOutOfTimeBunch = true; // True: skip events from nearest out-of-time bunch on either side (recommended).
		// Int_t locNumOutOfTimeBunchesToUse = locSkipNearestOutOfTimeBunch ? locNumOutOfTimeBunchesInTree-1:locNumOutOfTimeBunchesInTree;
		// Double_t locAccidentalScalingFactor = dAnalysisUtilities.Get_AccidentalScalingFactor(Get_RunNumber(), locBeamP4.E(), dIsMC); // Ideal value would be 1, but deviations require added factor, which is different for data and MC.
		// Double_t locAccidentalScalingFactorError = dAnalysisUtilities.Get_AccidentalScalingFactorError(Get_RunNumber(), locBeamP4.E()); // Ideal value would be 1, but deviations observed, need added factor.
		// Double_t locHistAccidWeightFactor = locRelBeamBucket==0 ? 1 : -locAccidentalScalingFactor/(2*locNumOutOfTimeBunchesToUse) ; // Weight by 1 for in-time events, ScalingFactor*(1/NBunches) for out-of-time
		// if(locSkipNearestOutOfTimeBunch && abs(locRelBeamBucket)==1) { // Skip nearest out-of-time bunch: tails of in-time distribution also leak in
		// 	dComboWrapper->Set_IsComboCut(true);
		// 	continue;
		// }

		/******************************************** Accidental Subtraction ***************************************************/
		double AccWeight = 0.;
		double locRFTime = dComboWrapper->Get_RFTime_Measured();
		// time difference between tagger and RF (corrected for production vertex position relative to target center)
		double locBeamDeltaT = locBeamX4_Measured.T() - (locRFTime + (locBeamX4_Measured.Z() - dTargetCenter.Z())/29.9792458);
		if(fabs(locBeamDeltaT) < 0.5*4.008) { // prompt signal recieves a weight of 1
			AccWeight = 1.;
	    }
	    else { // accidentals recieve a weight of 1/(# RF bunches included in TTree) which is 4 in this case
			AccWeight = -1./(2*locNumOutOfTimeBunchesInTree); //need to have the total number of buckets. Make the histogram wider to find out how many
	    }

		/*
		 cout << " Tagger Accidentals: dTargetCenter=" << dTargetCenter.Z()
		 << " locRFTime=" << locRFTime
		 << " locBeamDeltaT=" << locBeamDeltaT
		 << " AccWeight=" << AccWeight << endl;

		 cout << " locBeam_X4_Measured="; locBeam_X4_Measured.Print();
		*/

	
		//cout << "accidental subtraction" << "\n";



	  	/********************************************* COMBINE FOUR-MOMENTUM ********************************************/

		// Combine 4-vectors
		TLorentzVector locMissingP4_Measured = locBeamP4_Measured + dTargetP4; // Start with beam + target (missing p4 before subtracting final state)
		locMissingP4_Measured -= locPositiveP4_Measured + locNegativeP4_Measured; // Subtract measured final state p4's to get missing p4
		

		TLorentzVector loc2TrackP4 = locPositiveP4 + locNegativeP4;
		TLorentzVector loc2TrackP4_Measured = locPositiveP4_Measured + locNegativeP4_Measured;

		double M2TrackKin = loc2TrackP4.M();  // Generic name for 2-track invariant mass

		double positive_phi = locPositiveP4.Phi();
		double negative_phi = locNegativeP4.Phi();

		//Missing Mass Squared
		double locMissingMassSquared_Measured = locMissingP4_Measured.M2();
		double locMissingMassSquared_Residual = locMissingMassSquared_Measured - dTargetP4.M2();
		double locMissingEnergy_Residual = locMissingP4_Measured.E() - dTargetP4.E();

		

		// FCAL variables
		double FCAL_Energy_plus = dPositiveWrapper->Get_Energy_FCAL();
		double FCAL_Energy_minus = dNegativeWrapper->Get_Energy_FCAL();
		double FCAL_Elasticity = (FCAL_Energy_plus + FCAL_Energy_minus)/locBeamP4.E();
		double TrackFCAL_DOCA_plus = dPositiveWrapper->Get_TrackFCAL_DOCA();
		double TrackFCAL_DOCA_minus = dNegativeWrapper->Get_TrackFCAL_DOCA();
		double FCAL_E1E9_plus = dPositiveWrapper->Get_E1E9_FCAL();
		double FCAL_E1E9_minus = dNegativeWrapper->Get_E1E9_FCAL();
		double FCAL_E9E25_plus = dPositiveWrapper->Get_E9E25_FCAL();
		double FCAL_E9E25_minus = dNegativeWrapper->Get_E9E25_FCAL();
		double FCAL_kin_res_plus = (FCAL_Energy_plus - locPositiveP4.E())/locPositiveP4.E();
		double FCAL_kin_res_minus = (FCAL_Energy_minus - locNegativeP4.E())/locNegativeP4.E();
		double FCAL_meas_res_plus = (FCAL_Energy_plus - locPositiveP4_Measured.E())/locPositiveP4_Measured.E();
		double FCAL_meas_res_minus = (FCAL_Energy_minus - locNegativeP4_Measured.E())/locNegativeP4_Measured.E();
		// Note: The official FCAL_EoverP variable is declared later after potentially swapping positive/negative tracks if we're producing a tree to use in the training of neural nets. We write FCAL_Energy_plus/locPositiveP4.Vect().Mag() here for the histograms we make now.
		double FCAL_EoverPkin_plus_hist = FCAL_Energy_plus/locPositiveP4.Vect().Mag();
		double FCAL_EoverPkin_minus_hist = FCAL_Energy_minus/locNegativeP4.Vect().Mag();
		double FCAL_EoverPmeas_plus_hist = FCAL_Energy_plus/locPositiveP4_Measured.Vect().Mag();
		double FCAL_EoverPmeas_minus_hist = FCAL_Energy_minus/locNegativeP4_Measured.Vect().Mag();
		double FCAL_Emax_plus = FCAL_Energy_plus * FCAL_E9E25_plus * FCAL_E1E9_plus;
        double FCAL_Emax_minus = FCAL_Energy_minus * FCAL_E9E25_minus * FCAL_E1E9_minus;
        double FCAL_Asymmetry = fabs((FCAL_Energy_plus - FCAL_Energy_minus)/(FCAL_Energy_plus + FCAL_Energy_minus));
		// Geometric asymmetry from eq. 2.4: A_uv = |σ²_u - σ²_v| / |σ²_u + σ²_v|
		// SumU and SumV correspond to the second moments σ²_u and σ²_v
		double SumU_plus = dPositiveWrapper->Get_SumU_FCAL();
		double SumV_plus = dPositiveWrapper->Get_SumV_FCAL();
		double SumU_minus = dNegativeWrapper->Get_SumU_FCAL();
		double SumV_minus = dNegativeWrapper->Get_SumV_FCAL();
		double FCAL_UV_Asymmetry_plus = fabs(SumU_plus - SumV_plus) / fabs(SumU_plus + SumV_plus);
		double FCAL_UV_Asymmetry_minus = fabs(SumU_minus - SumV_minus) / fabs(SumU_minus + SumV_minus);

		// FCAL post-correction observables are always computed for MC side-by-side
		// with raw values (raw inputs are not overwritten).
		double FCAL_Energy_plus_postCorr = FCAL_Energy_plus;
		double FCAL_Energy_minus_postCorr = FCAL_Energy_minus;
		if(dIsMC) {
			const double pPlus = locPositiveP4.Vect().Mag();
			const double pMinus = locNegativeP4.Vect().Mag();
			if(pPlus <= 3.60) {
				FCAL_Energy_plus_postCorr -= ep_FCAL_cor1->Eval(pPlus);
				FCAL_Energy_plus_postCorr += ep_FCAL_cor2data1->Eval(pPlus);
			} else {
				FCAL_Energy_plus_postCorr -= ep_FCAL_cor2->Eval(pPlus);
				FCAL_Energy_plus_postCorr += ep_FCAL_cor2data2->Eval(pPlus);
			}
			if(pMinus <= 3.60) {
				FCAL_Energy_minus_postCorr -= em_FCAL_cor1->Eval(pMinus);
				FCAL_Energy_minus_postCorr += em_FCAL_cor2data1->Eval(pMinus);
			} else {
				FCAL_Energy_minus_postCorr -= em_FCAL_cor2->Eval(pMinus);
				FCAL_Energy_minus_postCorr += em_FCAL_cor2data2->Eval(pMinus);
			}
		}

		const double pPlusKin = locPositiveP4.Vect().Mag();
		const double pMinusKin = locNegativeP4.Vect().Mag();
		const double pPlusMeas = locPositiveP4_Measured.Vect().Mag();
		const double pMinusMeas = locNegativeP4_Measured.Vect().Mag();
		double FCAL_EoverPkin_plus_postCorr = (pPlusKin > 0.0) ? FCAL_Energy_plus_postCorr/pPlusKin : 0.0;
		double FCAL_EoverPkin_minus_postCorr = (pMinusKin > 0.0) ? FCAL_Energy_minus_postCorr/pMinusKin : 0.0;
		double FCAL_EoverPmeas_plus_postCorr = (pPlusMeas > 0.0) ? FCAL_Energy_plus_postCorr/pPlusMeas : 0.0;
		double FCAL_EoverPmeas_minus_postCorr = (pMinusMeas > 0.0) ? FCAL_Energy_minus_postCorr/pMinusMeas : 0.0;
		double FCAL_Elasticity_postCorr = (FCAL_Energy_plus_postCorr + FCAL_Energy_minus_postCorr)/locBeamP4.E();
		double FCAL_Asymmetry_postCorr = fabs((FCAL_Energy_plus_postCorr - FCAL_Energy_minus_postCorr)/(FCAL_Energy_plus_postCorr + FCAL_Energy_minus_postCorr));
		double deltaEfcalKinfit_plus_postCorr = FCAL_Energy_plus_postCorr - locPositiveP4.E();
		double deltaEfcalKinfit_minus_postCorr = FCAL_Energy_minus_postCorr - locNegativeP4.E();
		double kinRes_plus_abs_postCorr = (locPositiveP4.E() != 0.0) ? fabs(deltaEfcalKinfit_plus_postCorr/locPositiveP4.E()) : 0.0;
		double kinRes_minus_abs_postCorr = (locNegativeP4.E() != 0.0) ? fabs(deltaEfcalKinfit_minus_postCorr/locNegativeP4.E()) : 0.0;
		double measRes_plus_abs_postCorr = (locPositiveP4_Measured.E() != 0.0) ? fabs((FCAL_Energy_plus_postCorr - locPositiveP4_Measured.E())/locPositiveP4_Measured.E()) : 0.0;
		double measRes_minus_abs_postCorr = (locNegativeP4_Measured.E() != 0.0) ? fabs((FCAL_Energy_minus_postCorr - locNegativeP4_Measured.E())/locNegativeP4_Measured.E()) : 0.0;
		double FCAL_Emax_plus_postCorr = FCAL_Energy_plus_postCorr * FCAL_E9E25_plus * FCAL_E1E9_plus;
		double FCAL_Emax_minus_postCorr = FCAL_Energy_minus_postCorr * FCAL_E9E25_minus * FCAL_E1E9_minus;


		// Additional kinematic variables for histogram filling
		double theta_plus_deg = locPositiveP4.Theta()*TMath::RadToDeg();
		double theta_minus_deg = locNegativeP4.Theta()*TMath::RadToDeg();
		double pkin_plus = locPositiveP4.Vect().Mag();
		double pkin_minus = locNegativeP4.Vect().Mag();
		double deltaEfcalKinfit_plus = FCAL_Energy_plus - locPositiveP4.E();
		double deltaEfcalKinfit_minus = FCAL_Energy_minus - locNegativeP4.E();
		double kinRes_plus_abs = fabs(FCAL_kin_res_plus);
		double kinRes_minus_abs = fabs(FCAL_kin_res_minus);
		double measRes_plus_abs = fabs(FCAL_meas_res_plus);
		double measRes_minus_abs = fabs(FCAL_meas_res_minus);

		// dE/dx variables
		double TOF_dEdx_plus = dPositiveWrapper->Get_dEdx_TOF();
		double TOF_dEdx_minus = dNegativeWrapper->Get_dEdx_TOF();
		double FDC_dEdx_plus = dPositiveWrapper->Get_dEdx_FDC();
		double FDC_dEdx_minus = dNegativeWrapper->Get_dEdx_FDC();

		Int_t NumUnusedTracks = dComboWrapper->Get_NumUnusedTracks();
		double Energy_UnusedShowers = dComboWrapper->Get_Energy_UnusedShowers();
		// cout << "DEBUG 4.7: Finished variable declarations, starting cuts" << endl;

		// Calculate FCAL E/P variables for MLP/BDT inputs
		Double_t FCAL_EoverPkin_plus = FCAL_Energy_plus/locPositiveP4.Vect().Mag();
		Double_t FCAL_EoverPkin_minus = FCAL_Energy_minus/locNegativeP4.Vect().Mag();
		Double_t FCAL_EoverPmeas_plus = FCAL_Energy_plus/locPositiveP4_Measured.Vect().Mag();
		Double_t FCAL_EoverPmeas_minus = FCAL_Energy_minus/locNegativeP4_Measured.Vect().Mag();

		// TMVA input policy:
		//  - MC: use post-correction FCAL-derived observables
		//  - Data: use uncorrected observables
		const Double_t tmva_FCAL_EoverPkin_plus = dIsMC ? FCAL_EoverPkin_plus_postCorr : FCAL_EoverPkin_plus;
		const Double_t tmva_FCAL_EoverPkin_minus = dIsMC ? FCAL_EoverPkin_minus_postCorr : FCAL_EoverPkin_minus;
		const Double_t tmva_FCAL_Saturation_plus = dIsMC ? FCAL_Emax_plus_postCorr : FCAL_Emax_plus;
		const Double_t tmva_FCAL_Saturation_minus = dIsMC ? FCAL_Emax_minus_postCorr : FCAL_Emax_minus;

		// Declare variables for MLP/BDT responses (will be filled conditionally)
        Float_t mlp_response_minus = -999.0;  // Use sentinel value as default
        Float_t mlp_response_plus = -999.0;
		Float_t bdt_response_minus = -999.0;
		Float_t bdt_response_plus = -999.0;


		/******************************* STRUCT-DRIVEN COMBO CUTS ********************************/
		ComboCutInputs comboCutInputs;
		comboCutInputs.beamE = locBeamP4.E();
		comboCutInputs.ThisComboIsBestChiSq = ThisComboIsBestChiSq;
		comboCutInputs.theta1 = locPositiveP4.Theta()*TMath::RadToDeg();
		comboCutInputs.theta2 = locNegativeP4.Theta()*TMath::RadToDeg();
		comboCutInputs.Wkin = M2TrackKin;
		Double_t locKinFitChiSq = dComboWrapper->Get_ChiSq_KinFit("");
		Double_t locKinFitNDF = static_cast<Double_t>(dComboWrapper->Get_NDF_KinFit(""));
		comboCutInputs.chisqndf = (locKinFitNDF > 0.0) ? (locKinFitChiSq/locKinFitNDF) : 1.0e9;
		comboCutInputs.kinFitCL = dComboWrapper->Get_ConfidenceLevel_KinFit("");
		comboCutInputs.vertexZ = locBeamX4_Measured.Z();
		comboCutInputs.pMagPlus = locPositiveP4.Vect().Mag();
		comboCutInputs.pMagMinus = locNegativeP4.Vect().Mag();
		comboCutInputs.numUnusedTracks = NumUnusedTracks;
		comboCutInputs.energyUnusedShowers = Energy_UnusedShowers;
		comboCutInputs.FCAL_Energy_plus = FCAL_Energy_plus;
		comboCutInputs.FCAL_Energy_minus = FCAL_Energy_minus;
		comboCutInputs.TOF_dEdx_plus = TOF_dEdx_plus;
		comboCutInputs.TOF_dEdx_minus = TOF_dEdx_minus;
		comboCutInputs.FCAL_EoverPkin_plus = FCAL_EoverPkin_plus;
		comboCutInputs.FCAL_EoverPkin_minus = FCAL_EoverPkin_minus;
		comboCutInputs.FCAL_EoverPmeas_plus = FCAL_EoverPmeas_plus;
		comboCutInputs.FCAL_EoverPmeas_minus = FCAL_EoverPmeas_minus;
		comboCutInputs.FCAL_Elasticity = FCAL_Elasticity;
		comboCutInputs.TrackFCAL_DOCA_plus = TrackFCAL_DOCA_plus;
		comboCutInputs.TrackFCAL_DOCA_minus = TrackFCAL_DOCA_minus;
		comboCutInputs.MLP1 = mlp_response_plus;
		comboCutInputs.MLP2 = mlp_response_minus;
		comboCutInputs.BDT1 = bdt_response_plus;
		comboCutInputs.BDT2 = bdt_response_minus;


		// Special pre-fiducial diagnostics: fill before low-q2 exclusivity prefilter.
		if(dRuntimeFillSwitches.fillPreFidDiagnostics && ThisComboIsBestChiSq) {
			const bool keepUpToOneUnusedTrack = (NumUnusedTracks <= 1);
			if(keepUpToOneUnusedTrack) {
				dHist_MissingMassSquared_Measured_AllowOneUnused->Fill(locMissingMassSquared_Residual);
				dHist_preFid_MM2Residual_vs_RecoilP_UpTo1UnusedTrack->Fill(
					locRecoilP4.Vect().Mag(), locMissingMassSquared_Residual);
			}
		}

		// Evaluate MLP/BDT before low-q2 exclusivity prefilter so prefilter N-1 studies
		// can still keep ee MVA selection active. Sanitize all TMVA inputs and outputs.
		if(ApplyMLPClassification && dTMVAReader_MLP != NULL) {
			// ==================== Evaluate ELECTRON (negative track) ====================
			dTMVA_EoverPkin_FCAL = sanitize_f((Float_t)tmva_FCAL_EoverPkin_minus);
			dTMVA_TrackFCAL_DOCA = sanitize_f((Float_t)TrackFCAL_DOCA_minus);
			dTMVA_E9E25_FCAL = sanitize_f((Float_t)FCAL_E9E25_minus);
			dTMVA_E1E9_FCAL = sanitize_f((Float_t)FCAL_E1E9_minus);
			dTMVA_SumU_FCAL = sanitize_f((Float_t)SumU_minus);
			dTMVA_SumV_FCAL = sanitize_f((Float_t)SumV_minus);
			dTMVA_UV_Asymmetry_FCAL = sanitize_f((Float_t)FCAL_UV_Asymmetry_minus);
			dTMVA_Saturation_FCAL = sanitize_f((Float_t)tmva_FCAL_Saturation_minus);

			mlp_response_minus = (Float_t)dTMVAReader_MLP->EvaluateMVA("MLP::MLP");

			// ==================== Now evaluate POSITRON (positive track) ====================
			dTMVA_EoverPkin_FCAL = sanitize_f((Float_t)tmva_FCAL_EoverPkin_plus);
			dTMVA_TrackFCAL_DOCA = sanitize_f((Float_t)TrackFCAL_DOCA_plus);
			dTMVA_E9E25_FCAL = sanitize_f((Float_t)FCAL_E9E25_plus);
			dTMVA_E1E9_FCAL = sanitize_f((Float_t)FCAL_E1E9_plus);
			dTMVA_SumU_FCAL = sanitize_f((Float_t)SumU_plus);
			dTMVA_SumV_FCAL = sanitize_f((Float_t)SumV_plus);
			dTMVA_UV_Asymmetry_FCAL = sanitize_f((Float_t)FCAL_UV_Asymmetry_plus);
			dTMVA_Saturation_FCAL = sanitize_f((Float_t)tmva_FCAL_Saturation_plus);

			mlp_response_plus = (Float_t)dTMVAReader_MLP->EvaluateMVA("MLP::MLP");
		}

		if(ApplyBDTClassification && dTMVAReader_BDT != NULL) {
			// ==================== Evaluate ELECTRON (negative track) ====================
			dTMVA_EoverPkin_FCAL = sanitize_f((Float_t)tmva_FCAL_EoverPkin_minus);
			dTMVA_TrackFCAL_DOCA = sanitize_f((Float_t)TrackFCAL_DOCA_minus);
			dTMVA_E9E25_FCAL = sanitize_f((Float_t)FCAL_E9E25_minus);
			dTMVA_E1E9_FCAL = sanitize_f((Float_t)FCAL_E1E9_minus);
			dTMVA_SumU_FCAL = sanitize_f((Float_t)SumU_minus);
			dTMVA_SumV_FCAL = sanitize_f((Float_t)SumV_minus);
			dTMVA_UV_Asymmetry_FCAL = sanitize_f((Float_t)FCAL_UV_Asymmetry_minus);
			dTMVA_Saturation_FCAL = sanitize_f((Float_t)tmva_FCAL_Saturation_minus);

			bdt_response_minus = (Float_t)dTMVAReader_BDT->EvaluateMVA("BDT");

			// ==================== Now evaluate POSITRON (positive track) ====================
			dTMVA_EoverPkin_FCAL = sanitize_f((Float_t)tmva_FCAL_EoverPkin_plus);
			dTMVA_TrackFCAL_DOCA = sanitize_f((Float_t)TrackFCAL_DOCA_plus);
			dTMVA_E9E25_FCAL = sanitize_f((Float_t)FCAL_E9E25_plus);
			dTMVA_E1E9_FCAL = sanitize_f((Float_t)FCAL_E1E9_plus);
			dTMVA_SumU_FCAL = sanitize_f((Float_t)SumU_plus);
			dTMVA_SumV_FCAL = sanitize_f((Float_t)SumV_plus);
			dTMVA_UV_Asymmetry_FCAL = sanitize_f((Float_t)FCAL_UV_Asymmetry_plus);
			dTMVA_Saturation_FCAL = sanitize_f((Float_t)tmva_FCAL_Saturation_plus);

			bdt_response_plus = (Float_t)dTMVAReader_BDT->EvaluateMVA("BDT");
		}

		comboCutInputs.MLP1 = mlp_response_plus;
		comboCutInputs.MLP2 = mlp_response_minus;
		comboCutInputs.BDT1 = bdt_response_plus;
		comboCutInputs.BDT2 = bdt_response_minus;

		// Low-q2 exclusivity prefilter:
		//  - No unused charged tracks
		//  - No unused neutral shower energy
		//  - Both tracks with non-zero FCAL energy
		//  - Both tracks with non-zero TOF dE/dx
		CutConditions lowQ2ExclusivityCuts = {};
		lowQ2ExclusivityCuts.applyNoUnusedTracksCut = dAnalysisCutConditions.applyNoUnusedTracksCut;
		lowQ2ExclusivityCuts.maxNumUnusedTracks = dAnalysisCutConditions.maxNumUnusedTracks;
		lowQ2ExclusivityCuts.applyNoUnusedShowersCut = dAnalysisCutConditions.applyNoUnusedShowersCut;
		lowQ2ExclusivityCuts.maxUnusedShowersEnergy = 0.0;
		lowQ2ExclusivityCuts.applyFCALEnergyNonZeroCut = dAnalysisCutConditions.applyFCALEnergyNonZeroCut;
		lowQ2ExclusivityCuts.applyTOFdEdxNonZeroCut = dAnalysisCutConditions.applyTOFdEdxNonZeroCut;
		if(!ComboPassesCuts(comboCutInputs, lowQ2ExclusivityCuts)) {
			dComboWrapper->Set_IsComboCut(true);
			continue;
		}

		// Fill MLP response histograms (before systematics analysis)
		if(ApplyMLPClassification){
			dHist_MLPResponsePlus_vs_MLPResponseMinus->Fill(mlp_response_minus, mlp_response_plus);
			dHist_MLPResponse_minus->Fill(mlp_response_minus);
			dHist_MLPResponse_plus->Fill(mlp_response_plus);
		}

		// cout << "DEBUG 4.8: Finished cuts, starting FCAL corrections" << endl;

	// -------------------------------------------------------------
	//  FCAL ENERGY CORRECTIONS
	// Raw and corrected FCAL observables are both available above.
	// Keep ApplyFCALEnergyCorrections as a future behavior switch without overwriting raw values.
	if(ApplyFCALEnergyCorrections == true){
		// Intentionally no in-place overwrite of FCAL_Energy_plus/minus.
		// Use *_postCorr companion variables where corrected observables are desired.
		(void)FCAL_Elasticity_postCorr;
		(void)FCAL_Asymmetry_postCorr;
	}

		// cout << "DEBUG 4.9: Starting Bethe-Heitler polarization calculation" << endl;
		/**************************** Bethe-Heitler Polarization ****************************/
		// Calculate qvec2 and J_T polarization observable
		double qvec2 = (locBeamP4.Vect() - locPositiveP4.Vect() - locNegativeP4.Vect()).Mag2();
		double qvec2_meas = (locBeamP4_Measured.Vect() - locPositiveP4_Measured.Vect() - locNegativeP4_Measured.Vect()).Mag2();
		// cout << "DEBUG 4.9.1: Calculated qvec2=" << qvec2 << ", qvec2_meas=" << qvec2_meas << endl;

		/***** POLARIZATION OBSERVABLE phi of vector J_T *****/
		//yes, it's the magnitude of the full momentum, not just the transverse momentum, that goes in the formula
		Double_t p1_mag = locPositiveP4.Vect().Mag(); 
		Double_t p2_mag = locNegativeP4.Vect().Mag();
		Double_t p1_theta = locPositiveP4.Theta();
		Double_t p2_theta = locNegativeP4.Theta();
		// cout << "DEBUG 4.9.2: p1_mag=" << p1_mag << ", p2_mag=" << p2_mag << ", p1_theta=" << p1_theta << ", p2_theta=" << p2_theta << endl;

		Double_t denom1 = locPositiveP4.E() - p1_mag * cos(p1_theta);
		Double_t denom2 = locNegativeP4.E() - p2_mag * cos(p2_theta);
		// cout << "DEBUG 4.9.3: denom1=" << denom1 << ", denom2=" << denom2 << endl;
		
		Double_t JTx = locPositiveP4.X() * 2*locNegativeP4.E()/denom1 + locNegativeP4.X() * 2*locPositiveP4.E()/denom2;
		// cout << "DEBUG 4.9.4: JTx=" << JTx << endl;
		
		Double_t JTy = locPositiveP4.Y() * 2*locNegativeP4.E()/denom1 + locNegativeP4.Y() * 2*locPositiveP4.E()/denom2;
		// cout << "DEBUG 4.9.5: JTy=" << JTy << endl;

		Double_t JTphi = atan2(JTy, JTx)*TMath::RadToDeg();
		// cout << "DEBUG 4.9.6: JTphi=" << JTphi << endl;

		Double_t p1_mag_meas = locPositiveP4_Measured.Vect().Mag();
		Double_t p2_mag_meas = locNegativeP4_Measured.Vect().Mag();
		Double_t p1_theta_meas = locPositiveP4_Measured.Theta();
		Double_t p2_theta_meas = locNegativeP4_Measured.Theta();
		// cout << "DEBUG 4.9.7: Starting measured JT calculation" << endl;

		Double_t denom1_meas = locPositiveP4_Measured.E() - p1_mag_meas * cos(p1_theta_meas);
		Double_t denom2_meas = locNegativeP4_Measured.E() - p2_mag_meas * cos(p2_theta_meas);
		// cout << "DEBUG 4.9.8: denom1_meas=" << denom1_meas << ", denom2_meas=" << denom2_meas << endl;
		
		Double_t JTx_meas = locPositiveP4_Measured.X() * 2*locNegativeP4_Measured.E()/denom1_meas +
			locNegativeP4_Measured.X() * 2*locPositiveP4_Measured.E()/denom2_meas;
		Double_t JTy_meas = locPositiveP4_Measured.Y() * 2*locNegativeP4_Measured.E()/denom1_meas +
			locNegativeP4_Measured.Y() * 2*locPositiveP4_Measured.E()/denom2_meas;	
		Double_t JTphi_meas = atan2(JTy_meas, JTx_meas)*TMath::RadToDeg();
		// cout << "DEBUG 4.9.9: JTphi_meas=" << JTphi_meas << endl;
	  	/********************************* EXECUTE ANALYSIS ACTIONS *********************************/

	  	// Loop through the analysis actions, executing them in order for the active particle combo
	  	//dAnalyzeCutActions->Perform_Action(); // Must be executed before Execute_Actions()
		// cout << "DEBUG 5.0: About to call Execute_Actions()" << endl;
	  	if(!Execute_Actions()) //if the active combo fails a cut, IsComboCutFlag automatically set
	    	continue;
		// cout << "DEBUG 5.1: Execute_Actions() completed successfully" << endl;

	  	//if you manually execute any actions, and it fails a cut, be sure to call:
	  	//dComboWrapper->Set_IsComboCut(true);


		// Fill all diagnostic histograms now
		// cout << "DEBUG 5.2: About to fill diagnostic histograms" << endl;
		if(dRuntimeFillSwitches.fillPreFidDiagnostics) {
			dHist_BeamEnergy->Fill(locBeamP4.E());
			// cout << "DEBUG 5.2.1: Filled BeamEnergy" << endl;
			dHist_TaggerAccidentals->Fill(locBeamDeltaT);
			// cout << "DEBUG 5.2.2: Filled TaggerAccidentals" << endl;
			if(ThisComboIsBestChiSq){
				dHist_BeamEnergy_BestChiSq->Fill(locBeamP4.E());
				dHist_Wepem->Fill(M2TrackKin);
				dHist_MissingMassSquared_Measured->Fill(locMissingMassSquared_Residual);
				dHist_MissingEnergy_Measured->Fill(locMissingEnergy_Residual);
				dHist_VertexZ_BestChiSq->Fill(comboCutInputs.vertexZ);
				dHist_EoverP_measured_plus->Fill(FCAL_EoverPmeas_plus_hist);
				dHist_EoverP_measured_minus->Fill(FCAL_EoverPmeas_minus_hist);
				dHist_MLPResponse_plus->Fill(mlp_response_plus);
				dHist_MLPResponse_minus->Fill(mlp_response_minus);
				dHist_MLPResponsePlus_vs_MLPResponseMinus->Fill(mlp_response_plus, mlp_response_minus);			
			if(locKinFitNDF > 0.0){
				dHist_KinFitChiSq->Fill(locKinFitChiSq/locKinFitNDF);
			}
			dHist_KinFitCL->Fill(dComboWrapper->Get_ConfidenceLevel_KinFit(""));
			}
		}

		if(!kRunPostPreFidBlocks)
			continue;

		auto buildCategoryCuts = [&](const CutConditions& baseCuts, int categoryIndex) {
			CutConditions categoryCuts = baseCuts;
			switch(categoryIndex) {
				case kNoMVA_cutsBased:
					categoryCuts.applyMVACuts = false;
					categoryCuts.modelChoice = "none";
					break;
				case kPreMVASelection:
					categoryCuts.applyMVACuts = false;
					categoryCuts.modelChoice = "none";
					categoryCuts.applyPreselectionEoPCut = true;
					break;
				case kMLP_ee:
					categoryCuts.applyMVACuts = true;
					categoryCuts.modelChoice = "MLP";
					categoryCuts.particleChoice = "ee";
					break;
				case kMLP_pi:
					categoryCuts.applyMVACuts = true;
					categoryCuts.modelChoice = "MLP";
					categoryCuts.particleChoice = "pi";
					break;
				case kBDT_ee:
					categoryCuts.applyMVACuts = true;
					categoryCuts.modelChoice = "BDT";
					categoryCuts.particleChoice = "ee";
					break;
				case kBDT_pi:
					categoryCuts.applyMVACuts = true;
					categoryCuts.modelChoice = "BDT";
					categoryCuts.particleChoice = "pi";
					break;
				default:
					break;
			}
			return categoryCuts;
		};

		// Raw accidental-subtracted stream: fill for all surviving combos, not best-only.
		{
			auto fillRawAccSubdHistograms = [&](int beamWindowIndex, int categoryIndex) {
				CategoryHistogramSet& histSet = dRawAccSubdHistograms[beamWindowIndex][categoryIndex];
				histSet.BeamEnergy->Fill(locBeamP4.E(), AccWeight);
				histSet.RelBeamBucket->Fill(locRelBeamBucket);
				histSet.TaggerAccidentals->Fill(locBeamDeltaT, AccWeight);
				histSet.MissingMassSquared->Fill(locMissingMassSquared_Residual, AccWeight);
				histSet.MissingEnergy->Fill(locMissingEnergy_Residual, AccWeight);
				histSet.Wepem->Fill(M2TrackKin, AccWeight);
				histSet.qvec2_varWidth->Fill(qvec2, AccWeight);
				histSet.qvec2->Fill(qvec2, AccWeight);
				histSet.theta1->Fill(theta_plus_deg, AccWeight);
				histSet.theta2->Fill(theta_minus_deg, AccWeight);
				histSet.theta2_vs_theta1->Fill(theta_plus_deg, theta_minus_deg, AccWeight);
				histSet.ep_Pmag->Fill(pkin_plus, AccWeight);
				histSet.em_Pmag->Fill(pkin_minus, AccWeight);
				histSet.RecoilThetaVsP->Fill(locRecoilP4.Vect().Mag(), locRecoilP4.Theta()*TMath::RadToDeg(), AccWeight);

				if(activeRunPeriodIndex >= 0 && activePolarizationIndex >= 0) {
					const int runPeriodIndex = activeRunPeriodIndex;
					const int polIndex = activePolarizationIndex;
					if (histSet.JTphi[runPeriodIndex][polIndex]) histSet.JTphi[runPeriodIndex][polIndex]->Fill(JTphi, AccWeight);
					if (histSet.ep_phi[runPeriodIndex][polIndex]) histSet.ep_phi[runPeriodIndex][polIndex]->Fill(locPositiveP4.Phi()*TMath::RadToDeg(), AccWeight);
					if (histSet.em_phi[runPeriodIndex][polIndex]) histSet.em_phi[runPeriodIndex][polIndex]->Fill(locNegativeP4.Phi()*TMath::RadToDeg(), AccWeight);
				}
			};

			CutConditions rawFullSpectrumCuts = activeFiducialConditions;
			rawFullSpectrumCuts.applyBestChiSqComboCut = false;
			rawFullSpectrumCuts.applyMVACuts = false;
			if (dIsCPPRunPeriod) {
				rawFullSpectrumCuts.minBeamE = 4.0;
				rawFullSpectrumCuts.maxBeamE = 7.6;
			} else {
				rawFullSpectrumCuts.minBeamE = 7.0;
				rawFullSpectrumCuts.maxBeamE = 11.4;
			}
			for(int categoryIndex = 0; categoryIndex < kNumBestChiSqCategories; ++categoryIndex) {
				CutConditions categoryCuts = buildCategoryCuts(rawFullSpectrumCuts, categoryIndex);
				if (ComboPassesCuts(comboCutInputs, categoryCuts)) {
					fillRawAccSubdHistograms(0, categoryIndex);
				}
			}

			CutConditions rawCoherentPeakCuts = activeFiducialConditions;
			rawCoherentPeakCuts.applyBestChiSqComboCut = false;
			rawCoherentPeakCuts.applyMVACuts = false;
			if (dIsCPPRunPeriod) {
				rawCoherentPeakCuts.minBeamE = 5.35;
				rawCoherentPeakCuts.maxBeamE = 5.75;
			} else {
				rawCoherentPeakCuts.minBeamE = 8.2;
				rawCoherentPeakCuts.maxBeamE = 8.8;
			}
			for(int categoryIndex = 0; categoryIndex < kNumBestChiSqCategories; ++categoryIndex) {
				CutConditions categoryCuts = buildCategoryCuts(rawCoherentPeakCuts, categoryIndex);
				if (ComboPassesCuts(comboCutInputs, categoryCuts)) {
					fillRawAccSubdHistograms(1, categoryIndex);
				}
			}
		}


		// Main Fiducial Scope: Fill histograms and perform systematics for combos that pass all cuts up to this point, which may include some or all of the fiducial cuts depending on how the cut conditions are structured.
		
		if(ThisComboIsBestChiSq){
			// cout << "DEBUG 5.3.2: This is the best in-time combo, about to call FillSystematics" << endl;
				
			// Fill systematics histograms for best chi2 combo only
			if(dRuntimeFillSwitches.fillSystematics) {
				const bool includeQ2ResSystematics = (dIsMC && foundThrownParticles);
				const Double_t q2kinRes = includeQ2ResSystematics ? (qvec2 - qvec2_thrown) : 0.0;
				FillSystematics(comboCutInputs, qvec2, JTphi, M2TrackKin, q2kinRes, includeQ2ResSystematics, activeRunPeriodIndex, activePolarizationIndex);
			}
			// cout << "DEBUG 5.3.3: FillSystematics completed" << endl;
				
			// Fill thrown histograms for MC - only for best combo
			// foundThrownParticles ensures we found both e+ and e- in thrown data
			// cout << "DEBUG 5.3.4: dIsMC=" << dIsMC << ", foundThrownParticles=" << foundThrownParticles << endl;
			if (dRuntimeFillSwitches.fillThrownResolution && dIsMC && foundThrownParticles) {
				// cout << "DEBUG 5.3.5: Entering MC thrown histogram section" << endl;
				// Fill Full Spectrum histograms (7.0-11.4 GeV for GlueX-I)
				CutConditions fullSpectrumCuts = activeFiducialConditions;
				fullSpectrumCuts.applyMVACuts = false;  // Turn off MVA cuts for MC thrown histograms
				if (dIsCPPRunPeriod) {
					fullSpectrumCuts.minBeamE = 4.0;
					fullSpectrumCuts.maxBeamE = 7.6;
				} else {
					fullSpectrumCuts.minBeamE = 7.0;
					fullSpectrumCuts.maxBeamE = 11.4;
				}
				
				// cout << "DEBUG 5.3.6: About to check fullSpectrumCuts" << endl;
				if (ComboPassesCuts(comboCutInputs, fullSpectrumCuts)) {
					// cout << "DEBUG 5.3.7: Passed fullSpectrumCuts, about to fill resolution histograms" << endl;
					int beamWindowIndex = 0; // FullSpectrum
					if (foundThrownProton && dHist_RecoilThetaVsP_Thrown[beamWindowIndex]) {
						dHist_RecoilThetaVsP_Thrown[beamWindowIndex]->Fill(pmag_recoil_thrown, theta_recoil_thrown);
					}
					// Resolution histograms - kinfit vs thrown
					// cout << "DEBUG 5.3.8: Filling qvec2_res_vs_q2kin" << endl;
					dHist_qvec2_res_vs_q2kin[beamWindowIndex]->Fill(qvec2, qvec2 - qvec2_thrown);
					// cout << "DEBUG 5.3.9: Filling theta_KinRes (positive)" << endl;
					dHist_theta_KinRes_vs_theta_Thrown[beamWindowIndex]->Fill(theta_plus_thrown, theta_plus_deg - theta_plus_thrown);
					// cout << "DEBUG 5.3.10: Filling theta_KinRes (negative)" << endl;
					dHist_theta_KinRes_vs_theta_Thrown[beamWindowIndex]->Fill(theta_minus_thrown, theta_minus_deg - theta_minus_thrown);
					dHist_ep_Pmag_KinRes_vs_ep_Pmag_Thrown[beamWindowIndex]->Fill(Pmag_plus_thrown, pkin_plus - Pmag_plus_thrown);
					dHist_em_Pmag_KinRes_vs_em_Pmag_Thrown[beamWindowIndex]->Fill(Pmag_minus_thrown, pkin_minus - Pmag_minus_thrown);
					dHist_Wepem_KinRes_vs_Wepem_Thrown[beamWindowIndex]->Fill(Wepem_thrown, M2TrackKin - Wepem_thrown);
					
					// Thrown distributions
					dHist_qvec2_varWidth_Thrown[beamWindowIndex]->Fill(qvec2_thrown);
					dHist_qvec2_Thrown[beamWindowIndex]->Fill(qvec2_thrown);
					dHist_theta1_Thrown[beamWindowIndex]->Fill(theta_plus_thrown);
					dHist_theta2_Thrown[beamWindowIndex]->Fill(theta_minus_thrown);
					dHist_ep_Pmag_Thrown[beamWindowIndex]->Fill(Pmag_plus_thrown);
					dHist_em_Pmag_Thrown[beamWindowIndex]->Fill(Pmag_minus_thrown);
					dHist_Wepem_Thrown[beamWindowIndex]->Fill(Wepem_thrown);
					
					// Kinematic resolutions - 1D
					dHist_ep_Pmag_KinRes[beamWindowIndex]->Fill(pkin_plus - Pmag_plus_thrown);
					dHist_em_Pmag_KinRes[beamWindowIndex]->Fill(pkin_minus - Pmag_minus_thrown);
					dHist_Wepem_KinRes[beamWindowIndex]->Fill(M2TrackKin - Wepem_thrown);
					dHist_phi_KinRes[beamWindowIndex]->Fill(locPositiveP4.Phi()*TMath::RadToDeg() - phi_plus_thrown);
					dHist_phi_KinRes[beamWindowIndex]->Fill(locNegativeP4.Phi()*TMath::RadToDeg() - phi_minus_thrown);
					
					// Polarization-dependent thrown histograms
					if(activeRunPeriodIndex >= 0 && activePolarizationIndex >= 0) {
						const int runPeriodIndex = activeRunPeriodIndex;
						const int polIndex = activePolarizationIndex;
						dHist_JTphi_Thrown[beamWindowIndex][runPeriodIndex][polIndex]->Fill(JTphi_thrown);
						dHist_JTphi_kinRes[beamWindowIndex][runPeriodIndex][polIndex]->Fill(JTphi - JTphi_thrown);
						dHist_ep_phi_Thrown[beamWindowIndex][runPeriodIndex][polIndex]->Fill(phi_plus_thrown);
						dHist_em_phi_Thrown[beamWindowIndex][runPeriodIndex][polIndex]->Fill(phi_minus_thrown);
					}
				}
					
				// Fill Coherent Peak histograms (8.2-8.8 GeV for GlueX-I)
				CutConditions coherentPeakCuts = activeFiducialConditions;
				coherentPeakCuts.applyMVACuts = false;  // Turn off MVA cuts for MC thrown histograms
				if (dIsCPPRunPeriod) {
					coherentPeakCuts.minBeamE = 5.35;
					coherentPeakCuts.maxBeamE = 5.75;
				} else {
					coherentPeakCuts.minBeamE = 8.2;
					coherentPeakCuts.maxBeamE = 8.8;
				}
					
				// cout << "DEBUG 5.3.14: About to call ComboPassesCuts for coherentPeakCuts" << endl;
				bool coherentPeakPassed = ComboPassesCuts(comboCutInputs, coherentPeakCuts);
				// cout << "DEBUG 5.3.14.5: ComboPassesCuts returned " << (coherentPeakPassed ? "TRUE" : "FALSE") << endl;
				
				if (coherentPeakPassed) {
					// cout << "DEBUG 5.3.15: Passed coherentPeakCuts" << endl;
					int beamWindowIndex = 1; // CoherentPeak
					if (foundThrownProton && dHist_RecoilThetaVsP_Thrown[beamWindowIndex]) {
						dHist_RecoilThetaVsP_Thrown[beamWindowIndex]->Fill(pmag_recoil_thrown, theta_recoil_thrown);
					}
					// Resolution histograms - kinfit vs thrown
					dHist_qvec2_res_vs_q2kin[beamWindowIndex]->Fill(qvec2, qvec2 - qvec2_thrown);
					dHist_theta_KinRes_vs_theta_Thrown[beamWindowIndex]->Fill(theta_plus_thrown, theta_plus_deg - theta_plus_thrown);
					dHist_theta_KinRes_vs_theta_Thrown[beamWindowIndex]->Fill(theta_minus_thrown, theta_minus_deg - theta_minus_thrown);
					dHist_ep_Pmag_KinRes_vs_ep_Pmag_Thrown[beamWindowIndex]->Fill(Pmag_plus_thrown, pkin_plus - Pmag_plus_thrown);
					dHist_em_Pmag_KinRes_vs_em_Pmag_Thrown[beamWindowIndex]->Fill(Pmag_minus_thrown, pkin_minus - Pmag_minus_thrown);
					dHist_Wepem_KinRes_vs_Wepem_Thrown[beamWindowIndex]->Fill(Wepem_thrown, M2TrackKin - Wepem_thrown);
					
					// Thrown distributions
					dHist_qvec2_varWidth_Thrown[beamWindowIndex]->Fill(qvec2_thrown);
					dHist_qvec2_Thrown[beamWindowIndex]->Fill(qvec2_thrown);
					dHist_theta1_Thrown[beamWindowIndex]->Fill(theta_plus_thrown);
					dHist_theta2_Thrown[beamWindowIndex]->Fill(theta_minus_thrown);
					dHist_ep_Pmag_Thrown[beamWindowIndex]->Fill(Pmag_plus_thrown);
					dHist_em_Pmag_Thrown[beamWindowIndex]->Fill(Pmag_minus_thrown);
					dHist_Wepem_Thrown[beamWindowIndex]->Fill(Wepem_thrown);
					
					// Kinematic resolutions - 1D
					dHist_ep_Pmag_KinRes[beamWindowIndex]->Fill(pkin_plus - Pmag_plus_thrown);
					dHist_em_Pmag_KinRes[beamWindowIndex]->Fill(pkin_minus - Pmag_minus_thrown);
					dHist_Wepem_KinRes[beamWindowIndex]->Fill(M2TrackKin - Wepem_thrown);
					dHist_phi_KinRes[beamWindowIndex]->Fill(locPositiveP4.Phi()*TMath::RadToDeg() - phi_plus_thrown);
					dHist_phi_KinRes[beamWindowIndex]->Fill(locNegativeP4.Phi()*TMath::RadToDeg() - phi_minus_thrown);
					
					// Polarization-dependent thrown histograms
					if(activeRunPeriodIndex >= 0 && activePolarizationIndex >= 0) {
						const int runPeriodIndex = activeRunPeriodIndex;
						const int polIndex = activePolarizationIndex;
						dHist_JTphi_Thrown[beamWindowIndex][runPeriodIndex][polIndex]->Fill(JTphi_thrown);
						dHist_JTphi_kinRes[beamWindowIndex][runPeriodIndex][polIndex]->Fill(JTphi - JTphi_thrown);
						dHist_ep_phi_Thrown[beamWindowIndex][runPeriodIndex][polIndex]->Fill(phi_plus_thrown);
						dHist_em_phi_Thrown[beamWindowIndex][runPeriodIndex][polIndex]->Fill(phi_minus_thrown);
					}
				}
					// cout << "DEBUG 5.3.16: After coherentPeakCuts check" << endl;
			}
		}

			if(ThisComboIsBestChiSq) {
				// ============================================================================
				// BEST CHI-SQUARED COMBO FILLS
				// ============================================================================
				
				// cout << "DEBUG 5.3.17: About to set up fcalDiagnosticCuts" << endl;
				// Set up cut conditions for FCAL diagnostic histograms
				// Apply fiducial cuts but turn off theta and momentum cuts
				// Use full spectrum beam energy (not coherent peak), MLP model, ee selection
				// cout << "DEBUG 5.3.18: Creating fcalDiagnosticCuts" << endl;
				CutConditions fcalDiagnosticCuts = activeFiducialConditions;
				// cout << "DEBUG 5.3.19: Assigned baseline fiducial conditions to fcalDiagnosticCuts" << endl;
				fcalDiagnosticCuts.applyMinThetaCut = false;
				fcalDiagnosticCuts.applyMaxThetaCut = false;
				fcalDiagnosticCuts.apply2DThetaCut = false;
				fcalDiagnosticCuts.applyMomentumRangeCut = false;
				fcalDiagnosticCuts.minBeamE = 7.0;  // Full spectrum
				fcalDiagnosticCuts.maxBeamE = 11.4; // Full spectrum
				// cout << "DEBUG 5.3.20: About to assign modelChoice" << endl;
				fcalDiagnosticCuts.modelChoice = "MLP";
				// cout << "DEBUG 5.3.21: About to assign particleChoice" << endl;
				fcalDiagnosticCuts.particleChoice = "ee";
				// cout << "DEBUG 5.3.22: Finished setting up fcalDiagnosticCuts, about to call ComboPassesCuts" << endl;
				
				bool passFCALDiagnosticCuts = ComboPassesCuts(comboCutInputs, fcalDiagnosticCuts);
				// cout << "DEBUG 5.3.23: ComboPassesCuts returned for fcalDiagnosticCuts" << endl;
			
				// Helper lambda to fill CategoryHistogramSet for a given beam window
				auto fillCategoryHistograms = [&](int beamWindowIndex, int categoryIndex) {
					CategoryHistogramSet& histSet = dBestChiSqHistograms[beamWindowIndex][categoryIndex];
				
				// Fill CategoryHistogramSet - Fiducial directory histograms
				histSet.BeamEnergy->Fill(locBeamP4.E());
				histSet.RelBeamBucket->Fill(locRelBeamBucket);
				histSet.TaggerAccidentals->Fill(locBeamDeltaT);
				histSet.MissingMassSquared->Fill(locMissingMassSquared_Residual);
				histSet.MissingEnergy->Fill(locMissingEnergy_Residual);
				histSet.Wepem->Fill(M2TrackKin);
				histSet.qvec2_varWidth->Fill(qvec2);
				histSet.qvec2->Fill(qvec2);
				histSet.theta1->Fill(theta_plus_deg);
				histSet.theta2->Fill(theta_minus_deg);
				histSet.theta2_vs_theta1->Fill(theta_plus_deg, theta_minus_deg);
				histSet.ep_Pmag->Fill(pkin_plus);
				histSet.em_Pmag->Fill(pkin_minus);
				histSet.RecoilThetaVsP->Fill(locRecoilP4.Vect().Mag(), locRecoilP4.Theta()*TMath::RadToDeg());
				
				// Fill polarization-dependent histograms
				if(activeRunPeriodIndex >= 0 && activePolarizationIndex >= 0) {
					const int runPeriodIndex = activeRunPeriodIndex;
					const int polIndex = activePolarizationIndex;
					if (histSet.JTphi[runPeriodIndex][polIndex]) histSet.JTphi[runPeriodIndex][polIndex]->Fill(JTphi);
					if (histSet.ep_phi[runPeriodIndex][polIndex]) histSet.ep_phi[runPeriodIndex][polIndex]->Fill(locPositiveP4.Phi()*TMath::RadToDeg());
					if (histSet.em_phi[runPeriodIndex][polIndex]) histSet.em_phi[runPeriodIndex][polIndex]->Fill(locNegativeP4.Phi()*TMath::RadToDeg());
				}
				
				// Fill CategoryHistogramSet - FCAL directory histograms (with cut guard)
				if (passFCALDiagnosticCuts) {
					if (histSet.FCAL_Energy_pip) histSet.FCAL_Energy_pip->Fill(FCAL_Energy_plus);
					if (histSet.FCAL_Energy_pim) histSet.FCAL_Energy_pim->Fill(FCAL_Energy_minus);
					if (histSet.FCAL_EoverP_pip) histSet.FCAL_EoverP_pip->Fill(FCAL_EoverPkin_plus_hist);
					if (histSet.FCAL_EoverP_pim) histSet.FCAL_EoverP_pim->Fill(FCAL_EoverPkin_minus_hist);
					if (histSet.FCAL_EoverPmeas_pip) histSet.FCAL_EoverPmeas_pip->Fill(FCAL_EoverPmeas_plus_hist);
					if (histSet.FCAL_EoverPmeas_pim) histSet.FCAL_EoverPmeas_pim->Fill(FCAL_EoverPmeas_minus_hist);
					if (histSet.Delta_Efcal_kinfitE_vs_kinPmag_plus) histSet.Delta_Efcal_kinfitE_vs_kinPmag_plus->Fill(pkin_plus, deltaEfcalKinfit_plus);
					if (histSet.Delta_Efcal_kinfitE_vs_kinPmag_minus) histSet.Delta_Efcal_kinfitE_vs_kinPmag_minus->Fill(pkin_minus, deltaEfcalKinfit_minus);
					if (histSet.FCAL_Elasticity) histSet.FCAL_Elasticity->Fill(FCAL_Elasticity);
					if (histSet.FCAL_Asymmetry) histSet.FCAL_Asymmetry->Fill(FCAL_Asymmetry);
					
					// Fill elasticity in asymmetry regions
					if(FCAL_Asymmetry < 0.2) {
						if (histSet.FCAL_Elasticity_Asym_L0pt2) histSet.FCAL_Elasticity_Asym_L0pt2->Fill(FCAL_Elasticity);
					} else if(FCAL_Asymmetry < 0.5) {
						if (histSet.FCAL_Elasticity_Asym_G0pt2_L0pt5) histSet.FCAL_Elasticity_Asym_G0pt2_L0pt5->Fill(FCAL_Elasticity);
					} else {
						if (histSet.FCAL_Elasticity_Asym_G0pt5) histSet.FCAL_Elasticity_Asym_G0pt5->Fill(FCAL_Elasticity);
					}
					
					if (histSet.TrackFCAL_DOCA_plus) histSet.TrackFCAL_DOCA_plus->Fill(TrackFCAL_DOCA_plus);
					if (histSet.TrackFCAL_DOCA_minus) histSet.TrackFCAL_DOCA_minus->Fill(TrackFCAL_DOCA_minus);
					if (histSet.FCAL_E1E9_plus) histSet.FCAL_E1E9_plus->Fill(FCAL_E1E9_plus);
					if (histSet.FCAL_E1E9_minus) histSet.FCAL_E1E9_minus->Fill(FCAL_E1E9_minus);
					if (histSet.FCAL_E9E25_plus) histSet.FCAL_E9E25_plus->Fill(FCAL_E9E25_plus);
					if (histSet.FCAL_E9E25_minus) histSet.FCAL_E9E25_minus->Fill(FCAL_E9E25_minus);
					if (histSet.FCAL_kin_res_plus) histSet.FCAL_kin_res_plus->Fill(kinRes_plus_abs);
					if (histSet.FCAL_kin_res_minus) histSet.FCAL_kin_res_minus->Fill(kinRes_minus_abs);
					if (histSet.FCAL_meas_res_plus) histSet.FCAL_meas_res_plus->Fill(measRes_plus_abs);
					if (histSet.FCAL_meas_res_minus) histSet.FCAL_meas_res_minus->Fill(measRes_minus_abs);
					if (histSet.FCAL_Saturation_vs_Eshower_plus) histSet.FCAL_Saturation_vs_Eshower_plus->Fill(FCAL_Energy_plus, FCAL_Emax_plus);
					if (histSet.FCAL_Saturation_vs_Eshower_minus) histSet.FCAL_Saturation_vs_Eshower_minus->Fill(FCAL_Energy_minus, FCAL_Emax_minus);
					if (histSet.FCAL_Saturation_plus) histSet.FCAL_Saturation_plus->Fill(FCAL_Emax_plus);
					if (histSet.FCAL_Saturation_minus) histSet.FCAL_Saturation_minus->Fill(FCAL_Emax_minus);
					if (histSet.FCAL_SumU_plus) histSet.FCAL_SumU_plus->Fill(SumU_plus);
					if (histSet.FCAL_SumU_minus) histSet.FCAL_SumU_minus->Fill(SumU_minus);
					if (histSet.FCAL_SumV_plus) histSet.FCAL_SumV_plus->Fill(SumV_plus);
					if (histSet.FCAL_SumV_minus) histSet.FCAL_SumV_minus->Fill(SumV_minus);
					if (histSet.FCAL_UV_Asymmetry_plus) histSet.FCAL_UV_Asymmetry_plus->Fill(FCAL_UV_Asymmetry_plus);
					if (histSet.FCAL_UV_Asymmetry_minus) histSet.FCAL_UV_Asymmetry_minus->Fill(FCAL_UV_Asymmetry_minus);
				}
				
				// Fill CategoryHistogramSet - Measured directory histograms
				histSet.qvec2_Meas->Fill(qvec2_meas);
				histSet.theta1_Meas->Fill(locPositiveP4_Measured.Theta()*TMath::RadToDeg());
				histSet.theta2_Meas->Fill(locNegativeP4_Measured.Theta()*TMath::RadToDeg());
				histSet.theta2_vs_theta1_Meas->Fill(locPositiveP4_Measured.Theta()*TMath::RadToDeg(), locNegativeP4_Measured.Theta()*TMath::RadToDeg());
				histSet.ep_Pmag_Meas->Fill(locPositiveP4_Measured.Vect().Mag());
				histSet.em_Pmag_Meas->Fill(locNegativeP4_Measured.Vect().Mag());
				
				if(activeRunPeriodIndex >= 0 && activePolarizationIndex >= 0) {
					const int runPeriodIndex = activeRunPeriodIndex;
					const int polIndex = activePolarizationIndex;
					if (histSet.JTphi_Meas[runPeriodIndex][polIndex]) histSet.JTphi_Meas[runPeriodIndex][polIndex]->Fill(JTphi_meas);
					if (histSet.ep_phi_Meas[runPeriodIndex][polIndex]) histSet.ep_phi_Meas[runPeriodIndex][polIndex]->Fill(locPositiveP4_Measured.Phi()*TMath::RadToDeg());
					if (histSet.em_phi_Meas[runPeriodIndex][polIndex]) histSet.em_phi_Meas[runPeriodIndex][polIndex]->Fill(locNegativeP4_Measured.Phi()*TMath::RadToDeg());
				}
			};

				// Fill Full Spectrum histograms
				CutConditions fullSpectrumCuts = activeFiducialConditions;
				fullSpectrumCuts.applyMVACuts = false;  // Turn off MVA cuts for diagnostic histograms
				if (dIsCPPRunPeriod) {
					fullSpectrumCuts.minBeamE = 4.0;
					fullSpectrumCuts.maxBeamE = 7.6;
				} else {
					fullSpectrumCuts.minBeamE = 7.0;
					fullSpectrumCuts.maxBeamE = 11.4;
				}
				if (dRuntimeFillSwitches.fillCategoryDirectories) {
					for(int categoryIndex = 0; categoryIndex < kNumBestChiSqCategories; ++categoryIndex) {
						CutConditions categoryCuts = buildCategoryCuts(fullSpectrumCuts, categoryIndex);
						if (ComboPassesCuts(comboCutInputs, categoryCuts))
							fillCategoryHistograms(0, categoryIndex); // FullSpectrum
					}
				}
			
			// Fill Coherent Peak histograms
			CutConditions coherentPeakCuts = activeFiducialConditions;
			coherentPeakCuts.applyMVACuts = false;  // Turn off MVA cuts for diagnostic histograms
			if (dIsCPPRunPeriod) {
				coherentPeakCuts.minBeamE = 5.35;
				coherentPeakCuts.maxBeamE = 5.75;
			} else {
				coherentPeakCuts.minBeamE = 8.2;
				coherentPeakCuts.maxBeamE = 8.8;
			}
			if (dRuntimeFillSwitches.fillCategoryDirectories) {
				for(int categoryIndex = 0; categoryIndex < kNumBestChiSqCategories; ++categoryIndex) {
					CutConditions categoryCuts = buildCategoryCuts(coherentPeakCuts, categoryIndex);
					if (ComboPassesCuts(comboCutInputs, categoryCuts))
						fillCategoryHistograms(1, categoryIndex); // CoherentPeak
				}
			}
			
			// ============================================================================
			// FCAL STUDY HISTOGRAMS (with guard conditions)
			// Standalone FCAL study block used to validate FCAL response across regimes
			// ============================================================================
			if (dRuntimeFillSwitches.fillFCALStudy && passFCALDiagnosticCuts) {
				dHist_BeamEnergy_BestChiSq->Fill(locBeamP4.E());
				dHist_Wepem_BestChiSq->Fill(M2TrackKin);
				
				// FCAL directory - Basic FCAL histograms
				dHist_FCAL_Energy_pip->Fill(FCAL_Energy_plus);
				dHist_FCAL_Energy_pim->Fill(FCAL_Energy_minus);
				dHist_FCAL_EoverP_pip->Fill(FCAL_EoverPkin_plus_hist);
				dHist_FCAL_EoverP_pim->Fill(FCAL_EoverPkin_minus_hist);
				dHist_FCAL_EoverPmeas_pip->Fill(FCAL_EoverPmeas_plus_hist);
				dHist_FCAL_EoverPmeas_pim->Fill(FCAL_EoverPmeas_minus_hist);
				dHist_Delta_Efcal_kinfitE_vs_kinPmag_plus->Fill(pkin_plus, deltaEfcalKinfit_plus);
				dHist_Delta_Efcal_kinfitE_vs_kinPmag_minus->Fill(pkin_minus, deltaEfcalKinfit_minus);
				dHist_FCAL_Elasticity->Fill(FCAL_Elasticity);
				dHist_FCAL_Asymmetry->Fill(FCAL_Asymmetry);
				
				// Fill elasticity in asymmetry regions
				if(FCAL_Asymmetry < 0.2) {
					dHist_FCAL_Elasticity_Asym_L0pt2->Fill(FCAL_Elasticity);
				} else if(FCAL_Asymmetry < 0.5) {
					dHist_FCAL_Elasticity_Asym_G0pt2_L0pt5->Fill(FCAL_Elasticity);
				} else {
					dHist_FCAL_Elasticity_Asym_G0pt5->Fill(FCAL_Elasticity);
				}
				
				dHist_TrackFCAL_DOCA_plus->Fill(TrackFCAL_DOCA_plus);
				dHist_TrackFCAL_DOCA_minus->Fill(TrackFCAL_DOCA_minus);
				dHist_FCAL_E1E9_plus->Fill(FCAL_E1E9_plus);
				dHist_FCAL_E1E9_minus->Fill(FCAL_E1E9_minus);
				dHist_FCAL_E9E25_plus->Fill(FCAL_E9E25_plus);
				dHist_FCAL_E9E25_minus->Fill(FCAL_E9E25_minus);
				dHist_FCAL_kin_res_plus->Fill(kinRes_plus_abs);
				dHist_FCAL_kin_res_minus->Fill(kinRes_minus_abs);
				dHist_FCAL_meas_res_plus->Fill(measRes_plus_abs);
				dHist_FCAL_meas_res_minus->Fill(measRes_minus_abs);
				dHist_FCAL_Saturation_vs_Eshower_plus->Fill(FCAL_Energy_plus, FCAL_Emax_plus);
				dHist_FCAL_Saturation_vs_Eshower_minus->Fill(FCAL_Energy_minus, FCAL_Emax_minus);
				dHist_FCAL_Saturation_plus->Fill(FCAL_Emax_plus);
				dHist_FCAL_Saturation_minus->Fill(FCAL_Emax_minus);
				dHist_FCAL_SumU_plus->Fill(SumU_plus);
				dHist_FCAL_SumU_minus->Fill(SumU_minus);
				dHist_FCAL_SumV_plus->Fill(SumV_plus);
				dHist_FCAL_SumV_minus->Fill(SumV_minus);
				dHist_FCAL_UV_Asymmetry_plus->Fill(FCAL_UV_Asymmetry_plus);
				dHist_FCAL_UV_Asymmetry_minus->Fill(FCAL_UV_Asymmetry_minus);

				if(dIsMC) {
					dHist_FCAL_Energy_pip_PostCorr->Fill(FCAL_Energy_plus_postCorr);
					dHist_FCAL_Energy_pim_PostCorr->Fill(FCAL_Energy_minus_postCorr);
					dHist_FCAL_EoverP_pip_PostCorr->Fill(FCAL_EoverPkin_plus_postCorr);
					dHist_FCAL_EoverP_pim_PostCorr->Fill(FCAL_EoverPkin_minus_postCorr);
					dHist_FCAL_EoverPmeas_pip_PostCorr->Fill(FCAL_EoverPmeas_plus_postCorr);
					dHist_FCAL_EoverPmeas_pim_PostCorr->Fill(FCAL_EoverPmeas_minus_postCorr);
					dHist_Delta_Efcal_kinfitE_vs_kinPmag_plus_PostCorr->Fill(pkin_plus, deltaEfcalKinfit_plus_postCorr);
					dHist_Delta_Efcal_kinfitE_vs_kinPmag_minus_PostCorr->Fill(pkin_minus, deltaEfcalKinfit_minus_postCorr);
					dHist_FCAL_Elasticity_PostCorr->Fill(FCAL_Elasticity_postCorr);
					dHist_FCAL_Asymmetry_PostCorr->Fill(FCAL_Asymmetry_postCorr);
					if(FCAL_Asymmetry_postCorr < 0.2) {
						dHist_FCAL_Elasticity_Asym_L0pt2_PostCorr->Fill(FCAL_Elasticity_postCorr);
					} else if(FCAL_Asymmetry_postCorr < 0.5) {
						dHist_FCAL_Elasticity_Asym_G0pt2_L0pt5_PostCorr->Fill(FCAL_Elasticity_postCorr);
					} else {
						dHist_FCAL_Elasticity_Asym_G0pt5_PostCorr->Fill(FCAL_Elasticity_postCorr);
					}
					dHist_FCAL_kin_res_plus_PostCorr->Fill(kinRes_plus_abs_postCorr);
					dHist_FCAL_kin_res_minus_PostCorr->Fill(kinRes_minus_abs_postCorr);
					dHist_FCAL_meas_res_plus_PostCorr->Fill(measRes_plus_abs_postCorr);
					dHist_FCAL_meas_res_minus_PostCorr->Fill(measRes_minus_abs_postCorr);
					dHist_FCAL_Saturation_vs_Eshower_plus_PostCorr->Fill(FCAL_Energy_plus_postCorr, FCAL_Emax_plus_postCorr);
					dHist_FCAL_Saturation_vs_Eshower_minus_PostCorr->Fill(FCAL_Energy_minus_postCorr, FCAL_Emax_minus_postCorr);
					dHist_FCAL_Saturation_plus_PostCorr->Fill(FCAL_Emax_plus_postCorr);
					dHist_FCAL_Saturation_minus_PostCorr->Fill(FCAL_Emax_minus_postCorr);
				}
				
				// FCALvsTheta directory
				dFCALvsTheta.Energy[kPlus]->Fill(theta_plus_deg, FCAL_Energy_plus);
				dFCALvsTheta.Energy[kMinus]->Fill(theta_minus_deg, FCAL_Energy_minus);
				dFCALvsTheta.EoverP[kPlus]->Fill(theta_plus_deg, FCAL_EoverPkin_plus_hist);
				dFCALvsTheta.EoverP[kMinus]->Fill(theta_minus_deg, FCAL_EoverPkin_minus_hist);
				dFCALvsTheta.EoverPmeas[kPlus]->Fill(theta_plus_deg, FCAL_EoverPmeas_plus_hist);
				dFCALvsTheta.EoverPmeas[kMinus]->Fill(theta_minus_deg, FCAL_EoverPmeas_minus_hist);
				dFCALvsTheta.DeltaEfcal_kinfitE[kPlus]->Fill(theta_plus_deg, deltaEfcalKinfit_plus);
				dFCALvsTheta.DeltaEfcal_kinfitE[kMinus]->Fill(theta_minus_deg, deltaEfcalKinfit_minus);
				dFCALvsTheta.TrackDOCA[kPlus]->Fill(theta_plus_deg, TrackFCAL_DOCA_plus);
				dFCALvsTheta.TrackDOCA[kMinus]->Fill(theta_minus_deg, TrackFCAL_DOCA_minus);
				dFCALvsTheta.E1E9[kPlus]->Fill(theta_plus_deg, FCAL_E1E9_plus);
				dFCALvsTheta.E1E9[kMinus]->Fill(theta_minus_deg, FCAL_E1E9_minus);
				dFCALvsTheta.E9E25[kPlus]->Fill(theta_plus_deg, FCAL_E9E25_plus);
				dFCALvsTheta.E9E25[kMinus]->Fill(theta_minus_deg, FCAL_E9E25_minus);
				dFCALvsTheta.KinRes[kPlus]->Fill(theta_plus_deg, kinRes_plus_abs);
				dFCALvsTheta.KinRes[kMinus]->Fill(theta_minus_deg, kinRes_minus_abs);
				dFCALvsTheta.Saturation[kPlus]->Fill(theta_plus_deg, FCAL_Emax_plus);
				dFCALvsTheta.Saturation[kMinus]->Fill(theta_minus_deg, FCAL_Emax_minus);
				dFCALvsTheta.SumU[kPlus]->Fill(theta_plus_deg, SumU_plus);
				dFCALvsTheta.SumU[kMinus]->Fill(theta_minus_deg, SumU_minus);
				dFCALvsTheta.SumV[kPlus]->Fill(theta_plus_deg, SumV_plus);
				dFCALvsTheta.SumV[kMinus]->Fill(theta_minus_deg, SumV_minus);
				dFCALvsTheta.UVAsymmetry[kPlus]->Fill(theta_plus_deg, FCAL_UV_Asymmetry_plus);
				dFCALvsTheta.UVAsymmetry[kMinus]->Fill(theta_minus_deg, FCAL_UV_Asymmetry_minus);
				dFCALvsTheta.MeasRes[kPlus]->Fill(theta_plus_deg, measRes_plus_abs);
				dFCALvsTheta.MeasRes[kMinus]->Fill(theta_minus_deg, measRes_minus_abs);
				if(dIsMC) {
					dFCALvsTheta.Energy_PostCorr[kPlus]->Fill(theta_plus_deg, FCAL_Energy_plus_postCorr);
					dFCALvsTheta.Energy_PostCorr[kMinus]->Fill(theta_minus_deg, FCAL_Energy_minus_postCorr);
					dFCALvsTheta.EoverP_PostCorr[kPlus]->Fill(theta_plus_deg, FCAL_EoverPkin_plus_postCorr);
					dFCALvsTheta.EoverP_PostCorr[kMinus]->Fill(theta_minus_deg, FCAL_EoverPkin_minus_postCorr);
					dFCALvsTheta.EoverPmeas_PostCorr[kPlus]->Fill(theta_plus_deg, FCAL_EoverPmeas_plus_postCorr);
					dFCALvsTheta.EoverPmeas_PostCorr[kMinus]->Fill(theta_minus_deg, FCAL_EoverPmeas_minus_postCorr);
					dFCALvsTheta.DeltaEfcal_kinfitE_PostCorr[kPlus]->Fill(theta_plus_deg, deltaEfcalKinfit_plus_postCorr);
					dFCALvsTheta.DeltaEfcal_kinfitE_PostCorr[kMinus]->Fill(theta_minus_deg, deltaEfcalKinfit_minus_postCorr);
					dFCALvsTheta.KinRes_PostCorr[kPlus]->Fill(theta_plus_deg, kinRes_plus_abs_postCorr);
					dFCALvsTheta.KinRes_PostCorr[kMinus]->Fill(theta_minus_deg, kinRes_minus_abs_postCorr);
					dFCALvsTheta.MeasRes_PostCorr[kPlus]->Fill(theta_plus_deg, measRes_plus_abs_postCorr);
					dFCALvsTheta.MeasRes_PostCorr[kMinus]->Fill(theta_minus_deg, measRes_minus_abs_postCorr);
					dFCALvsTheta.Saturation_PostCorr[kPlus]->Fill(theta_plus_deg, FCAL_Emax_plus_postCorr);
					dFCALvsTheta.Saturation_PostCorr[kMinus]->Fill(theta_minus_deg, FCAL_Emax_minus_postCorr);
				}
				
				// FCALvsPkin directory
				dFCALvsPkin.Energy[kPlus]->Fill(pkin_plus, FCAL_Energy_plus);
				dFCALvsPkin.Energy[kMinus]->Fill(pkin_minus, FCAL_Energy_minus);
				dFCALvsPkin.EoverP[kPlus]->Fill(pkin_plus, FCAL_EoverPkin_plus_hist);
				dFCALvsPkin.EoverP[kMinus]->Fill(pkin_minus, FCAL_EoverPkin_minus_hist);
				dFCALvsPkin.EoverPmeas[kPlus]->Fill(pkin_plus, FCAL_EoverPmeas_plus_hist);
				dFCALvsPkin.EoverPmeas[kMinus]->Fill(pkin_minus, FCAL_EoverPmeas_minus_hist);
				dFCALvsPkin.DeltaEfcal_kinfitE[kPlus]->Fill(pkin_plus, deltaEfcalKinfit_plus);
				dFCALvsPkin.DeltaEfcal_kinfitE[kMinus]->Fill(pkin_minus, deltaEfcalKinfit_minus);
				dFCALvsPkin.TrackDOCA[kPlus]->Fill(pkin_plus, TrackFCAL_DOCA_plus);
				dFCALvsPkin.TrackDOCA[kMinus]->Fill(pkin_minus, TrackFCAL_DOCA_minus);
				dFCALvsPkin.E1E9[kPlus]->Fill(pkin_plus, FCAL_E1E9_plus);
				dFCALvsPkin.E1E9[kMinus]->Fill(pkin_minus, FCAL_E1E9_minus);
				dFCALvsPkin.E9E25[kPlus]->Fill(pkin_plus, FCAL_E9E25_plus);
				dFCALvsPkin.E9E25[kMinus]->Fill(pkin_minus, FCAL_E9E25_minus);
				dFCALvsPkin.KinRes[kPlus]->Fill(pkin_plus, kinRes_plus_abs);
				dFCALvsPkin.KinRes[kMinus]->Fill(pkin_minus, kinRes_minus_abs);
				dFCALvsPkin.MeasRes[kPlus]->Fill(pkin_plus, measRes_plus_abs);
				dFCALvsPkin.MeasRes[kMinus]->Fill(pkin_minus, measRes_minus_abs);
				dFCALvsPkin.Saturation[kPlus]->Fill(pkin_plus, FCAL_Emax_plus);
				dFCALvsPkin.Saturation[kMinus]->Fill(pkin_minus, FCAL_Emax_minus);
				dFCALvsPkin.SumU[kPlus]->Fill(pkin_plus, SumU_plus);
				dFCALvsPkin.SumU[kMinus]->Fill(pkin_minus, SumU_minus);
				dFCALvsPkin.SumV[kPlus]->Fill(pkin_plus, SumV_plus);
				dFCALvsPkin.SumV[kMinus]->Fill(pkin_minus, SumV_minus);
				dFCALvsPkin.UVAsymmetry[kPlus]->Fill(pkin_plus, FCAL_UV_Asymmetry_plus);
				dFCALvsPkin.UVAsymmetry[kMinus]->Fill(pkin_minus, FCAL_UV_Asymmetry_minus);
				if(dIsMC) {
					dFCALvsPkin.Energy_PostCorr[kPlus]->Fill(pkin_plus, FCAL_Energy_plus_postCorr);
					dFCALvsPkin.Energy_PostCorr[kMinus]->Fill(pkin_minus, FCAL_Energy_minus_postCorr);
					dFCALvsPkin.EoverP_PostCorr[kPlus]->Fill(pkin_plus, FCAL_EoverPkin_plus_postCorr);
					dFCALvsPkin.EoverP_PostCorr[kMinus]->Fill(pkin_minus, FCAL_EoverPkin_minus_postCorr);
					dFCALvsPkin.EoverPmeas_PostCorr[kPlus]->Fill(pkin_plus, FCAL_EoverPmeas_plus_postCorr);
					dFCALvsPkin.EoverPmeas_PostCorr[kMinus]->Fill(pkin_minus, FCAL_EoverPmeas_minus_postCorr);
					dFCALvsPkin.DeltaEfcal_kinfitE_PostCorr[kPlus]->Fill(pkin_plus, deltaEfcalKinfit_plus_postCorr);
					dFCALvsPkin.DeltaEfcal_kinfitE_PostCorr[kMinus]->Fill(pkin_minus, deltaEfcalKinfit_minus_postCorr);
					dFCALvsPkin.KinRes_PostCorr[kPlus]->Fill(pkin_plus, kinRes_plus_abs_postCorr);
					dFCALvsPkin.KinRes_PostCorr[kMinus]->Fill(pkin_minus, kinRes_minus_abs_postCorr);
					dFCALvsPkin.MeasRes_PostCorr[kPlus]->Fill(pkin_plus, measRes_plus_abs_postCorr);
					dFCALvsPkin.MeasRes_PostCorr[kMinus]->Fill(pkin_minus, measRes_minus_abs_postCorr);
					dFCALvsPkin.Saturation_PostCorr[kPlus]->Fill(pkin_plus, FCAL_Emax_plus_postCorr);
					dFCALvsPkin.Saturation_PostCorr[kMinus]->Fill(pkin_minus, FCAL_Emax_minus_postCorr);
				}

				// FCALvsqvec2 directory: fill plus/minus views directly; kBoth is composed in Finalize.
				dFCALvsqvec2.Energy[kPlus]->Fill(qvec2, FCAL_Energy_plus);
				dFCALvsqvec2.Energy[kMinus]->Fill(qvec2, FCAL_Energy_minus);
				dFCALvsqvec2.EoverP[kPlus]->Fill(qvec2, FCAL_EoverPkin_plus_hist);
				dFCALvsqvec2.EoverP[kMinus]->Fill(qvec2, FCAL_EoverPkin_minus_hist);
				dFCALvsqvec2.EoverPmeas[kPlus]->Fill(qvec2, FCAL_EoverPmeas_plus_hist);
				dFCALvsqvec2.EoverPmeas[kMinus]->Fill(qvec2, FCAL_EoverPmeas_minus_hist);
				dFCALvsqvec2.DeltaEfcal_kinfitE[kPlus]->Fill(qvec2, deltaEfcalKinfit_plus);
				dFCALvsqvec2.DeltaEfcal_kinfitE[kMinus]->Fill(qvec2, deltaEfcalKinfit_minus);
				dFCALvsqvec2.TrackDOCA[kPlus]->Fill(qvec2, TrackFCAL_DOCA_plus);
				dFCALvsqvec2.TrackDOCA[kMinus]->Fill(qvec2, TrackFCAL_DOCA_minus);
				dFCALvsqvec2.E1E9[kPlus]->Fill(qvec2, FCAL_E1E9_plus);
				dFCALvsqvec2.E1E9[kMinus]->Fill(qvec2, FCAL_E1E9_minus);
				dFCALvsqvec2.E9E25[kPlus]->Fill(qvec2, FCAL_E9E25_plus);
				dFCALvsqvec2.E9E25[kMinus]->Fill(qvec2, FCAL_E9E25_minus);
				dFCALvsqvec2.KinRes[kPlus]->Fill(qvec2, kinRes_plus_abs);
				dFCALvsqvec2.KinRes[kMinus]->Fill(qvec2, kinRes_minus_abs);
				dFCALvsqvec2.MeasRes[kPlus]->Fill(qvec2, measRes_plus_abs);
				dFCALvsqvec2.MeasRes[kMinus]->Fill(qvec2, measRes_minus_abs);
				dFCALvsqvec2.Saturation[kPlus]->Fill(qvec2, FCAL_Emax_plus);
				dFCALvsqvec2.Saturation[kMinus]->Fill(qvec2, FCAL_Emax_minus);
				dFCALvsqvec2.SumU[kPlus]->Fill(qvec2, SumU_plus);
				dFCALvsqvec2.SumU[kMinus]->Fill(qvec2, SumU_minus);
				dFCALvsqvec2.SumV[kPlus]->Fill(qvec2, SumV_plus);
				dFCALvsqvec2.SumV[kMinus]->Fill(qvec2, SumV_minus);
				dFCALvsqvec2.UVAsymmetry[kPlus]->Fill(qvec2, FCAL_UV_Asymmetry_plus);
				dFCALvsqvec2.UVAsymmetry[kMinus]->Fill(qvec2, FCAL_UV_Asymmetry_minus);

				if(dIsMC) {
					dFCALvsqvec2.Energy_PostCorr[kPlus]->Fill(qvec2, FCAL_Energy_plus_postCorr);
					dFCALvsqvec2.Energy_PostCorr[kMinus]->Fill(qvec2, FCAL_Energy_minus_postCorr);
					dFCALvsqvec2.EoverP_PostCorr[kPlus]->Fill(qvec2, FCAL_EoverPkin_plus_postCorr);
					dFCALvsqvec2.EoverP_PostCorr[kMinus]->Fill(qvec2, FCAL_EoverPkin_minus_postCorr);
					dFCALvsqvec2.EoverPmeas_PostCorr[kPlus]->Fill(qvec2, FCAL_EoverPmeas_plus_postCorr);
					dFCALvsqvec2.EoverPmeas_PostCorr[kMinus]->Fill(qvec2, FCAL_EoverPmeas_minus_postCorr);
					dFCALvsqvec2.DeltaEfcal_kinfitE_PostCorr[kPlus]->Fill(qvec2, deltaEfcalKinfit_plus_postCorr);
					dFCALvsqvec2.DeltaEfcal_kinfitE_PostCorr[kMinus]->Fill(qvec2, deltaEfcalKinfit_minus_postCorr);
					dFCALvsqvec2.KinRes_PostCorr[kPlus]->Fill(qvec2, kinRes_plus_abs_postCorr);
					dFCALvsqvec2.KinRes_PostCorr[kMinus]->Fill(qvec2, kinRes_minus_abs_postCorr);
					dFCALvsqvec2.MeasRes_PostCorr[kPlus]->Fill(qvec2, measRes_plus_abs_postCorr);
					dFCALvsqvec2.MeasRes_PostCorr[kMinus]->Fill(qvec2, measRes_minus_abs_postCorr);
					dFCALvsqvec2.Saturation_PostCorr[kPlus]->Fill(qvec2, FCAL_Emax_plus_postCorr);
					dFCALvsqvec2.Saturation_PostCorr[kMinus]->Fill(qvec2, FCAL_Emax_minus_postCorr);
				}
			}

			// ============================================================================
			// FID_Nminus1_Plots: Best in-time combo only
			// These are validation plots showing cut variable distributions with N-1 cuts applied
			// ============================================================================
			
			// Only fill N-1 plots for the best in-time chi-squared combo
			if (dRuntimeFillSwitches.fillNminus1 && ThisComboIsBestChiSq) {
				// Beam-energy N-1 should not be beam-window gated: disable beam cuts and
				// fill both window histograms so FullSpectrum/CoherentPeak views are identical.
				CutConditions fidN1_BeamE_all = activeFiducialConditions;
				fidN1_BeamE_all.applyMinBeamEcut = false;
				fidN1_BeamE_all.applyMaxBeamEcut = false;
				if (ComboPassesCuts(comboCutInputs, fidN1_BeamE_all)) {
					dHist_FID_Nminus1_BeamEnergy[kFullSpectrum]->Fill(comboCutInputs.beamE);
					dHist_FID_Nminus1_BeamEnergy[kCoherentPeak]->Fill(comboCutInputs.beamE);
				}

				// Helper lambda to fill N-1 histograms for a given beam window.
				// Use the caller-provided base cuts so FullSpectrum/CoherentPeak are treated correctly.
				auto fillNminus1Histograms = [&](int beamWindowIndex, const CutConditions& baseCuts) {
				// Standard N-1: All fiducial cuts ON except the one being plotted
				CutConditions fidN1_W = baseCuts;
				fidN1_W.applyMinWkinCut = false;
				fidN1_W.applyMaxWkinCut = false;
				if (ComboPassesCuts(comboCutInputs, fidN1_W)) {
					dHist_FID_Nminus1_Wepem[beamWindowIndex]->Fill(comboCutInputs.Wkin);
				}

				// Momentum N-1: turn off momentum-range cut only
				CutConditions fidN1_Pmag = baseCuts;
				fidN1_Pmag.applyMomentumRangeCut = false;
				if (ComboPassesCuts(comboCutInputs, fidN1_Pmag)) {
					dHist_FID_Nminus1_p_ep[beamWindowIndex]->Fill(comboCutInputs.pMagPlus);
					dHist_FID_Nminus1_p_em[beamWindowIndex]->Fill(comboCutInputs.pMagMinus);
				}

				// Theta N-1: Turn off min/max theta (NOT 2D) for theta1 and theta2
				CutConditions fidN1_Theta1 = baseCuts;
				fidN1_Theta1.applyMinThetaCut = false;
				fidN1_Theta1.applyMaxThetaCut = false;
				// Keep 2D cut ON
				if (ComboPassesCuts(comboCutInputs, fidN1_Theta1)) {
					dHist_FID_Nminus1_theta1[beamWindowIndex]->Fill(comboCutInputs.theta1);
				}

				CutConditions fidN1_Theta2 = baseCuts;
				fidN1_Theta2.applyMinThetaCut = false;
				fidN1_Theta2.applyMaxThetaCut = false;
				// Keep 2D cut ON
				if (ComboPassesCuts(comboCutInputs, fidN1_Theta2)) {
					dHist_FID_Nminus1_theta2[beamWindowIndex]->Fill(comboCutInputs.theta2);
				}

				// Theta2 vs Theta1 N-1: Turn off 2D theta cut ONLY
				CutConditions fidN1_Theta2D = baseCuts;
				fidN1_Theta2D.apply2DThetaCut = false;
				// Keep min/max cuts ON
				if (ComboPassesCuts(comboCutInputs, fidN1_Theta2D)) {
					dHist_FID_Nminus1_theta2_vs_theta1[beamWindowIndex]->Fill(comboCutInputs.theta1, comboCutInputs.theta2);
				}

				// Vertex-Z N-1: turn off vertex-Z cut only
				CutConditions fidN1_VertexZ = baseCuts;
				fidN1_VertexZ.applyVertexZCut = false;
				if (ComboPassesCuts(comboCutInputs, fidN1_VertexZ)) {
					dHist_FID_Nminus1_VertexZ[beamWindowIndex]->Fill(comboCutInputs.vertexZ);
				}

				// Theta variant N-1: Turn off ALL theta cuts (min, max, 2D)
				CutConditions fidN1_noThetaCuts = baseCuts;
				fidN1_noThetaCuts.applyMinThetaCut = false;
				fidN1_noThetaCuts.applyMaxThetaCut = false;
				fidN1_noThetaCuts.apply2DThetaCut = false;
				if (ComboPassesCuts(comboCutInputs, fidN1_noThetaCuts)) {
					dHist_FID_Nminus1_theta1_noThetaCuts[beamWindowIndex]->Fill(comboCutInputs.theta1);
					dHist_FID_Nminus1_theta2_noThetaCuts[beamWindowIndex]->Fill(comboCutInputs.theta2);
					dHist_FID_Nminus1_theta2_vs_theta1_noThetaCuts[beamWindowIndex]->Fill(comboCutInputs.theta1, comboCutInputs.theta2);
				}

				// Preselection E/p N-1: Turn off preselection E/p cut
				CutConditions fidN1_EoverP = baseCuts;
				fidN1_EoverP.applyPreselectionEoPCut = false;
				if (ComboPassesCuts(comboCutInputs, fidN1_EoverP)) {
					dHist_FID_Nminus1_EoverP_ep[beamWindowIndex]->Fill(comboCutInputs.FCAL_EoverPmeas_plus);
					dHist_FID_Nminus1_EoverP_em[beamWindowIndex]->Fill(comboCutInputs.FCAL_EoverPmeas_minus);
				}

				// Exclusivity N-1: Turn off unused tracks/showers cuts
				CutConditions fidN1_UnusedTracks = baseCuts;
				fidN1_UnusedTracks.applyNoUnusedTracksCut = false;
				if (ComboPassesCuts(comboCutInputs, fidN1_UnusedTracks)) {
					dHist_FID_Nminus1_NumUnusedTracks[beamWindowIndex]->Fill(comboCutInputs.numUnusedTracks);
				}

				CutConditions fidN1_UnusedShowers = baseCuts;
				fidN1_UnusedShowers.applyNoUnusedShowersCut = false;
				if (ComboPassesCuts(comboCutInputs, fidN1_UnusedShowers)) {
					dHist_FID_Nminus1_UnusedShowerEnergy[beamWindowIndex]->Fill(comboCutInputs.energyUnusedShowers);
				}

				// Forward detector N-1: Turn off FCAL/TOF cuts
				CutConditions fidN1_FCALEnergy = baseCuts;
				fidN1_FCALEnergy.applyFCALEnergyNonZeroCut = false;
				if (ComboPassesCuts(comboCutInputs, fidN1_FCALEnergy)) {
					dHist_FID_Nminus1_FCALEnergy_ep[beamWindowIndex]->Fill(comboCutInputs.FCAL_Energy_plus);
					dHist_FID_Nminus1_FCALEnergy_em[beamWindowIndex]->Fill(comboCutInputs.FCAL_Energy_minus);
				}

				CutConditions fidN1_TOFdEdx = baseCuts;
				fidN1_TOFdEdx.applyTOFdEdxNonZeroCut = false;
				if (ComboPassesCuts(comboCutInputs, fidN1_TOFdEdx)) {
					dHist_FID_Nminus1_TOFdEdx_ep[beamWindowIndex]->Fill(comboCutInputs.TOF_dEdx_plus);
					dHist_FID_Nminus1_TOFdEdx_em[beamWindowIndex]->Fill(comboCutInputs.TOF_dEdx_minus);
				}

				// MVA N-1: Show MLP/BDT distributions with MVA cuts OFF (all other fidcuts ON)
				CutConditions fidN1_MLP = baseCuts;
				fidN1_MLP.applyMVACuts = false;
				if (ComboPassesCuts(comboCutInputs, fidN1_MLP)) {
					dHist_FID_Nminus1_MLP_ep[beamWindowIndex]->Fill(comboCutInputs.MLP1);
					dHist_FID_Nminus1_MLP_em[beamWindowIndex]->Fill(comboCutInputs.MLP2);
					dHist_FID_Nminus1_MLP_ep_vs_em[beamWindowIndex]->Fill(comboCutInputs.MLP2, comboCutInputs.MLP1);
				}

				CutConditions fidN1_BDT = baseCuts;
				fidN1_BDT.modelChoice = "BDT";
				fidN1_BDT.applyMVACuts = false;
				if (ComboPassesCuts(comboCutInputs, fidN1_BDT)) {
					dHist_FID_Nminus1_BDT_ep[beamWindowIndex]->Fill(comboCutInputs.BDT1);
					dHist_FID_Nminus1_BDT_em[beamWindowIndex]->Fill(comboCutInputs.BDT2);
					dHist_FID_Nminus1_BDT_ep_vs_em[beamWindowIndex]->Fill(comboCutInputs.BDT2, comboCutInputs.BDT1);
				}

				// Special "and_noMVA" N-1: Cuts-based fallback (modelChoice="none", FCAL cuts for pion rejection)
				CutConditions fidN1_EoverP_noMVA = baseCuts;
				fidN1_EoverP_noMVA.modelChoice = "none";
				fidN1_EoverP_noMVA.applyPreselectionEoPCut = false;
				// Keep FCAL Elasticity cuts ON
				if (ComboPassesCuts(comboCutInputs, fidN1_EoverP_noMVA)) {
					dHist_FID_Nminus1_EoverP_ep_and_noMVA[beamWindowIndex]->Fill(comboCutInputs.FCAL_EoverPmeas_plus);
				}

				CutConditions fidN1_FCALElasticity_noMVA = baseCuts;
				fidN1_FCALElasticity_noMVA.modelChoice = "none";
				fidN1_FCALElasticity_noMVA.applyMinFCALElasticityCut = false;
				fidN1_FCALElasticity_noMVA.applyMaxFCALElasticityCut = false;
				// Keep E/p cuts ON
				if (ComboPassesCuts(comboCutInputs, fidN1_FCALElasticity_noMVA)) {
					dHist_FID_Nminus1_FCALElasticity_and_noMVA[beamWindowIndex]->Fill(comboCutInputs.FCAL_Elasticity);
				}
				}; // End lambda
				
				// Fill Full Spectrum N-1 histograms
				CutConditions fullSpectrumCuts = activeFiducialConditions;
				if (dIsCPPRunPeriod) {
					fullSpectrumCuts.minBeamE = 4.0;
					fullSpectrumCuts.maxBeamE = 7.6;
				} else {
					fullSpectrumCuts.minBeamE = 7.0;
					fullSpectrumCuts.maxBeamE = 11.4;
				}
				if (comboCutInputs.beamE >= fullSpectrumCuts.minBeamE && comboCutInputs.beamE <= fullSpectrumCuts.maxBeamE) {
					fillNminus1Histograms(0, fullSpectrumCuts); // FullSpectrum
				}
				
				// Fill Coherent Peak N-1 histograms
				CutConditions coherentPeakCuts = activeFiducialConditions;
				if (dIsCPPRunPeriod) {
					coherentPeakCuts.minBeamE = 5.35;
					coherentPeakCuts.maxBeamE = 5.75;
				} 
				else {
					coherentPeakCuts.minBeamE = 8.2;
					coherentPeakCuts.maxBeamE = 8.8;
				}
				if (comboCutInputs.beamE >= coherentPeakCuts.minBeamE && comboCutInputs.beamE <= coherentPeakCuts.maxBeamE) {
					fillNminus1Histograms(1, coherentPeakCuts); // CoherentPeak
				}
		}
			}
		
	

		//E.g. Cut on missing mass squared if desired
		//if((locMissingMassSquared < -0.04) || (locMissingMassSquared > 0.04))
		//{
		//	dComboWrapper->Set_IsComboCut(true);
		//	continue;
		//}


		/****************************************** FILL FLAT TREE (IF DESIRED) ******************************************/
		// It is not desired.
		//Fill_FlatTree(); //for all surviving combos
		
	} // end of combo loop

	//FILL HISTOGRAMS: Num combos / events surviving actions
	//Fill_NumCombosSurvivedHists();


	/************************************ EXAMPLE: FILL CLONE OF TTREE HERE WITH CUTS APPLIED ************************************/
	/*
	Bool_t locIsEventCut = true;
	for(UInt_t loc_i = 0; loc_i < Get_NumCombos(); ++loc_i) {
		//Set branch array indices for combo and all combo particles
		dComboWrapper->Set_ComboIndex(loc_i);
		// Is used to indicate when combos have been cut
		if(dComboWrapper->Get_IsComboCut())
			continue;
		locIsEventCut = false; // At least one combo succeeded
		break;
	}
	if(!locIsEventCut && dOutputTreeFileName != ""){
	  Fill_OutputTree();
		
	}
	*/

	return kTRUE;
}

void DSelector_2eMissingProton_Systematics::BookSystematics(void)
{
	// Ensure newly created TH1/TH2 objects are owned by the current TDirectory
	// so they are persisted in the output ROOT file.
	TH1::AddDirectory(kTRUE);

	TDirectory* mainDir = gDirectory;
	TDirectory* systDir = dDir_Systematics ? dDir_Systematics : mainDir->mkdir("systematics");
	TDirectory* phiJTDir = dDir_Systematics_JTphi ? dDir_Systematics_JTphi : systDir->mkdir("JTphi");
	TDirectory* q2Dir = dDir_Systematics_q2 ? dDir_Systematics_q2 : systDir->mkdir("q2");
	TDirectory* thetaSysDir = dDir_Systematics_theta ? dDir_Systematics_theta : systDir->mkdir("theta_sys");
	TDirectory* invmassSysDir = dDir_Systematics_invmass ? dDir_Systematics_invmass : systDir->mkdir("invmass_sys");
	TDirectory* chisqSysDir = dDir_Systematics_chisq ? dDir_Systematics_chisq : systDir->mkdir("chisq_sys");
	TDirectory* thrownDir = nullptr;
	if(dIsMC)
		thrownDir = dDir_Systematics_thrown ? dDir_Systematics_thrown : systDir->mkdir("thrown");
	TDirectory* mvaResponseDir = dDir_MVAResponsePlots ? dDir_MVAResponsePlots : mainDir->mkdir("MVA_response_plots");

	dMethods[0] = "MLP";
	dMethods[1] = "BDT";
	dPidChoices[0] = "ee";
	dPidChoices[1] = "pi";
	dPidChoices[2] = "none";

	const double minAngles[kNumAngles] = {0.0, 0.5, 1.0, 1.5, 1.75, 2.0};
	for (int i = 0; i < kNumAngles; ++i) {
		dMinAngles[i] = minAngles[i];
	}

	const double beamEnergyRegions[kNumBeamRegions][2] = {
		{7.8, 8.0}, {8.0, 8.2}, {8.2, 8.4}, {8.4, 8.6},
		{8.6, 8.8}, {8.8, 9.0}, {9.0, 9.2}, {9.2, 9.4},
		{9.4, 9.6}, {9.6, 9.8}, {9.8, 10.0}, {10.0, 10.2},
		{10.2, 10.4}, {10.4, 10.6}, {10.6, 10.8}, {10.8, 11.0}
	};
	for (int i = 0; i < kNumBeamRegions; ++i) {
		dBeamEnergyRegions[i][0] = beamEnergyRegions[i][0];
		dBeamEnergyRegions[i][1] = beamEnergyRegions[i][1];
	}

	const double massRegions[kNumMassRegions][2] = {
		{0.250, 0.275}, {0.275, 0.300}, {0.300, 0.325}, {0.325, 0.350},
		{0.350, 0.375}, {0.375, 0.400}, {0.400, 0.45}, {0.45, 0.50},
		{0.50, 0.55}, {0.55, 0.60}
	};
	for (int i = 0; i < kNumMassRegions; ++i) {
		dMassRegions[i][0] = massRegions[i][0];
		dMassRegions[i][1] = massRegions[i][1];
	}

	const double minInvMass[kNumWminCuts] = {0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4};
	for (int i = 0; i < kNumWminCuts; ++i) {
		dMinInvMass[i] = minInvMass[i];
	}

	const double maxChiSq[kNumChiSqCuts] = {15., 14., 13., 12., 11., 10., 9., 8., 7., 6., 5., 4., 3., 2., 1., 0.5};
	for (int i = 0; i < kNumChiSqCuts; ++i) {
		dMaxChiSq[i] = maxChiSq[i];
	}

	const float mlpEe[kNumMLPeeThresholds] = {0.55f, 0.60f, 0.65f, 0.70f, 0.75f, 0.80f, 0.85f, 0.90f, 0.95f};
	const float mlpPi[kNumMLPpiThresholds] = {0.45f, 0.40f, 0.35f, 0.30f, 0.25f, 0.20f, 0.15f, 0.10f, 0.05f};
	for (int i = 0; i < kNumMLPeeThresholds; ++i) {
		dMLPeeThresholds[i] = mlpEe[i];
	}
	for (int i = 0; i < kNumMLPpiThresholds; ++i) {
		dMLPpiThresholds[i] = mlpPi[i];
	}

	const float bdtEe[kNumBDTeeThresholds] = {-0.010f, 0.010f, 0.030f, 0.050f, 0.070f, 0.090f, 0.110f,
		0.130f, 0.150f, 0.170f, 0.190f, 0.210f, 0.230f, 0.250f};
	const float bdtPi[kNumBDTpiThresholds] = {-0.010f, -0.030f, -0.050f, -0.070f, -0.090f, -0.110f, -0.130f,
		-0.150f, -0.170f, -0.190f, -0.210f, -0.230f, -0.250f, -0.270f,
		-0.290f, -0.310f, -0.330f, -0.350f};
	for (int i = 0; i < kNumBDTeeThresholds; ++i) {
		dBDTeeThresholds[i] = bdtEe[i];
	}
	for (int i = 0; i < kNumBDTpiThresholds; ++i) {
		dBDTpiThresholds[i] = bdtPi[i];
	}

	const Int_t kNumBins = 200;
	const Double_t edges[kNumBins + 1] = {0.0000, 0.0002, 0.0006, 0.0013, 0.002, 0.003, 0.0042, 0.0055, 0.007, 0.009, 0.011,
		0.013, 0.015, 0.017, 0.019, 0.021, 0.023, 0.025, 0.027, 0.029, 0.031, 0.033, 0.035, 0.037, 0.039, 0.041, 0.043,
		0.045, 0.047, 0.049, 0.051, 0.053, 0.055, 0.057, 0.059, 0.061, 0.063, 0.065, 0.067, 0.069, 0.071, 0.073, 0.075,
		0.077, 0.079, 0.081, 0.083, 0.085, 0.087, 0.089, 0.091, 0.093, 0.095, 0.097, 0.099, 0.101, 0.103, 0.105, 0.107,
		0.109, 0.111, 0.113, 0.115, 0.117, 0.119, 0.121, 0.123, 0.125, 0.127, 0.129, 0.131, 0.133, 0.135, 0.137, 0.139,
		0.141, 0.143, 0.145, 0.147, 0.149, 0.151, 0.153, 0.155, 0.157, 0.159, 0.161, 0.163, 0.165, 0.167, 0.169, 0.171,
		0.173, 0.175, 0.177, 0.179, 0.181, 0.183, 0.185, 0.187, 0.189, 0.191, 0.193, 0.195, 0.197, 0.199, 0.201, 0.203,
		0.205, 0.207, 0.209, 0.211, 0.213, 0.215, 0.217, 0.219, 0.221, 0.223, 0.225, 0.227, 0.229, 0.231, 0.233, 0.235,
		0.237, 0.239, 0.241, 0.243, 0.245, 0.247, 0.249, 0.251, 0.253, 0.255, 0.257, 0.259, 0.261, 0.263, 0.265, 0.267,
		0.269, 0.271, 0.273, 0.275, 0.277, 0.279, 0.281, 0.283, 0.285, 0.287, 0.289, 0.291, 0.293, 0.295, 0.297, 0.299,
		0.301, 0.303, 0.305, 0.307, 0.309, 0.311, 0.313, 0.315, 0.317, 0.319, 0.321, 0.323, 0.325, 0.327, 0.329, 0.331,
		0.333, 0.335, 0.337, 0.339, 0.341, 0.343, 0.345, 0.347, 0.349, 0.351, 0.353, 0.355, 0.357, 0.359, 0.361, 0.363,
		0.365, 0.367, 0.369, 0.371, 0.373, 0.375, 0.377, 0.379, 0.381, 0.383, 0.385, 0.387, 0.389, 0.391};

	const bool enableThetaSystematics = true;
	const bool enableInvMassSystematics = true;
	const bool enableChiSqSystematics = true;

	if(enableThetaSystematics)
		thetaSysDir->cd();

	if(enableThetaSystematics) for (int i = 0; i < kNumMethods; i++) {
		for (int j = 0; j < kNumPid; j++) {
			thetaSysDir->cd();
			dHist_theta2_vs_theta1_2DCutOFF[i][j] = new TH2D(
				Form("theta2_vs_theta1_%s_%s_fidCuts_2DCutOFF", dMethods[i], dPidChoices[j]),
				";lab #theta_{1} (deg);lab #theta_{2} (deg)", 100, 0, 15, 100, 0, 15);
			dHist_theta1_noCuts[i][j] = new TH1D(
				Form("theta1_%s_%s_noCuts", dMethods[i], dPidChoices[j]),
				";lab #theta_{1} (deg)", 100, 0, 15);
			dHist_theta1_fidCuts_noTheta[i][j] = new TH1D(
				Form("theta1_%s_%s_fidCuts_noTheta", dMethods[i], dPidChoices[j]),
				"; #theta_{1} (deg)", 100, 0, 15);
			dHist_theta2_noCuts[i][j] = new TH1D(
				Form("theta2_%s_%s_noCuts", dMethods[i], dPidChoices[j]),
				";lab #theta_{2} (deg)", 100, 0, 15);
			dHist_theta2_fidCuts_noTheta[i][j] = new TH1D(
				Form("theta2_%s_%s_fidCuts_noTheta", dMethods[i], dPidChoices[j]),
				"; #theta_{2} (deg)", 100, 0, 15);

			phiJTDir->cd();
			dHist_JTphi_FID_2DCutOFF[i][j] = new TH1D(
				Form("JTphi_%s_%s_fidCuts_2DCutOFF", dMethods[i], dPidChoices[j]),
				";J_{T}.#phi (deg)", 90, -180, 180);

			q2Dir->cd();
			dHist_q2kin_FID_2DCutOFF[i][j] = new TH1D(
				Form("q2_%s_%s_fidCuts_2DCutOFF", dMethods[i], dPidChoices[j]),
				";q^{2} (GeV/c)^{2}", 30, 0.000000, 0.015);
			dHist_q2kin_varWidth_FID_2DCutOFF[i][j] = new TH1D(
				Form("q2_varWidth_%s_%s_fidCuts_2DCutOFF", dMethods[i], dPidChoices[j]),
				";q^{2} (GeV/c)^{2}", kNumBins, edges);
			if(dIsMC) {
				dHist_qvec2_res_vs_q2kin_FID_2DCutOFF[i][j] = new TH2D(
					Form("q2kinRes_vs_q2kin_%s_%s_fidCuts_2DCutOFF", dMethods[i], dPidChoices[j]),
					";#vec{q}^{2}_{kin} (GeV/c)^{2};#vec{q}^{2}_{kin} - #vec{q}^{2}_{thrown} (GeV/c)^{2}",
					200, 0.000000, 0.015, 200, -0.0015, 0.0015);
			}

			for (int k = 0; k < kNumAngles; k++) {
				phiJTDir->cd();
				dHist_JTphi_angles[i][j][k] = new TH1D(
					Form("JTphi_%s_%s_thetaCut_%.2f", dMethods[i], dPidChoices[j], dMinAngles[k]),
					";J_{T}.#phi (deg)", 90, -180, 180);

				q2Dir->cd();
				dHist_q2kin_angles[i][j][k] = new TH1D(
					Form("q2_%s_%s_thetaCut_%.2f", dMethods[i], dPidChoices[j], dMinAngles[k]),
					";q^{2} (GeV/c)^{2}", 30, 0.000000, 0.015);
				dHist_q2kin_varWidth_angles[i][j][k] = new TH1D(
					Form("q2_varWidth_%s_%s_thetaCut_%.2f", dMethods[i], dPidChoices[j], dMinAngles[k]),
					";q^{2} (GeV/c)^{2}", kNumBins, edges);
				if(dIsMC) {
					dHist_qvec2_res_vs_q2kin_angles[i][j][k] = new TH2D(
						Form("q2kinRes_vs_q2kin_%s_%s_thetaCut_%.2f", dMethods[i], dPidChoices[j], dMinAngles[k]),
						";#vec{q}^{2}_{kin} (GeV/c)^{2};#vec{q}^{2}_{kin} - #vec{q}^{2}_{thrown} (GeV/c)^{2}",
						200, 0.000000, 0.015, 200, -0.0015, 0.0015);
				}
			}
		}
	}

	for (int i = 0; i < kNumMethods; i++) {
		for (int j = 0; j < kNumPid; j++) {
			for (int k = 0; k < kNumBeamRegions; k++) {
				phiJTDir->cd();
				dHist_JTphi_BeamRegions[i][j][k] = new TH1D(
					Form("JTphi_%s_%s_beamRegion_%.3f_%.3f", dMethods[i], dPidChoices[j],
						dBeamEnergyRegions[k][0], dBeamEnergyRegions[k][1]),
					";J_{T}.#phi (deg)", 90, -180, 180);
				q2Dir->cd();
				dHist_q2kin_BeamRegions[i][j][k] = new TH1D(
					Form("q2_%s_%s_beamRegion_%.3f_%.3f", dMethods[i], dPidChoices[j],
						dBeamEnergyRegions[k][0], dBeamEnergyRegions[k][1]),
					";q^{2} (GeV/c)^{2}", 30, 0.000000, 0.015);
				dHist_q2kin_varWidth_BeamRegions[i][j][k] = new TH1D(
					Form("q2_varWidth_%s_%s_beamRegion_%.3f_%.3f", dMethods[i], dPidChoices[j],
						dBeamEnergyRegions[k][0], dBeamEnergyRegions[k][1]),
					";q^{2} (GeV/c)^{2}", kNumBins, edges);
				if(dIsMC) {
					dHist_qvec2_res_vs_q2kin_BeamRegions[i][j][k] = new TH2D(
						Form("q2kinRes_vs_q2kin_%s_%s_beamRegion_%.3f_%.3f", dMethods[i], dPidChoices[j],
							dBeamEnergyRegions[k][0], dBeamEnergyRegions[k][1]),
						";#vec{q}^{2}_{kin} (GeV/c)^{2};#vec{q}^{2}_{kin} - #vec{q}^{2}_{thrown} (GeV/c)^{2}",
						200, 0.000000, 0.015, 200, -0.0015, 0.0015);
				}
			}
		}
	}

	if(enableInvMassSystematics)
		invmassSysDir->cd();
	if(enableInvMassSystematics) for (int i = 0; i < kNumMethods; i++) {
		for (int j = 0; j < kNumPid; j++) {
			invmassSysDir->cd();
			dHist_Wepem_kin_noCuts[i][j] = new TH1D(
				Form("Wepem_kin_%s_%s_noCuts", dMethods[i], dPidChoices[j]),
				";M_{e^{+}e^{-}} Kin (GeV/c^{2})", 238, 0.0, 1.4);
			dHist_Wepem_kin_fidCuts_noWcut[i][j] = new TH1D(
				Form("Wepem_kin_%s_%s_fidCuts_noWcut", dMethods[i], dPidChoices[j]),
				";M_{e^{+}e^{-}} Kin (GeV/c^{2})", 238, 0.0, 1.4);
			dHist_Wepem_Measured_noCuts[i][j] = new TH1D(
				Form("Wepem_Measured_%s_%s_noCuts", dMethods[i], dPidChoices[j]),
				";M_{e^{+}e^{-}} Gen (GeV/c^{2})", 238, 0.0, 1.4);
			dHist_Wepem_Measured_fidCuts_noWcut[i][j] = new TH1D(
				Form("Wepem_Measured_%s_%s_fidCuts_noWcut", dMethods[i], dPidChoices[j]),
				";M_{e^{+}e^{-}} Gen (GeV/c^{2})", 238, 0.0, 1.4);

			for (int k = 0; k < kNumMassRegions; k++) {
				invmassSysDir->cd();
				dHist_Wkin_MassRegions[i][j][k] = new TH1D(
					Form("W_%s_%s_massRegion_%.3f_%.3f", dMethods[i], dPidChoices[j],
						dMassRegions[k][0], dMassRegions[k][1]),
					";M_{e^{+}e^{-}} Gen (GeV/c^{2})", 238, 0.0, 1.4);

				phiJTDir->cd();
				dHist_JTphi_MassRegions[i][j][k] = new TH1D(
					Form("JTphi_%s_%s_massRegion_%.3f_%.3f", dMethods[i], dPidChoices[j],
						dMassRegions[k][0], dMassRegions[k][1]),
					";J_{T}.#phi (deg)", 90, -180, 180);

				q2Dir->cd();
				dHist_q2kin_MassRegions[i][j][k] = new TH1D(
					Form("q2_%s_%s_massRegion_%.3f_%.3f", dMethods[i], dPidChoices[j],
						dMassRegions[k][0], dMassRegions[k][1]),
					";q^{2} (GeV/c)^{2}", 30, 0.000000, 0.015);
				dHist_q2kin_varWidth_MassRegions[i][j][k] = new TH1D(
					Form("q2_varWidth_%s_%s_massRegion_%.3f_%.3f", dMethods[i], dPidChoices[j],
						dMassRegions[k][0], dMassRegions[k][1]),
					";q^{2} (GeV/c)^{2}", kNumBins, edges);
				if(dIsMC) {
					dHist_qvec2_res_vs_q2kin_MassRegions[i][j][k] = new TH2D(
						Form("q2kinRes_vs_q2kin_%s_%s_massRegion_%.3f_%.3f", dMethods[i], dPidChoices[j],
							dMassRegions[k][0], dMassRegions[k][1]),
						";#vec{q}^{2}_{kin} (GeV/c)^{2};#vec{q}^{2}_{kin} - #vec{q}^{2}_{thrown} (GeV/c)^{2}",
						200, 0.000000, 0.015, 200, -0.0015, 0.0015);
				}
			}

			for (int k = 0; k < kNumWminCuts; k++) {
				phiJTDir->cd();
				dHist_JTphi_WminCuts[i][j][k] = new TH1D(
					Form("JTphi_%s_%s_Wmin_%.2f", dMethods[i], dPidChoices[j], dMinInvMass[k]),
					";J_{T}.#phi (deg)", 90, -180, 180);
				q2Dir->cd();
				dHist_q2kin_WminCuts[i][j][k] = new TH1D(
					Form("q2_%s_%s_Wmin_%.2f", dMethods[i], dPidChoices[j], dMinInvMass[k]),
					";q^{2} (GeV/c)^{2}", 30, 0.000000, 0.015);
				dHist_q2kin_varWidth_WminCuts[i][j][k] = new TH1D(
					Form("q2_varWidth_%s_%s_Wmin_%.2f", dMethods[i], dPidChoices[j], dMinInvMass[k]),
					";q^{2} (GeV/c)^{2}", kNumBins, edges);
				if(dIsMC) {
					dHist_qvec2_res_vs_q2kin_WminCuts[i][j][k] = new TH2D(
						Form("q2kinRes_vs_q2kin_%s_%s_Wmin_%.2f", dMethods[i], dPidChoices[j], dMinInvMass[k]),
						";#vec{q}^{2}_{kin} (GeV/c)^{2};#vec{q}^{2}_{kin} - #vec{q}^{2}_{thrown} (GeV/c)^{2}",
						200, 0.000000, 0.015, 200, -0.0015, 0.0015);
				}
			}
		}
	}

	if(enableChiSqSystematics)
		chisqSysDir->cd();
	if(enableChiSqSystematics) for (int i = 0; i < kNumMethods; i++) {
		for (int j = 0; j < kNumPid; j++) {
			chisqSysDir->cd();
			dHist_KinFitChiSq_noCuts[i][j] = new TH1D(
				Form("KinFitChiSq_%s_%s_noCuts", dMethods[i], dPidChoices[j]),
				";Kinematic Fit #chi^{2}/NDF", 250, 0., 25.);
			dHist_KinFitChiSq_fidCuts[i][j] = new TH1D(
				Form("KinFitChiSq_%s_%s_fidCuts_noChiSq", dMethods[i], dPidChoices[j]),
				";Kinematic Fit #chi^{2}/NDF", 250, 0., 25.);

			for (int k = 0; k < kNumChiSqCuts; k++) {
				phiJTDir->cd();
				dHist_JTphi_MaxChiSq[i][j][k] = new TH1D(
					Form("JTphi_%s_%s_maxChiSq_%.2f", dMethods[i], dPidChoices[j], dMaxChiSq[k]),
					";J_{T}.#phi (deg)", 90, -180, 180);

				q2Dir->cd();
				dHist_q2kin_MaxChiSq[i][j][k] = new TH1D(
					Form("q2_%s_%s_maxChiSq_%.2f", dMethods[i], dPidChoices[j], dMaxChiSq[k]),
					";q^{2} (GeV/c)^{2}", 30, 0.000000, 0.015);
				dHist_q2kin_varWidth_MaxChiSq[i][j][k] = new TH1D(
					Form("q2_varWidth_%s_%s_maxChiSq_%.2f", dMethods[i], dPidChoices[j], dMaxChiSq[k]),
					";q^{2} (GeV/c)^{2}", kNumBins, edges);
				if(dIsMC) {
					dHist_qvec2_res_vs_q2kin_MaxChiSq[i][j][k] = new TH2D(
						Form("q2kinRes_vs_q2kin_%s_%s_maxChiSq_%.2f", dMethods[i], dPidChoices[j], dMaxChiSq[k]),
						";#vec{q}^{2}_{kin} (GeV/c)^{2};#vec{q}^{2}_{kin} - #vec{q}^{2}_{thrown} (GeV/c)^{2}",
						200, 0.000000, 0.015, 200, -0.0015, 0.0015);
				}
			}
		}
	}

	for (int i = 0; i < kNumMLPeeThresholds; i++) {
		phiJTDir->cd();
		dHist_JTphi_MLP_ee[i] = new TH1D(
			Form("JTphi_MLP_ee_MVAcut_%.2f", dMLPeeThresholds[i]),
			";J_{T}.#phi (deg)", 90, -180, 180);
		q2Dir->cd();
		dHist_q2kin_MLP_ee[i] = new TH1D(
			Form("q2_MLP_ee_MVAcut_%.2f", dMLPeeThresholds[i]),
			";q^{2} (GeV/c)^{2}", 30, 0.000000, 0.015);
		dHist_q2kin_varWidth_MLP_ee[i] = new TH1D(
			Form("q2_varWidth_MLP_ee_MVAcut_%.2f", dMLPeeThresholds[i]),
			";q^{2} (GeV/c)^{2}", kNumBins, edges);
		if(dIsMC) {
			dHist_qvec2_res_vs_q2kin_MLP_ee[i] = new TH2D(
				Form("q2kinRes_vs_q2kin_MLP_ee_MVAcut_%.2f", dMLPeeThresholds[i]),
				";#vec{q}^{2}_{kin} (GeV/c)^{2};#vec{q}^{2}_{kin} - #vec{q}^{2}_{thrown} (GeV/c)^{2}",
				200, 0.000000, 0.015, 200, -0.0015, 0.0015);
		}
	}
	for (int i = 0; i < kNumMLPpiThresholds; i++) {
		phiJTDir->cd();
		dHist_JTphi_MLP_pi[i] = new TH1D(
			Form("JTphi_MLP_pi_MVAcut_%.2f", dMLPpiThresholds[i]),
			";J_{T}.#phi (deg)", 90, -180, 180);
		q2Dir->cd();
		dHist_q2kin_MLP_pi[i] = new TH1D(
			Form("q2_MLP_pi_MVAcut_%.2f", dMLPpiThresholds[i]),
			";q^{2} (GeV/c)^{2}", 30, 0.000000, 0.015);
		dHist_q2kin_varWidth_MLP_pi[i] = new TH1D(
			Form("q2_varWidth_MLP_pi_MVAcut_%.2f", dMLPpiThresholds[i]),
			";q^{2} (GeV/c)^{2}", kNumBins, edges);
		if(dIsMC) {
			dHist_qvec2_res_vs_q2kin_MLP_pi[i] = new TH2D(
				Form("q2kinRes_vs_q2kin_MLP_pi_MVAcut_%.2f", dMLPpiThresholds[i]),
				";#vec{q}^{2}_{kin} (GeV/c)^{2};#vec{q}^{2}_{kin} - #vec{q}^{2}_{thrown} (GeV/c)^{2}",
				200, 0.000000, 0.015, 200, -0.0015, 0.0015);
		}
	}

	for (int i = 0; i < kNumBDTeeThresholds; i++) {
		phiJTDir->cd();
		dHist_JTphi_BDT_ee[i] = new TH1D(
			Form("JTphi_BDT_ee_MVAcut_%.3f", dBDTeeThresholds[i]),
			";J_{T}.#phi (deg)", 90, -180, 180);
		q2Dir->cd();
		dHist_q2kin_BDT_ee[i] = new TH1D(
			Form("q2_BDT_ee_MVAcut_%.3f", dBDTeeThresholds[i]),
			";q^{2} (GeV/c)^{2}", 30, 0.000000, 0.015);
		dHist_q2kin_varWidth_BDT_ee[i] = new TH1D(
			Form("q2_varWidth_BDT_ee_MVAcut_%.3f", dBDTeeThresholds[i]),
			";q^{2} (GeV/c)^{2}", kNumBins, edges);
		if(dIsMC) {
			dHist_qvec2_res_vs_q2kin_BDT_ee[i] = new TH2D(
				Form("q2kinRes_vs_q2kin_BDT_ee_MVAcut_%.3f", dBDTeeThresholds[i]),
				";#vec{q}^{2}_{kin} (GeV/c)^{2};#vec{q}^{2}_{kin} - #vec{q}^{2}_{thrown} (GeV/c)^{2}",
				200, 0.000000, 0.015, 200, -0.0015, 0.0015);
		}
	}
	for (int i = 0; i < kNumBDTpiThresholds; i++) {
		phiJTDir->cd();
		dHist_JTphi_BDT_pi[i] = new TH1D(
			Form("JTphi_BDT_pi_MVAcut_%.3f", dBDTpiThresholds[i]),
			";J_{T}.#phi (deg)", 90, -180, 180);
		q2Dir->cd();
		dHist_q2kin_BDT_pi[i] = new TH1D(
			Form("q2_BDT_pi_MVAcut_%.3f", dBDTpiThresholds[i]),
			";q^{2} (GeV/c)^{2}", 30, 0.000000, 0.015);
		dHist_q2kin_varWidth_BDT_pi[i] = new TH1D(
			Form("q2_varWidth_BDT_pi_MVAcut_%.3f", dBDTpiThresholds[i]),
			";q^{2} (GeV/c)^{2}", kNumBins, edges);
		if(dIsMC) {
			dHist_qvec2_res_vs_q2kin_BDT_pi[i] = new TH2D(
				Form("q2kinRes_vs_q2kin_BDT_pi_MVAcut_%.3f", dBDTpiThresholds[i]),
				";#vec{q}^{2}_{kin} (GeV/c)^{2};#vec{q}^{2}_{kin} - #vec{q}^{2}_{thrown} (GeV/c)^{2}",
				200, 0.000000, 0.015, 200, -0.0015, 0.0015);
		}
	}

	mvaResponseDir->cd();
	dHist_MLP_responsePIP = new TH1D("MLP_responsePIP", ";MLP_responsePIP", 100, 0.0, 1.05);
	dHist_MLP_responsePIM = new TH1D("MLP_responsePIM", ";MLP_responsePIM", 100, 0.0, 1.05);
	dHist_MLP_response_PIP_vs_PIM = new TH2D("MLP_PIP_vs_PIM", ";MLP PIM ;MLP PIP", 100, 0, 1, 100, 0, 1);

	dHist_BDT_responsePIP = new TH1D("BDT_responsePIP", ";BDT_responsePIP", 100, -0.5, 0.5);
	dHist_BDT_responsePIM = new TH1D("BDT_responsePIM", ";BDT_responsePIM", 100, -0.5, 0.5);
	dHist_BDT_response_PIP_vs_PIM = new TH2D("BDT_PIP_vs_PIM", ";BDT PIM ;BDT PIP", 100, -0.5, 0.5, 100, -0.5, 0.5);

	mainDir->cd();
}

bool DSelector_2eMissingProton_Systematics::ComboPassesCuts(const ComboCutInputs& inputs, const CutConditions& conditions)
{
	// std::cout << "DEBUG COMBOCUTS 1: Entered ComboPassesCuts, applyMaxBeamEcut=" << conditions.applyMaxBeamEcut << ", maxBeamE=" << conditions.maxBeamE << std::endl;
	std::cout.flush();
	
	bool passedAllCuts = true;  // Use single return point
	
	if (passedAllCuts && conditions.applyPreselectionEoPCut) {
		if (inputs.FCAL_EoverPmeas_plus < conditions.preselectionMinEoverP || inputs.FCAL_EoverPmeas_minus < conditions.preselectionMinEoverP) {
			// std::cout << "DEBUG COMBOCUTS 1.5: Failed preselection" << std::endl;
			passedAllCuts = false;
		}
	}

	if (passedAllCuts) {
		// std::cout << "DEBUG COMBOCUTS 2: Passed preselection" << std::endl;
		std::cout.flush();
	}
	
	if (passedAllCuts && conditions.applyNoUnusedTracksCut) {
		if (inputs.numUnusedTracks > conditions.maxNumUnusedTracks) {
			// std::cout << "DEBUG COMBOCUTS 2.5: Failed unused tracks" << std::endl;
			passedAllCuts = false;
		}
	}

	if (passedAllCuts) {
		// std::cout << "DEBUG COMBOCUTS 3: Passed unused tracks" << std::endl;
		std::cout.flush();
	}
	
	if (passedAllCuts && conditions.applyNoUnusedShowersCut) {
		if (inputs.energyUnusedShowers > conditions.maxUnusedShowersEnergy) {
			// std::cout << "DEBUG COMBOCUTS 3.5: Failed unused showers" << std::endl;
			passedAllCuts = false;
		}
	}

	if (passedAllCuts) {
		// std::cout << "DEBUG COMBOCUTS 4: Passed unused showers" << std::endl;
		std::cout.flush();
	}

	if (passedAllCuts && conditions.applyFCALEnergyNonZeroCut) {
		// std::cout << "DEBUG COMBOCUTS 4.1: Checking FCAL energy non-zero" << std::endl;
		std::cout.flush();
		if (inputs.FCAL_Energy_plus <= 0.0 || inputs.FCAL_Energy_minus <= 0.0) {
			// std::cout << "DEBUG COMBOCUTS 4.1.5: Failed FCAL energy" << std::endl;
			passedAllCuts = false;
		}
	}

	if (passedAllCuts && conditions.applyTOFdEdxNonZeroCut) {
		// std::cout << "DEBUG COMBOCUTS 4.2: Checking TOF dEdx non-zero" << std::endl;
		std::cout.flush();
		if (inputs.TOF_dEdx_plus <= 0.0 || inputs.TOF_dEdx_minus <= 0.0) {
			// std::cout << "DEBUG COMBOCUTS 4.2.5: Failed TOF dEdx" << std::endl;
			passedAllCuts = false;
		}
	}

	if (passedAllCuts && conditions.applyMinBeamEcut) {
		// std::cout << "DEBUG COMBOCUTS 4.3: Checking minBeamE" << std::endl;
		std::cout.flush();
		if (inputs.beamE < conditions.minBeamE) {
			// std::cout << "DEBUG COMBOCUTS 4.3.5: Failed minBeamE" << std::endl;
			passedAllCuts = false;
		}
	}

	if (passedAllCuts && conditions.applyMaxBeamEcut) {
		// std::cout << "DEBUG COMBOCUTS 4.4: Checking maxBeamE, beamE=" << inputs.beamE << ", maxBeamE=" << conditions.maxBeamE << std::endl;
		std::cout.flush();
		if (inputs.beamE > conditions.maxBeamE) {
			// std::cout << "DEBUG COMBOCUTS 4.4.A: Beam energy check failed, setting passedAllCuts=false" << std::endl;
			std::cout.flush();
			passedAllCuts = false;
		} else {
			// std::cout << "DEBUG COMBOCUTS 4.4.0: Passed maxBeamE check" << std::endl;
			std::cout.flush();
		}
	}

	// std::cout << "DEBUG COMBOCUTS 4.4.1: About to access applyMinThetaCut" << std::endl;
	std::cout.flush();
	bool minThetaFlag = conditions.applyMinThetaCut;
	// std::cout << "DEBUG COMBOCUTS 4.4.2: applyMinThetaCut = " << minThetaFlag << std::endl;
	std::cout.flush();
	
	if (passedAllCuts && minThetaFlag) {
		// cout << "DEBUG COMBOCUTS 4.5: Checking minTheta" << endl;
		if (inputs.theta1 < conditions.minTheta || inputs.theta2 < conditions.minTheta)
			passedAllCuts = false;
	}

	if (passedAllCuts && conditions.applyMaxThetaCut) {
		// cout << "DEBUG COMBOCUTS 4.6: Checking maxTheta" << endl;
		if (inputs.theta1 > conditions.maxTheta || inputs.theta2 > conditions.maxTheta)
			passedAllCuts = false;
	}

	if (passedAllCuts && conditions.apply2DThetaCut) {
		// cout << "DEBUG COMBOCUTS 4.7: About to compute 2D theta cut, theta1=" << inputs.theta1 << ", theta2=" << inputs.theta2 << endl;
		if (14 * inputs.theta1 / (TMath::Power(inputs.theta1, 2.5) + 0.2) < inputs.theta2)
			passedAllCuts = false;
		// else
			// cout << "DEBUG COMBOCUTS 4.8: Passed 2D theta cut" << endl;
	}

	// cout << "DEBUG COMBOCUTS 4.9: About to check momentum range" << endl;
	if (passedAllCuts && conditions.applyMomentumRangeCut) {
		// cout << "DEBUG COMBOCUTS 4.9.1: Checking momentum range" << endl;
		if (inputs.pMagPlus < conditions.minPmag || inputs.pMagPlus > conditions.maxPmag)
			passedAllCuts = false;
		if (passedAllCuts && (inputs.pMagMinus < conditions.minPmag || inputs.pMagMinus > conditions.maxPmag))
			passedAllCuts = false;
	}

	// cout << "DEBUG COMBOCUTS 4.10: About to check Wkin cuts" << endl;
	if (passedAllCuts && conditions.applyMinWkinCut) {
		if (inputs.Wkin < conditions.minWkin)
			passedAllCuts = false;
	}

	if (passedAllCuts && conditions.applyMaxWkinCut) {
		if (inputs.Wkin > conditions.maxWkin)
			passedAllCuts = false;
	}

	// cout << "DEBUG COMBOCUTS 4.11: About to check vertex Z" << endl;
	if (passedAllCuts && conditions.applyVertexZCut) {
		if (inputs.vertexZ < conditions.minVertexZ || inputs.vertexZ > conditions.maxVertexZ)
			passedAllCuts = false;
	}

	// cout << "DEBUG COMBOCUTS 4.12: About to check chi-squared" << endl;
	if (passedAllCuts && conditions.applyChiSqCut) {
		if (inputs.chisqndf > conditions.maxKFChiSq)
			passedAllCuts = false;
	}

	// cout << "DEBUG COMBOCUTS 4.13: About to check kinfit CL" << endl;
	if (passedAllCuts && conditions.applyKinFitCLCut) {
		if (inputs.kinFitCL <= conditions.minKinFitCL)
			passedAllCuts = false;
	}

	if (passedAllCuts && conditions.applyBestChiSqComboCut) {
		if (!inputs.ThisComboIsBestChiSq)
			passedAllCuts = false;
	}

	if (passedAllCuts && conditions.applyMinEoverP1Cut) {
		if (inputs.FCAL_EoverPkin_plus < conditions.minEoverP1)
			passedAllCuts = false;
	}

	if (passedAllCuts && conditions.applyMaxEoverP1Cut) {
		if (inputs.FCAL_EoverPkin_plus > conditions.maxEoverP1)
			passedAllCuts = false;
	}

	if (passedAllCuts && conditions.applyMinEoverP2Cut) {
		if (inputs.FCAL_EoverPkin_minus < conditions.minEoverP2)
			passedAllCuts = false;
	}

	if (passedAllCuts && conditions.applyMaxEoverP2Cut) {
		if (inputs.FCAL_EoverPkin_minus > conditions.maxEoverP2)
			passedAllCuts = false;
	}

	if (passedAllCuts && conditions.applyMinFCALElasticityCut) {
		if (inputs.FCAL_Elasticity < conditions.minFCALElasticity)
			passedAllCuts = false;
	}

	if (passedAllCuts && conditions.applyMaxFCALElasticityCut) {
		if (inputs.FCAL_Elasticity > conditions.maxFCALElasticity)
			passedAllCuts = false;
	}

	if (passedAllCuts && conditions.applyMaxFCALDOCACut) {
		if (inputs.TrackFCAL_DOCA_plus > conditions.maxFCALDOCA || inputs.TrackFCAL_DOCA_minus > conditions.maxFCALDOCA)
			passedAllCuts = false;
	}

	if (passedAllCuts && conditions.applyMaxTOFdEdxCut) {
		if (inputs.TOF_dEdx_plus > conditions.maxTOFdEdx || inputs.TOF_dEdx_minus > conditions.maxTOFdEdx)
			passedAllCuts = false;
	}

	// cout << "DEBUG COMBOCUTS 5: About to check applyMVACuts=" << conditions.applyMVACuts << endl;
	if (passedAllCuts && conditions.applyMVACuts) {
		// cout << "DEBUG COMBOCUTS 6: modelChoice=" << conditions.modelChoice << ", particleChoice=" << conditions.particleChoice << endl;
		// cout << "DEBUG COMBOCUTS 7: MLP1=" << inputs.MLP1 << ", MLP2=" << inputs.MLP2 << ", BDT1=" << inputs.BDT1 << ", BDT2=" << inputs.BDT2 << endl;
		
		if (conditions.modelChoice == "BDT") {
			// cout << "DEBUG COMBOCUTS 8: Checking BDT cuts" << endl;
			if (conditions.particleChoice == "ee") {
				if (inputs.BDT1 < conditions.BDT_ee || inputs.BDT2 < conditions.BDT_ee)
					passedAllCuts = false;
			} else if (conditions.particleChoice == "pi") {
				if (inputs.BDT1 > conditions.BDT_pi || inputs.BDT2 > conditions.BDT_pi)
					passedAllCuts = false;
			}
		} else if (conditions.modelChoice == "MLP") {
			// cout << "DEBUG COMBOCUTS 9: Checking MLP cuts, MLP_ee=" << conditions.MLP_ee << ", MLP_pi=" << conditions.MLP_pi << endl;
			if (conditions.particleChoice == "ee") {
				// cout << "DEBUG COMBOCUTS 10: ee particle choice" << endl;
				if (inputs.MLP1 < conditions.MLP_ee || inputs.MLP2 < conditions.MLP_ee)
					passedAllCuts = false;
			} else if (conditions.particleChoice == "pi") {
				// cout << "DEBUG COMBOCUTS 11: pi particle choice" << endl;
				if (inputs.MLP1 > conditions.MLP_pi || inputs.MLP2 > conditions.MLP_pi)
					passedAllCuts = false;
			}
		}
	}

	if (passedAllCuts && conditions.applyMLPResponseCut) {
		if (inputs.MLP1 < conditions.minMLPResponse || inputs.MLP2 < conditions.minMLPResponse)
			passedAllCuts = false;
	}

	// cout << "DEBUG COMBOCUTS 12: Returning " << (passedAllCuts ? "TRUE" : "FALSE") << " from ComboPassesCuts" << endl;
	std::cout.flush();
	return passedAllCuts;
}

void DSelector_2eMissingProton_Systematics::FillSystematics(const ComboCutInputs& inputs, Double_t q2kin, Double_t JTphi, Double_t Wmeas, Double_t q2kinRes, bool includeQ2ResSystematics, int runPeriodIndex, int polarizationIndex)
{
	// cout << "DEBUG FILLSYST 1: Entered FillSystematics function" << endl;
	const bool fillQ2ResSystematics = includeQ2ResSystematics && dIsMC;
	const bool enableThetaSystematics = true;
	const bool enableInvMassSystematics = true;
	const bool enableChiSqSystematics = true;
	const bool fillJTphiByRunPol = (runPeriodIndex >= 0 && runPeriodIndex < kNumRunPeriods && polarizationIndex >= 0 && polarizationIndex < nPolarizations);

	const CutConditions activeFiducialConditions = BuildActiveFiducialConditions();

	auto setJTphiBeamWindow = [&](CutConditions& conditions) {
		if(dIsCPPRunPeriod) {
			conditions.minBeamE = 5.35;
			conditions.maxBeamE = 5.75;
		} else {
			conditions.minBeamE = 8.2;
			conditions.maxBeamE = 8.8;
		}
	};

	auto setQ2BeamWindow = [&](CutConditions& conditions) {
		if(dIsCPPRunPeriod) {
			conditions.minBeamE = 4.0;
			conditions.maxBeamE = 7.6;
		} else {
			conditions.minBeamE = 7.0;
			conditions.maxBeamE = 11.4;
		}
	};

	CutConditions fidCutsOff = activeFiducialConditions;
	setQ2BeamWindow(fidCutsOff);
	fidCutsOff.applyMinThetaCut = false;
	fidCutsOff.applyMaxThetaCut = false;
	fidCutsOff.apply2DThetaCut = false;
	fidCutsOff.applyMinWkinCut = false;
	fidCutsOff.applyMaxWkinCut = false;
	fidCutsOff.applyChiSqCut = false;

	// cout << "DEBUG FILLSYST 2: About to start first loop (noCuts)" << endl;
	if(enableThetaSystematics || enableInvMassSystematics || enableChiSqSystematics) for (int i = 0; i < kNumMethods; i++) {
		fidCutsOff.modelChoice = dMethods[i];
		for (int j = 0; j < kNumPid; j++) {
			fidCutsOff.particleChoice = dPidChoices[j];
			if (ComboPassesCuts(inputs, fidCutsOff)) {
				dHist_theta1_noCuts[i][j]->Fill(inputs.theta1);
				dHist_theta2_noCuts[i][j]->Fill(inputs.theta2);
				dHist_Wepem_kin_noCuts[i][j]->Fill(inputs.Wkin);
				dHist_Wepem_Measured_noCuts[i][j]->Fill(Wmeas);
				dHist_KinFitChiSq_noCuts[i][j]->Fill(inputs.chisqndf);
			}
		}
	}
	// cout << "DEBUG FILLSYST 3: Completed first loop" << endl;

	CutConditions fidButThetaOff = activeFiducialConditions;
	setQ2BeamWindow(fidButThetaOff);
	fidButThetaOff.applyMinThetaCut = false;
	fidButThetaOff.applyMaxThetaCut = false;
	fidButThetaOff.apply2DThetaCut = false;
	if(enableThetaSystematics) for (int i = 0; i < kNumMethods; i++) {
		fidButThetaOff.modelChoice = dMethods[i];
		for (int j = 0; j < kNumPid; j++) {
			fidButThetaOff.particleChoice = dPidChoices[j];
			if (ComboPassesCuts(inputs, fidButThetaOff)) {
				dHist_theta1_fidCuts_noTheta[i][j]->Fill(inputs.theta1);
				dHist_theta2_fidCuts_noTheta[i][j]->Fill(inputs.theta2);
			}
		}
	}

	CutConditions fidButWOff = activeFiducialConditions;
	setQ2BeamWindow(fidButWOff);
	fidButWOff.applyMinWkinCut = false;
	fidButWOff.applyMaxWkinCut = false;
	if(enableInvMassSystematics) for (int i = 0; i < kNumMethods; i++) {
		fidButWOff.modelChoice = dMethods[i];
		for (int j = 0; j < kNumPid; j++) {
			fidButWOff.particleChoice = dPidChoices[j];
			if (ComboPassesCuts(inputs, fidButWOff)) {
				dHist_Wepem_kin_fidCuts_noWcut[i][j]->Fill(inputs.Wkin);
				dHist_Wepem_Measured_fidCuts_noWcut[i][j]->Fill(Wmeas);
			}
		}
	}

	CutConditions thetaSystConditions = activeFiducialConditions;
	setQ2BeamWindow(thetaSystConditions);
	thetaSystConditions.apply2DThetaCut = false;
	if(enableThetaSystematics) for (int i = 0; i < kNumMethods; i++) {
		thetaSystConditions.modelChoice = dMethods[i];
		for (int j = 0; j < kNumPid; j++) {
			thetaSystConditions.particleChoice = dPidChoices[j];
			if (ComboPassesCuts(inputs, thetaSystConditions)) {
				dHist_theta2_vs_theta1_2DCutOFF[i][j]->Fill(inputs.theta1, inputs.theta2);
				dHist_JTphi_FID_2DCutOFF[i][j]->Fill(JTphi);
				dHist_q2kin_FID_2DCutOFF[i][j]->Fill(q2kin);
				dHist_q2kin_varWidth_FID_2DCutOFF[i][j]->Fill(q2kin);
				if(fillQ2ResSystematics)
					dHist_qvec2_res_vs_q2kin_FID_2DCutOFF[i][j]->Fill(q2kin, q2kinRes);
			}
		}
	}

	CutConditions conditions = activeFiducialConditions;
	for (int i = 0; i < kNumMethods; i++) {
		conditions.modelChoice = dMethods[i];
		for (int j = 0; j < kNumPid; j++) {
			conditions.particleChoice = dPidChoices[j];
			conditions.minTheta = activeFiducialConditions.minTheta;
			conditions.minWkin = activeFiducialConditions.minWkin;
			conditions.maxWkin = activeFiducialConditions.maxWkin;
			conditions.maxKFChiSq = activeFiducialConditions.maxKFChiSq;
			conditions.applyChiSqCut = activeFiducialConditions.applyChiSqCut;

			if(enableThetaSystematics) for (int k = 0; k < kNumAngles; k++) {
				conditions.minTheta = dMinAngles[k];
				setJTphiBeamWindow(conditions);
				if (ComboPassesCuts(inputs, conditions))
					dHist_JTphi_angles[i][j][k]->Fill(JTphi);
				setQ2BeamWindow(conditions);
				if (ComboPassesCuts(inputs, conditions)) {
					dHist_q2kin_angles[i][j][k]->Fill(q2kin);
					dHist_q2kin_varWidth_angles[i][j][k]->Fill(q2kin);
					if(fillQ2ResSystematics)
						dHist_qvec2_res_vs_q2kin_angles[i][j][k]->Fill(q2kin, q2kinRes);
				}
			}
			conditions.minTheta = activeFiducialConditions.minTheta;

			if(enableInvMassSystematics) for (int k = 0; k < kNumWminCuts; k++) {
				conditions.minWkin = dMinInvMass[k];
				setJTphiBeamWindow(conditions);
				if (ComboPassesCuts(inputs, conditions))
					dHist_JTphi_WminCuts[i][j][k]->Fill(JTphi);
				setQ2BeamWindow(conditions);
				if (ComboPassesCuts(inputs, conditions)) {
					dHist_q2kin_WminCuts[i][j][k]->Fill(q2kin);
					dHist_q2kin_varWidth_WminCuts[i][j][k]->Fill(q2kin);
					if(fillQ2ResSystematics)
						dHist_qvec2_res_vs_q2kin_WminCuts[i][j][k]->Fill(q2kin, q2kinRes);
				}
			}
			conditions.minWkin = activeFiducialConditions.minWkin;

			for (int k = 0; k < kNumBeamRegions; k++) {
				conditions.minBeamE = dBeamEnergyRegions[k][0];
				conditions.maxBeamE = dBeamEnergyRegions[k][1];
				if (ComboPassesCuts(inputs, conditions)) {
					dHist_JTphi_BeamRegions[i][j][k]->Fill(JTphi);
					dHist_q2kin_BeamRegions[i][j][k]->Fill(q2kin);
					dHist_q2kin_varWidth_BeamRegions[i][j][k]->Fill(q2kin);
					if(fillQ2ResSystematics)
						dHist_qvec2_res_vs_q2kin_BeamRegions[i][j][k]->Fill(q2kin, q2kinRes);
				}
			}
			setQ2BeamWindow(conditions);

			if(enableInvMassSystematics) for (int k = 0; k < kNumMassRegions; k++) {
				conditions.minWkin = dMassRegions[k][0];
				conditions.maxWkin = dMassRegions[k][1];
				setJTphiBeamWindow(conditions);
				if (ComboPassesCuts(inputs, conditions))
					dHist_JTphi_MassRegions[i][j][k]->Fill(JTphi);
				setQ2BeamWindow(conditions);
				if (ComboPassesCuts(inputs, conditions)) {
					dHist_Wkin_MassRegions[i][j][k]->Fill(inputs.Wkin);
					dHist_q2kin_MassRegions[i][j][k]->Fill(q2kin);
					dHist_q2kin_varWidth_MassRegions[i][j][k]->Fill(q2kin);
					if(fillQ2ResSystematics)
						dHist_qvec2_res_vs_q2kin_MassRegions[i][j][k]->Fill(q2kin, q2kinRes);
				}
			}
			conditions.minWkin = activeFiducialConditions.minWkin;
			conditions.maxWkin = activeFiducialConditions.maxWkin;

			if(enableChiSqSystematics) {
				conditions.applyChiSqCut = true;
				for (int k = 0; k < kNumChiSqCuts; k++) {
					conditions.maxKFChiSq = dMaxChiSq[k];
					setJTphiBeamWindow(conditions);
					if (ComboPassesCuts(inputs, conditions))
						dHist_JTphi_MaxChiSq[i][j][k]->Fill(JTphi);
					setQ2BeamWindow(conditions);
					if (ComboPassesCuts(inputs, conditions)) {
						dHist_q2kin_MaxChiSq[i][j][k]->Fill(q2kin);
						dHist_q2kin_varWidth_MaxChiSq[i][j][k]->Fill(q2kin);
						if(fillQ2ResSystematics)
							dHist_qvec2_res_vs_q2kin_MaxChiSq[i][j][k]->Fill(q2kin, q2kinRes);
					}
				}
				conditions.applyChiSqCut = activeFiducialConditions.applyChiSqCut;
				conditions.maxKFChiSq = activeFiducialConditions.maxKFChiSq;
			}
		}
	}

	CutConditions mvaSystConds = activeFiducialConditions;
	mvaSystConds.applyMVACuts = true;
	mvaSystConds.modelChoice = "MLP";
	mvaSystConds.particleChoice = "ee";
	for (int i = 0; i < kNumMLPeeThresholds; i++) {
		mvaSystConds.MLP_ee = dMLPeeThresholds[i];
		setJTphiBeamWindow(mvaSystConds);
		if (ComboPassesCuts(inputs, mvaSystConds)) {
			dHist_JTphi_MLP_ee[i]->Fill(JTphi);
		}
		setQ2BeamWindow(mvaSystConds);
		if (ComboPassesCuts(inputs, mvaSystConds)) {
			dHist_q2kin_MLP_ee[i]->Fill(q2kin);
			dHist_q2kin_varWidth_MLP_ee[i]->Fill(q2kin);
			if(fillQ2ResSystematics)
				dHist_qvec2_res_vs_q2kin_MLP_ee[i]->Fill(q2kin, q2kinRes);
		}
	}

	mvaSystConds.particleChoice = "pi";
	for (int i = 0; i < kNumMLPpiThresholds; i++) {
		mvaSystConds.MLP_pi = dMLPpiThresholds[i];
		setJTphiBeamWindow(mvaSystConds);
		if (ComboPassesCuts(inputs, mvaSystConds)) {
			dHist_JTphi_MLP_pi[i]->Fill(JTphi);
		}
		setQ2BeamWindow(mvaSystConds);
		if (ComboPassesCuts(inputs, mvaSystConds)) {
			dHist_q2kin_MLP_pi[i]->Fill(q2kin);
			dHist_q2kin_varWidth_MLP_pi[i]->Fill(q2kin);
			if(fillQ2ResSystematics)
				dHist_qvec2_res_vs_q2kin_MLP_pi[i]->Fill(q2kin, q2kinRes);
		}
	}

	mvaSystConds.modelChoice = "BDT";
	for (int i = 0; i < kNumBDTpiThresholds; i++) {
		mvaSystConds.BDT_pi = dBDTpiThresholds[i];
		setJTphiBeamWindow(mvaSystConds);
		if (ComboPassesCuts(inputs, mvaSystConds)) {
			dHist_JTphi_BDT_pi[i]->Fill(JTphi);
		}
		setQ2BeamWindow(mvaSystConds);
		if (ComboPassesCuts(inputs, mvaSystConds)) {
			dHist_q2kin_BDT_pi[i]->Fill(q2kin);
			dHist_q2kin_varWidth_BDT_pi[i]->Fill(q2kin);
			if(fillQ2ResSystematics)
				dHist_qvec2_res_vs_q2kin_BDT_pi[i]->Fill(q2kin, q2kinRes);
		}
	}

	mvaSystConds.particleChoice = "ee";
	for (int i = 0; i < kNumBDTeeThresholds; i++) {
		mvaSystConds.BDT_ee = dBDTeeThresholds[i];
		setJTphiBeamWindow(mvaSystConds);
		if (ComboPassesCuts(inputs, mvaSystConds)) {
			dHist_JTphi_BDT_ee[i]->Fill(JTphi);
		}
		setQ2BeamWindow(mvaSystConds);
		if (ComboPassesCuts(inputs, mvaSystConds)) {
			dHist_q2kin_BDT_ee[i]->Fill(q2kin);
			dHist_q2kin_varWidth_BDT_ee[i]->Fill(q2kin);
			if(fillQ2ResSystematics)
				dHist_qvec2_res_vs_q2kin_BDT_ee[i]->Fill(q2kin, q2kinRes);
		}
	}

	dHist_MLP_responsePIP->Fill(inputs.MLP1);
	dHist_MLP_responsePIM->Fill(inputs.MLP2);
	dHist_MLP_response_PIP_vs_PIM->Fill(inputs.MLP2, inputs.MLP1);

	dHist_BDT_responsePIP->Fill(inputs.BDT1);
	dHist_BDT_responsePIM->Fill(inputs.BDT2);
	dHist_BDT_response_PIP_vs_PIM->Fill(inputs.BDT2, inputs.BDT1);
}

/* Commented out - ApplyFCALCorrections not implemented in systematics selector
void DSelector_2eMissingProton_Systematics::ApplyFCALCorrections(Double_t positiveP, Double_t negativeP,
                                                      Double_t& FCAL_Energy_plus, Double_t& FCAL_Energy_minus)
{
	// cout << "DEBUG FCAL_FUNC: Inside ApplyFCALCorrections, positiveP=" << positiveP << ", negativeP=" << negativeP << endl;
	// Apply momentum-dependent non-linear FCAL energy corrections
	// Corrects MC simulation to match data E/p distributions
	//
	// For positive track (e+):
	if(positiveP <= 3.60){
		FCAL_Energy_plus -= ep_FCAL_cor1->Eval(positiveP);
		FCAL_Energy_plus += ep_FCAL_cor2data1->Eval(positiveP);
	} else {
		FCAL_Energy_plus -= ep_FCAL_cor2->Eval(positiveP);
		FCAL_Energy_plus += ep_FCAL_cor2data2->Eval(positiveP);
	}
	
	// For negative track (e-):
	if(negativeP <= 3.60){
		FCAL_Energy_minus -= em_FCAL_cor1->Eval(negativeP);
		FCAL_Energy_minus += em_FCAL_cor2data1->Eval(negativeP);
	} else {
		FCAL_Energy_minus -= em_FCAL_cor2->Eval(negativeP);
		FCAL_Energy_minus += em_FCAL_cor2data2->Eval(negativeP);
	}
}
*/

int DSelector_2eMissingProton_Systematics::GetBeamWindowIndex(double beamE)
{
	// Select energy ranges based on experiment/run period
	double fullSpectrumMin, fullSpectrumMax, coherentPeakMin, coherentPeakMax;
	
	if (dIsCPPRunPeriod) {
		// CPP (Lead target): Runs 110621-112001
		fullSpectrumMin = 4.0;
		fullSpectrumMax = 7.6;
		coherentPeakMin = 5.35;
		coherentPeakMax = 5.75;
	} else {
		// GlueX-I (Diamond target)
		fullSpectrumMin = 7.0;
		fullSpectrumMax = 11.4;
		coherentPeakMin = 8.2;
		coherentPeakMax = 8.8;
	}
	
	// Check coherent peak first (narrower range)
	if (beamE >= coherentPeakMin && beamE <= coherentPeakMax) {
		return 1;  // CoherentPeak
	}
	// Check full spectrum
	else if (beamE >= fullSpectrumMin && beamE <= fullSpectrumMax) {
		return 0;  // FullSpectrum
	}
	// Outside both ranges
	return -1;
}

static bool EqualsIgnoreCase(const TString& a, const char* b)
{
	TString aa = a;
	aa.ToUpper();
	TString bb = b;
	bb.ToUpper();
	return aa == bb;
}

static bool IsAutoOrBlank(const TString& value)
{
	if(value.IsNull() || value == "")
		return true;
	return EqualsIgnoreCase(value, "auto");
}

static bool IsTokenMode(const TString& value)
{
	return EqualsIgnoreCase(value, "tokens") || EqualsIgnoreCase(value, "token") || EqualsIgnoreCase(value, "autotag") || EqualsIgnoreCase(value, "parsed");
}

// Token-based helpers retained for optional richer auto-tagging logic.
static bool IsAllDigits(const TString& s)
{
	if(s.IsNull() || s.Length() == 0)
		return false;
	for(Int_t i = 0; i < s.Length(); ++i) {
		if(!std::isdigit((unsigned char)s[i]))
			return false;
	}
	return true;
}

static TString ToUpperToken(TString s)
{
	s.ToUpper();
	return s;
}

static TString ExtractBRunFileTag(TObjArray* tokens)
{
	for(Int_t i = 0; i + 2 < tokens->GetEntries(); ++i) {
		TString t0 = ((TObjString*)tokens->At(i))->GetString();
		TString t1 = ((TObjString*)tokens->At(i + 1))->GetString();
		TString t2 = ((TObjString*)tokens->At(i + 2))->GetString();

		TString t0u = ToUpperToken(t0);
		if(t0u.Length() >= 2 && t0u[0] == 'B') {
			TString bnum = t0u(1, t0u.Length() - 1);
			if(IsAllDigits(bnum) && IsAllDigits(t1) && IsAllDigits(t2) && t1.Length() == 6 && t2.Length() == 3)
				return t0u + "_" + t1 + "_" + t2;  // Return full tag with run/file numbers
		}
	}
	return "";
}

static TString BuildStructuredTagFromTokens(TObjArray* tokens)
{
	TString study = "";
	TString dataset = "";
	TString runPeriod = "";
	TString pol = "";
	TString formFactor = "";
	TString radiation = "";
	TString runNumber = "";
	TString fileNo = "";

	for(Int_t i = 0; i < tokens->GetEntries(); ++i) {
		TString tok = ((TObjString*)tokens->At(i))->GetString();
		if(tok == "")
			continue;
		TString up = ToUpperToken(tok);

		if(study == "" && (up.BeginsWith("FFS") || up.BeginsWith("PS")))
			study = up;
		if(dataset == "" && (up == "BCFFN" || up == "BCFF1" || up == "QDATAQ" || up == "SIM"))
			dataset = up;
		if(runPeriod == "" && (up == "1801" || up == "1808" || up == "2205"))
			runPeriod = up;
		if(formFactor == "" && (up == "FF1" || up == "FFN"))
			formFactor = up;
		if(radiation == "" && (up == "RADOFF" || up == "NORAD" || up == "ONERAD" || up == "SINGLERAD" || up == "DBLRAD"))
			radiation = up;
		if(pol == "") {
			if(up == "AMO")
				pol = "AMO";
			else if(up.EndsWith("DEG")) {
				TString deg = up(0, up.Length() - 3);
				if(deg == "0" || deg == "45" || deg == "90" || deg == "135")
					pol = deg + "DEG";
			} else if(up == "0" || up == "45" || up == "90" || up == "135")
				pol = up + "DEG";
		}

		if(IsAllDigits(up)) {
			if(runNumber == "" && up.Length() == 6)
				runNumber = up;
			else if(fileNo == "" && up.Length() == 3)
				fileNo = up;
		}
	}

	if(runPeriod == "" && runNumber != "") {
		int runNum = runNumber.Atoi();
		if(runNum >= 40856 && runNum <= 42550)
			runPeriod = "1801";
		else if(runNum >= 50685 && runNum <= 51768)
			runPeriod = "1808";
		else if(runNum >= 110621 && runNum <= 112001)
			runPeriod = "CPP";
	}

	if(study == "")
		return "";
	if(runNumber == "" || fileNo == "")
		return "";

	TString tag = study;
	if(study.BeginsWith("FFS") && dataset != "")
		tag += "_" + dataset;
	if(runPeriod != "")
		tag += "_" + runPeriod;
	if(pol != "")
		tag += "_" + pol;
	if(formFactor != "")
		tag += "_" + formFactor;
	if(radiation != "")
		tag += "_" + radiation;
	tag += "_" + runNumber + "_" + fileNo;
	return tag;
}

void DSelector_2eMissingProton_Systematics::Finalize(void)
{
	//Save anything to output here that you do not want to be in the default DSelector output ROOT file.

	//Otherwise, don't do anything else (especially if you are using PROOF).
		//If you are using PROOF, this function is called on each thread,
		//so anything you do will not have the combined information from the various threads.
		//Besides, it is best-practice to do post-processing (e.g. fitting) separately, in case there is a problem.

	//DO YOUR STUFF HERE
	
	// Sum Plus and Minus histograms into Both histograms
	// dFCALvsTheta histograms
	if (dFCALvsTheta.Energy[kPlus] && dFCALvsTheta.Energy[kMinus] && dFCALvsTheta.Energy[kBoth]) {
		dFCALvsTheta.Energy[kBoth]->Add(dFCALvsTheta.Energy[kPlus]);
		dFCALvsTheta.Energy[kBoth]->Add(dFCALvsTheta.Energy[kMinus]);
	}
	if (dFCALvsTheta.EoverP[kPlus] && dFCALvsTheta.EoverP[kMinus] && dFCALvsTheta.EoverP[kBoth]) {
		dFCALvsTheta.EoverP[kBoth]->Add(dFCALvsTheta.EoverP[kPlus]);
		dFCALvsTheta.EoverP[kBoth]->Add(dFCALvsTheta.EoverP[kMinus]);
	}
	if (dFCALvsTheta.EoverPmeas[kPlus] && dFCALvsTheta.EoverPmeas[kMinus] && dFCALvsTheta.EoverPmeas[kBoth]) {
		dFCALvsTheta.EoverPmeas[kBoth]->Add(dFCALvsTheta.EoverPmeas[kPlus]);
		dFCALvsTheta.EoverPmeas[kBoth]->Add(dFCALvsTheta.EoverPmeas[kMinus]);
	}
	if (dFCALvsTheta.DeltaEfcal_kinfitE[kPlus] && dFCALvsTheta.DeltaEfcal_kinfitE[kMinus] && dFCALvsTheta.DeltaEfcal_kinfitE[kBoth]) {
		dFCALvsTheta.DeltaEfcal_kinfitE[kBoth]->Add(dFCALvsTheta.DeltaEfcal_kinfitE[kPlus]);
		dFCALvsTheta.DeltaEfcal_kinfitE[kBoth]->Add(dFCALvsTheta.DeltaEfcal_kinfitE[kMinus]);
	}
	if (dFCALvsTheta.TrackDOCA[kPlus] && dFCALvsTheta.TrackDOCA[kMinus] && dFCALvsTheta.TrackDOCA[kBoth]) {
		dFCALvsTheta.TrackDOCA[kBoth]->Add(dFCALvsTheta.TrackDOCA[kPlus]);
		dFCALvsTheta.TrackDOCA[kBoth]->Add(dFCALvsTheta.TrackDOCA[kMinus]);
	}
	if (dFCALvsTheta.E1E9[kPlus] && dFCALvsTheta.E1E9[kMinus] && dFCALvsTheta.E1E9[kBoth]) {
		dFCALvsTheta.E1E9[kBoth]->Add(dFCALvsTheta.E1E9[kPlus]);
		dFCALvsTheta.E1E9[kBoth]->Add(dFCALvsTheta.E1E9[kMinus]);
	}
	if (dFCALvsTheta.E9E25[kPlus] && dFCALvsTheta.E9E25[kMinus] && dFCALvsTheta.E9E25[kBoth]) {
		dFCALvsTheta.E9E25[kBoth]->Add(dFCALvsTheta.E9E25[kPlus]);
		dFCALvsTheta.E9E25[kBoth]->Add(dFCALvsTheta.E9E25[kMinus]);
	}
	if (dFCALvsTheta.KinRes[kPlus] && dFCALvsTheta.KinRes[kMinus] && dFCALvsTheta.KinRes[kBoth]) {
		dFCALvsTheta.KinRes[kBoth]->Add(dFCALvsTheta.KinRes[kPlus]);
		dFCALvsTheta.KinRes[kBoth]->Add(dFCALvsTheta.KinRes[kMinus]);
	}
	if (dFCALvsTheta.MeasRes[kPlus] && dFCALvsTheta.MeasRes[kMinus] && dFCALvsTheta.MeasRes[kBoth]) {
		dFCALvsTheta.MeasRes[kBoth]->Add(dFCALvsTheta.MeasRes[kPlus]);
		dFCALvsTheta.MeasRes[kBoth]->Add(dFCALvsTheta.MeasRes[kMinus]);
	}
	if (dFCALvsTheta.Saturation[kPlus] && dFCALvsTheta.Saturation[kMinus] && dFCALvsTheta.Saturation[kBoth]) {
		dFCALvsTheta.Saturation[kBoth]->Add(dFCALvsTheta.Saturation[kPlus]);
		dFCALvsTheta.Saturation[kBoth]->Add(dFCALvsTheta.Saturation[kMinus]);
	}
	if (dFCALvsTheta.SumU[kPlus] && dFCALvsTheta.SumU[kMinus] && dFCALvsTheta.SumU[kBoth]) {
		dFCALvsTheta.SumU[kBoth]->Add(dFCALvsTheta.SumU[kPlus]);
		dFCALvsTheta.SumU[kBoth]->Add(dFCALvsTheta.SumU[kMinus]);
	}
	if (dFCALvsTheta.SumV[kPlus] && dFCALvsTheta.SumV[kMinus] && dFCALvsTheta.SumV[kBoth]) {
		dFCALvsTheta.SumV[kBoth]->Add(dFCALvsTheta.SumV[kPlus]);
		dFCALvsTheta.SumV[kBoth]->Add(dFCALvsTheta.SumV[kMinus]);
	}
	if (dFCALvsTheta.UVAsymmetry[kPlus] && dFCALvsTheta.UVAsymmetry[kMinus] && dFCALvsTheta.UVAsymmetry[kBoth]) {
		dFCALvsTheta.UVAsymmetry[kBoth]->Add(dFCALvsTheta.UVAsymmetry[kPlus]);
		dFCALvsTheta.UVAsymmetry[kBoth]->Add(dFCALvsTheta.UVAsymmetry[kMinus]);
	}
	if (dFCALvsTheta.Energy_PostCorr[kPlus] && dFCALvsTheta.Energy_PostCorr[kMinus] && dFCALvsTheta.Energy_PostCorr[kBoth]) {
		dFCALvsTheta.Energy_PostCorr[kBoth]->Add(dFCALvsTheta.Energy_PostCorr[kPlus]);
		dFCALvsTheta.Energy_PostCorr[kBoth]->Add(dFCALvsTheta.Energy_PostCorr[kMinus]);
	}
	if (dFCALvsTheta.EoverP_PostCorr[kPlus] && dFCALvsTheta.EoverP_PostCorr[kMinus] && dFCALvsTheta.EoverP_PostCorr[kBoth]) {
		dFCALvsTheta.EoverP_PostCorr[kBoth]->Add(dFCALvsTheta.EoverP_PostCorr[kPlus]);
		dFCALvsTheta.EoverP_PostCorr[kBoth]->Add(dFCALvsTheta.EoverP_PostCorr[kMinus]);
	}
	if (dFCALvsTheta.EoverPmeas_PostCorr[kPlus] && dFCALvsTheta.EoverPmeas_PostCorr[kMinus] && dFCALvsTheta.EoverPmeas_PostCorr[kBoth]) {
		dFCALvsTheta.EoverPmeas_PostCorr[kBoth]->Add(dFCALvsTheta.EoverPmeas_PostCorr[kPlus]);
		dFCALvsTheta.EoverPmeas_PostCorr[kBoth]->Add(dFCALvsTheta.EoverPmeas_PostCorr[kMinus]);
	}
	if (dFCALvsTheta.DeltaEfcal_kinfitE_PostCorr[kPlus] && dFCALvsTheta.DeltaEfcal_kinfitE_PostCorr[kMinus] && dFCALvsTheta.DeltaEfcal_kinfitE_PostCorr[kBoth]) {
		dFCALvsTheta.DeltaEfcal_kinfitE_PostCorr[kBoth]->Add(dFCALvsTheta.DeltaEfcal_kinfitE_PostCorr[kPlus]);
		dFCALvsTheta.DeltaEfcal_kinfitE_PostCorr[kBoth]->Add(dFCALvsTheta.DeltaEfcal_kinfitE_PostCorr[kMinus]);
	}
	if (dFCALvsTheta.KinRes_PostCorr[kPlus] && dFCALvsTheta.KinRes_PostCorr[kMinus] && dFCALvsTheta.KinRes_PostCorr[kBoth]) {
		dFCALvsTheta.KinRes_PostCorr[kBoth]->Add(dFCALvsTheta.KinRes_PostCorr[kPlus]);
		dFCALvsTheta.KinRes_PostCorr[kBoth]->Add(dFCALvsTheta.KinRes_PostCorr[kMinus]);
	}
	if (dFCALvsTheta.MeasRes_PostCorr[kPlus] && dFCALvsTheta.MeasRes_PostCorr[kMinus] && dFCALvsTheta.MeasRes_PostCorr[kBoth]) {
		dFCALvsTheta.MeasRes_PostCorr[kBoth]->Add(dFCALvsTheta.MeasRes_PostCorr[kPlus]);
		dFCALvsTheta.MeasRes_PostCorr[kBoth]->Add(dFCALvsTheta.MeasRes_PostCorr[kMinus]);
	}
	if (dFCALvsTheta.Saturation_PostCorr[kPlus] && dFCALvsTheta.Saturation_PostCorr[kMinus] && dFCALvsTheta.Saturation_PostCorr[kBoth]) {
		dFCALvsTheta.Saturation_PostCorr[kBoth]->Add(dFCALvsTheta.Saturation_PostCorr[kPlus]);
		dFCALvsTheta.Saturation_PostCorr[kBoth]->Add(dFCALvsTheta.Saturation_PostCorr[kMinus]);
	}
	
	// dFCALvsPkin histograms
	if (dFCALvsPkin.Energy[kPlus] && dFCALvsPkin.Energy[kMinus] && dFCALvsPkin.Energy[kBoth]) {
		dFCALvsPkin.Energy[kBoth]->Add(dFCALvsPkin.Energy[kPlus]);
		dFCALvsPkin.Energy[kBoth]->Add(dFCALvsPkin.Energy[kMinus]);
	}
	if (dFCALvsPkin.EoverP[kPlus] && dFCALvsPkin.EoverP[kMinus] && dFCALvsPkin.EoverP[kBoth]) {
		dFCALvsPkin.EoverP[kBoth]->Add(dFCALvsPkin.EoverP[kPlus]);
		dFCALvsPkin.EoverP[kBoth]->Add(dFCALvsPkin.EoverP[kMinus]);
	}
	if (dFCALvsPkin.EoverPmeas[kPlus] && dFCALvsPkin.EoverPmeas[kMinus] && dFCALvsPkin.EoverPmeas[kBoth]) {
		dFCALvsPkin.EoverPmeas[kBoth]->Add(dFCALvsPkin.EoverPmeas[kPlus]);
		dFCALvsPkin.EoverPmeas[kBoth]->Add(dFCALvsPkin.EoverPmeas[kMinus]);
	}
	if (dFCALvsPkin.DeltaEfcal_kinfitE[kPlus] && dFCALvsPkin.DeltaEfcal_kinfitE[kMinus] && dFCALvsPkin.DeltaEfcal_kinfitE[kBoth]) {
		dFCALvsPkin.DeltaEfcal_kinfitE[kBoth]->Add(dFCALvsPkin.DeltaEfcal_kinfitE[kPlus]);
		dFCALvsPkin.DeltaEfcal_kinfitE[kBoth]->Add(dFCALvsPkin.DeltaEfcal_kinfitE[kMinus]);
	}
	if (dFCALvsPkin.TrackDOCA[kPlus] && dFCALvsPkin.TrackDOCA[kMinus] && dFCALvsPkin.TrackDOCA[kBoth]) {
		dFCALvsPkin.TrackDOCA[kBoth]->Add(dFCALvsPkin.TrackDOCA[kPlus]);
		dFCALvsPkin.TrackDOCA[kBoth]->Add(dFCALvsPkin.TrackDOCA[kMinus]);
	}
	if (dFCALvsPkin.E1E9[kPlus] && dFCALvsPkin.E1E9[kMinus] && dFCALvsPkin.E1E9[kBoth]) {
		dFCALvsPkin.E1E9[kBoth]->Add(dFCALvsPkin.E1E9[kPlus]);
		dFCALvsPkin.E1E9[kBoth]->Add(dFCALvsPkin.E1E9[kMinus]);
	}
	if (dFCALvsPkin.E9E25[kPlus] && dFCALvsPkin.E9E25[kMinus] && dFCALvsPkin.E9E25[kBoth]) {
		dFCALvsPkin.E9E25[kBoth]->Add(dFCALvsPkin.E9E25[kPlus]);
		dFCALvsPkin.E9E25[kBoth]->Add(dFCALvsPkin.E9E25[kMinus]);
	}
	if (dFCALvsPkin.KinRes[kPlus] && dFCALvsPkin.KinRes[kMinus] && dFCALvsPkin.KinRes[kBoth]) {
		dFCALvsPkin.KinRes[kBoth]->Add(dFCALvsPkin.KinRes[kPlus]);
		dFCALvsPkin.KinRes[kBoth]->Add(dFCALvsPkin.KinRes[kMinus]);
	}
	if (dFCALvsPkin.MeasRes[kPlus] && dFCALvsPkin.MeasRes[kMinus] && dFCALvsPkin.MeasRes[kBoth]) {
		dFCALvsPkin.MeasRes[kBoth]->Add(dFCALvsPkin.MeasRes[kPlus]);
		dFCALvsPkin.MeasRes[kBoth]->Add(dFCALvsPkin.MeasRes[kMinus]);
	}
	if (dFCALvsPkin.Saturation[kPlus] && dFCALvsPkin.Saturation[kMinus] && dFCALvsPkin.Saturation[kBoth]) {
		dFCALvsPkin.Saturation[kBoth]->Add(dFCALvsPkin.Saturation[kPlus]);
		dFCALvsPkin.Saturation[kBoth]->Add(dFCALvsPkin.Saturation[kMinus]);
	}
	if (dFCALvsPkin.SumU[kPlus] && dFCALvsPkin.SumU[kMinus] && dFCALvsPkin.SumU[kBoth]) {
		dFCALvsPkin.SumU[kBoth]->Add(dFCALvsPkin.SumU[kPlus]);
		dFCALvsPkin.SumU[kBoth]->Add(dFCALvsPkin.SumU[kMinus]);
	}
	if (dFCALvsPkin.SumV[kPlus] && dFCALvsPkin.SumV[kMinus] && dFCALvsPkin.SumV[kBoth]) {
		dFCALvsPkin.SumV[kBoth]->Add(dFCALvsPkin.SumV[kPlus]);
		dFCALvsPkin.SumV[kBoth]->Add(dFCALvsPkin.SumV[kMinus]);
	}
	if (dFCALvsPkin.UVAsymmetry[kPlus] && dFCALvsPkin.UVAsymmetry[kMinus] && dFCALvsPkin.UVAsymmetry[kBoth]) {
		dFCALvsPkin.UVAsymmetry[kBoth]->Add(dFCALvsPkin.UVAsymmetry[kPlus]);
		dFCALvsPkin.UVAsymmetry[kBoth]->Add(dFCALvsPkin.UVAsymmetry[kMinus]);
	}
	if (dFCALvsPkin.Energy_PostCorr[kPlus] && dFCALvsPkin.Energy_PostCorr[kMinus] && dFCALvsPkin.Energy_PostCorr[kBoth]) {
		dFCALvsPkin.Energy_PostCorr[kBoth]->Add(dFCALvsPkin.Energy_PostCorr[kPlus]);
		dFCALvsPkin.Energy_PostCorr[kBoth]->Add(dFCALvsPkin.Energy_PostCorr[kMinus]);
	}
	if (dFCALvsPkin.EoverP_PostCorr[kPlus] && dFCALvsPkin.EoverP_PostCorr[kMinus] && dFCALvsPkin.EoverP_PostCorr[kBoth]) {
		dFCALvsPkin.EoverP_PostCorr[kBoth]->Add(dFCALvsPkin.EoverP_PostCorr[kPlus]);
		dFCALvsPkin.EoverP_PostCorr[kBoth]->Add(dFCALvsPkin.EoverP_PostCorr[kMinus]);
	}
	if (dFCALvsPkin.EoverPmeas_PostCorr[kPlus] && dFCALvsPkin.EoverPmeas_PostCorr[kMinus] && dFCALvsPkin.EoverPmeas_PostCorr[kBoth]) {
		dFCALvsPkin.EoverPmeas_PostCorr[kBoth]->Add(dFCALvsPkin.EoverPmeas_PostCorr[kPlus]);
		dFCALvsPkin.EoverPmeas_PostCorr[kBoth]->Add(dFCALvsPkin.EoverPmeas_PostCorr[kMinus]);
	}
	if (dFCALvsPkin.DeltaEfcal_kinfitE_PostCorr[kPlus] && dFCALvsPkin.DeltaEfcal_kinfitE_PostCorr[kMinus] && dFCALvsPkin.DeltaEfcal_kinfitE_PostCorr[kBoth]) {
		dFCALvsPkin.DeltaEfcal_kinfitE_PostCorr[kBoth]->Add(dFCALvsPkin.DeltaEfcal_kinfitE_PostCorr[kPlus]);
		dFCALvsPkin.DeltaEfcal_kinfitE_PostCorr[kBoth]->Add(dFCALvsPkin.DeltaEfcal_kinfitE_PostCorr[kMinus]);
	}
	if (dFCALvsPkin.KinRes_PostCorr[kPlus] && dFCALvsPkin.KinRes_PostCorr[kMinus] && dFCALvsPkin.KinRes_PostCorr[kBoth]) {
		dFCALvsPkin.KinRes_PostCorr[kBoth]->Add(dFCALvsPkin.KinRes_PostCorr[kPlus]);
		dFCALvsPkin.KinRes_PostCorr[kBoth]->Add(dFCALvsPkin.KinRes_PostCorr[kMinus]);
	}
	if (dFCALvsPkin.MeasRes_PostCorr[kPlus] && dFCALvsPkin.MeasRes_PostCorr[kMinus] && dFCALvsPkin.MeasRes_PostCorr[kBoth]) {
		dFCALvsPkin.MeasRes_PostCorr[kBoth]->Add(dFCALvsPkin.MeasRes_PostCorr[kPlus]);
		dFCALvsPkin.MeasRes_PostCorr[kBoth]->Add(dFCALvsPkin.MeasRes_PostCorr[kMinus]);
	}
	if (dFCALvsPkin.Saturation_PostCorr[kPlus] && dFCALvsPkin.Saturation_PostCorr[kMinus] && dFCALvsPkin.Saturation_PostCorr[kBoth]) {
		dFCALvsPkin.Saturation_PostCorr[kBoth]->Add(dFCALvsPkin.Saturation_PostCorr[kPlus]);
		dFCALvsPkin.Saturation_PostCorr[kBoth]->Add(dFCALvsPkin.Saturation_PostCorr[kMinus]);
	}

	// dFCALvsqvec2 histograms
	if (dFCALvsqvec2.Energy[kPlus] && dFCALvsqvec2.Energy[kMinus] && dFCALvsqvec2.Energy[kBoth]) {
		dFCALvsqvec2.Energy[kBoth]->Add(dFCALvsqvec2.Energy[kPlus]);
		dFCALvsqvec2.Energy[kBoth]->Add(dFCALvsqvec2.Energy[kMinus]);
	}
	if (dFCALvsqvec2.EoverP[kPlus] && dFCALvsqvec2.EoverP[kMinus] && dFCALvsqvec2.EoverP[kBoth]) {
		dFCALvsqvec2.EoverP[kBoth]->Add(dFCALvsqvec2.EoverP[kPlus]);
		dFCALvsqvec2.EoverP[kBoth]->Add(dFCALvsqvec2.EoverP[kMinus]);
	}
	if (dFCALvsqvec2.EoverPmeas[kPlus] && dFCALvsqvec2.EoverPmeas[kMinus] && dFCALvsqvec2.EoverPmeas[kBoth]) {
		dFCALvsqvec2.EoverPmeas[kBoth]->Add(dFCALvsqvec2.EoverPmeas[kPlus]);
		dFCALvsqvec2.EoverPmeas[kBoth]->Add(dFCALvsqvec2.EoverPmeas[kMinus]);
	}
	if (dFCALvsqvec2.DeltaEfcal_kinfitE[kPlus] && dFCALvsqvec2.DeltaEfcal_kinfitE[kMinus] && dFCALvsqvec2.DeltaEfcal_kinfitE[kBoth]) {
		dFCALvsqvec2.DeltaEfcal_kinfitE[kBoth]->Add(dFCALvsqvec2.DeltaEfcal_kinfitE[kPlus]);
		dFCALvsqvec2.DeltaEfcal_kinfitE[kBoth]->Add(dFCALvsqvec2.DeltaEfcal_kinfitE[kMinus]);
	}
	if (dFCALvsqvec2.TrackDOCA[kPlus] && dFCALvsqvec2.TrackDOCA[kMinus] && dFCALvsqvec2.TrackDOCA[kBoth]) {
		dFCALvsqvec2.TrackDOCA[kBoth]->Add(dFCALvsqvec2.TrackDOCA[kPlus]);
		dFCALvsqvec2.TrackDOCA[kBoth]->Add(dFCALvsqvec2.TrackDOCA[kMinus]);
	}
	if (dFCALvsqvec2.E1E9[kPlus] && dFCALvsqvec2.E1E9[kMinus] && dFCALvsqvec2.E1E9[kBoth]) {
		dFCALvsqvec2.E1E9[kBoth]->Add(dFCALvsqvec2.E1E9[kPlus]);
		dFCALvsqvec2.E1E9[kBoth]->Add(dFCALvsqvec2.E1E9[kMinus]);
	}
	if (dFCALvsqvec2.E9E25[kPlus] && dFCALvsqvec2.E9E25[kMinus] && dFCALvsqvec2.E9E25[kBoth]) {
		dFCALvsqvec2.E9E25[kBoth]->Add(dFCALvsqvec2.E9E25[kPlus]);
		dFCALvsqvec2.E9E25[kBoth]->Add(dFCALvsqvec2.E9E25[kMinus]);
	}
	if (dFCALvsqvec2.KinRes[kPlus] && dFCALvsqvec2.KinRes[kMinus] && dFCALvsqvec2.KinRes[kBoth]) {
		dFCALvsqvec2.KinRes[kBoth]->Add(dFCALvsqvec2.KinRes[kPlus]);
		dFCALvsqvec2.KinRes[kBoth]->Add(dFCALvsqvec2.KinRes[kMinus]);
	}
	if (dFCALvsqvec2.MeasRes[kPlus] && dFCALvsqvec2.MeasRes[kMinus] && dFCALvsqvec2.MeasRes[kBoth]) {
		dFCALvsqvec2.MeasRes[kBoth]->Add(dFCALvsqvec2.MeasRes[kPlus]);
		dFCALvsqvec2.MeasRes[kBoth]->Add(dFCALvsqvec2.MeasRes[kMinus]);
	}
	if (dFCALvsqvec2.Saturation[kPlus] && dFCALvsqvec2.Saturation[kMinus] && dFCALvsqvec2.Saturation[kBoth]) {
		dFCALvsqvec2.Saturation[kBoth]->Add(dFCALvsqvec2.Saturation[kPlus]);
		dFCALvsqvec2.Saturation[kBoth]->Add(dFCALvsqvec2.Saturation[kMinus]);
	}
	if (dFCALvsqvec2.SumU[kPlus] && dFCALvsqvec2.SumU[kMinus] && dFCALvsqvec2.SumU[kBoth]) {
		dFCALvsqvec2.SumU[kBoth]->Add(dFCALvsqvec2.SumU[kPlus]);
		dFCALvsqvec2.SumU[kBoth]->Add(dFCALvsqvec2.SumU[kMinus]);
	}
	if (dFCALvsqvec2.SumV[kPlus] && dFCALvsqvec2.SumV[kMinus] && dFCALvsqvec2.SumV[kBoth]) {
		dFCALvsqvec2.SumV[kBoth]->Add(dFCALvsqvec2.SumV[kPlus]);
		dFCALvsqvec2.SumV[kBoth]->Add(dFCALvsqvec2.SumV[kMinus]);
	}
	if (dFCALvsqvec2.UVAsymmetry[kPlus] && dFCALvsqvec2.UVAsymmetry[kMinus] && dFCALvsqvec2.UVAsymmetry[kBoth]) {
		dFCALvsqvec2.UVAsymmetry[kBoth]->Add(dFCALvsqvec2.UVAsymmetry[kPlus]);
		dFCALvsqvec2.UVAsymmetry[kBoth]->Add(dFCALvsqvec2.UVAsymmetry[kMinus]);
	}

	if (dFCALvsqvec2.Energy_PostCorr[kPlus] && dFCALvsqvec2.Energy_PostCorr[kMinus] && dFCALvsqvec2.Energy_PostCorr[kBoth]) {
		dFCALvsqvec2.Energy_PostCorr[kBoth]->Add(dFCALvsqvec2.Energy_PostCorr[kPlus]);
		dFCALvsqvec2.Energy_PostCorr[kBoth]->Add(dFCALvsqvec2.Energy_PostCorr[kMinus]);
	}
	if (dFCALvsqvec2.EoverP_PostCorr[kPlus] && dFCALvsqvec2.EoverP_PostCorr[kMinus] && dFCALvsqvec2.EoverP_PostCorr[kBoth]) {
		dFCALvsqvec2.EoverP_PostCorr[kBoth]->Add(dFCALvsqvec2.EoverP_PostCorr[kPlus]);
		dFCALvsqvec2.EoverP_PostCorr[kBoth]->Add(dFCALvsqvec2.EoverP_PostCorr[kMinus]);
	}
	if (dFCALvsqvec2.EoverPmeas_PostCorr[kPlus] && dFCALvsqvec2.EoverPmeas_PostCorr[kMinus] && dFCALvsqvec2.EoverPmeas_PostCorr[kBoth]) {
		dFCALvsqvec2.EoverPmeas_PostCorr[kBoth]->Add(dFCALvsqvec2.EoverPmeas_PostCorr[kPlus]);
		dFCALvsqvec2.EoverPmeas_PostCorr[kBoth]->Add(dFCALvsqvec2.EoverPmeas_PostCorr[kMinus]);
	}
	if (dFCALvsqvec2.DeltaEfcal_kinfitE_PostCorr[kPlus] && dFCALvsqvec2.DeltaEfcal_kinfitE_PostCorr[kMinus] && dFCALvsqvec2.DeltaEfcal_kinfitE_PostCorr[kBoth]) {
		dFCALvsqvec2.DeltaEfcal_kinfitE_PostCorr[kBoth]->Add(dFCALvsqvec2.DeltaEfcal_kinfitE_PostCorr[kPlus]);
		dFCALvsqvec2.DeltaEfcal_kinfitE_PostCorr[kBoth]->Add(dFCALvsqvec2.DeltaEfcal_kinfitE_PostCorr[kMinus]);
	}
	if (dFCALvsqvec2.KinRes_PostCorr[kPlus] && dFCALvsqvec2.KinRes_PostCorr[kMinus] && dFCALvsqvec2.KinRes_PostCorr[kBoth]) {
		dFCALvsqvec2.KinRes_PostCorr[kBoth]->Add(dFCALvsqvec2.KinRes_PostCorr[kPlus]);
		dFCALvsqvec2.KinRes_PostCorr[kBoth]->Add(dFCALvsqvec2.KinRes_PostCorr[kMinus]);
	}
	if (dFCALvsqvec2.MeasRes_PostCorr[kPlus] && dFCALvsqvec2.MeasRes_PostCorr[kMinus] && dFCALvsqvec2.MeasRes_PostCorr[kBoth]) {
		dFCALvsqvec2.MeasRes_PostCorr[kBoth]->Add(dFCALvsqvec2.MeasRes_PostCorr[kPlus]);
		dFCALvsqvec2.MeasRes_PostCorr[kBoth]->Add(dFCALvsqvec2.MeasRes_PostCorr[kMinus]);
	}
	if (dFCALvsqvec2.Saturation_PostCorr[kPlus] && dFCALvsqvec2.Saturation_PostCorr[kMinus] && dFCALvsqvec2.Saturation_PostCorr[kBoth]) {
		dFCALvsqvec2.Saturation_PostCorr[kBoth]->Add(dFCALvsqvec2.Saturation_PostCorr[kPlus]);
		dFCALvsqvec2.Saturation_PostCorr[kBoth]->Add(dFCALvsqvec2.Saturation_PostCorr[kMinus]);
	}
	
	// Clean up TMVA Readers
	if(dTMVAReader_MLP != NULL) {
		delete dTMVAReader_MLP;
		dTMVAReader_MLP = NULL;
	}
	if(dTMVAReader_BDT != NULL) {
		delete dTMVAReader_BDT;
		dTMVAReader_BDT = NULL;
	}
	
	//CALL THIS LAST
	DSelector::Finalize(); //Saves results to the output file
}
