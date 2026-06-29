#include "DSelector_pippimmissp_psubMatrix.h"
#include "TMath.h"
#include <cctype>

static inline Float_t sanitize_f(Float_t v, Float_t sentinel = -999.f) {
    return std::isfinite((double)v) ? v : sentinel;
}


void DSelector_pippimmissp_psubMatrix::Init(TTree *locTree)
{
	// USERS: IN THIS FUNCTION, ONLY MODIFY SECTIONS WITH A "USER" OR "EXAMPLE" LABEL. LEAVE THE REST ALONE.

	// The Init() function is called when the selector needs to initialize a new tree or chain.
	// Typically here the branch addresses and branch pointers of the tree will be set.
	// Init() will be called many times when running on PROOF (once per file to be processed).

	//USERS: SET OUTPUT FILE NAME //can be overriden by user in PROOF
	dOutputFileName = "pippimmissp_psubMatrix.root"; //"" for none
	dOutputTreeFileName = ""; //"" for none
	dFlatTreeFileName = ""; //output flat tree (one combo per tree entry), "" for none
	dFlatTreeName = ""; //if blank, default name will be chosen
	//dSaveDefaultFlatBranches = true; // False: don't save default branches, reduce disk footprint.
	//dSaveTLorentzVectorsAsFundamentaFlatTree = false; // Default (or false): save particles as TLorentzVector objects. True: save as four doubles instead.

	//Because this function gets called for each TTree in the TChain, we must be careful:
		//We need to re-initialize the tree interface & branch wrappers, but don't want to recreate histograms
	bool locInitializedPriorFlag = dInitializedFlag; //save whether have been initialized previously
	DSelector::Init(locTree); //This must be called to initialize wrappers for each new TTree
	//gDirectory now points to the output file with name dOutputFileName (if any)
	if(locInitializedPriorFlag)
		return; //have already created histograms, etc. below: exit

	Get_ComboWrappers();
	dPreviousRunNumber = 0;

	/*********************************** EXAMPLE USER INITIALIZATION: ANALYSIS ACTIONS **********************************/

	// EXAMPLE: Create deque for histogramming particle masses:
	// // For histogramming the phi mass in phi -> K+ K-
	// // Be sure to change this and dAnalyzeCutActions to match reaction
	std::deque<Particle_t> MyPhi;
	MyPhi.push_back(KPlus); MyPhi.push_back(KMinus);

	//ANALYSIS ACTIONS: //Executed in order if added to dAnalysisActions
	//false/true below: use measured/kinfit data

	//PID
	dAnalysisActions.push_back(new DHistogramAction_ParticleID(dComboWrapper, false));
	//below: value: +/- N ns, UnknownParticle: All PIDs, SYS_NULL: all timing systems
	//dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 0.5, KPlus, SYS_BCAL));

	//PIDFOM (for charged tracks)
	dAnalysisActions.push_back(new DHistogramAction_PIDFOM(dComboWrapper));
	//dAnalysisActions.push_back(new DCutAction_PIDFOM(dComboWrapper, KPlus, 0.1));
	//dAnalysisActions.push_back(new DCutAction_EachPIDFOM(dComboWrapper, 0.1));

	//MASSES
	//dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, false, Lambda, 1000, 1.0, 1.2, "Lambda"));
	//dAnalysisActions.push_back(new DHistogramAction_MissingMassSquared(dComboWrapper, false, 1000, -0.1, 0.1));

	//KINFIT RESULTS
	dAnalysisActions.push_back(new DHistogramAction_KinFitResults(dComboWrapper));

	//CUT MISSING MASS
	//dAnalysisActions.push_back(new DCutAction_MissingMassSquared(dComboWrapper, false, -0.03, 0.02));

	//CUT ON SHOWER QUALITY
	//dAnalysisActions.push_back(new DCutAction_ShowerQuality(dComboWrapper, SYS_FCAL, 0.5));

	//BEAM ENERGY
	dAnalysisActions.push_back(new DHistogramAction_BeamEnergy(dComboWrapper, false));
	//dAnalysisActions.push_back(new DCutAction_BeamEnergy(dComboWrapper, false, 8.2, 8.8));  // Coherent peak for runs in the range 30000-59999

	//KINEMATICS
	dAnalysisActions.push_back(new DHistogramAction_ParticleComboKinematics(dComboWrapper, false));

	// ANALYZE CUT ACTIONS
	// // Change MyPhi to match reaction
	dAnalyzeCutActions = new DHistogramAction_AnalyzeCutActions( dAnalysisActions, dComboWrapper, false, 0, MyPhi, 1000, 0.9, 2.4, "CutActionEffect" );

	//INITIALIZE ACTIONS
	//If you create any actions that you want to run manually (i.e. don't add to dAnalysisActions), be sure to initialize them here as well
	Initialize_Actions();
	dAnalyzeCutActions->Initialize(); // manual action, must call Initialize()

	/******************************** EXAMPLE USER INITIALIZATION: STAND-ALONE HISTOGRAMS *******************************/

	//EXAMPLE MANUAL HISTOGRAMS:
	dHist_MissingMassSquared = new TH1I("MissingMassSquared", ";Missing Mass Squared (GeV/c^{2})^{2}", 600, -0.06, 0.06);
	dHist_BeamEnergy = new TH1I("BeamEnergy", ";Beam Energy (GeV)", 600, 0.0, 12.0);
	dHist_BeamEnergy_BestChiSq = new TH1I("BeamEnergy_BestChiSq", ";Beam Energy (GeV)", 600, 0.0, 12.0);
	dHist_InvMass_BestChiSq_total = new TH1D("InvMass_BestChiSq_total", ";Invariant Mass (GeV/c^{2});Counts", 100, 0.6, 0.8);
	dHist_InvMass_BestChiSq_EoPprecut = new TH1D("InvMass_BestChiSq_EoPprecut", ";Invariant Mass (GeV/c^{2});Counts", 100, 0.6, 0.8);
	dHist_MLP1_vs_MLP2_BestChiSq_total = new TH2D("MLP1_vs_MLP2_BestChiSq_total", "MLP1 vs MLP2 (Best ChiSq Combo);MLP1 Output;MLP2 Output", 100, 0.0, 1.0, 100, 0.0, 1.0);
	dHist_BDT1_vs_BDT2_BestChiSq_total = new TH2D("BDT1_vs_BDT2_BestChiSq_total", "BDT1 vs BDT2 (Best ChiSq Combo);BDT1 Output;BDT2 Output", 100, -1.0, 1.0, 100, -1.0, 1.0);

	dHist_MLP1_vs_MLP2_BestChiSq_EoPprecut = new TH2D("MLP1_vs_MLP2_BestChiSq_EoPprecut", "MLP1 vs MLP2 (Best ChiSq Combo);MLP1 Output;MLP2 Output", 100, 0.0, 1.0, 100, 0.0, 1.0);

	dHist_BDT1_vs_BDT2_BestChiSq_EoPprecut = new TH2D("BDT1_vs_BDT2_BestChiSq_EoPprecut", "BDT1 vs BDT2 (Best ChiSq Combo);BDT1 Output;BDT2 Output", 100, -1.0, 1.0, 100, -1.0, 1.0);
	/************************************** USER INITIALIZATION: CUT THRESHOLDS *************************************/

	dMinKinFitCL = 1e-6;
	dMinBeamEnergy = 7.0;
	dMaxBeamEnergy = 11.4;
	dMaxNumUnusedTracks = 0;

	// Cut flags: Set to true to enable, false to disable
	// Section 0: Preselection Cut Emulation:
	ApplyPreselectionEoPCut = false;         // 0a: E/p > 0.4 preselection (off by default)

	// Section 1: Exclusivity cuts
	ApplyNoUnusedTracksCut = true;           // 1a: No unused charged tracks
	ApplyNoUnusedShowersCut = true;          // 1b: No unused neutral shower energy (NOTE: Currently cuts if <= 0, may want to review)
	
	// Section 2: Forward detector cuts
	ApplyFCALEnergyNonZeroCut = false;       // 2a: FCAL E_1 and E_2 > 0, turn off for pion sub sample, need to first calc total rho0 yield
	ApplyTOFdEdxNonZeroCut = true;           // 2b: TOF dE/dx > 0 for both tracks
	
	// Section 3: Fiducial cuts (turn off for systematics studies)
	ApplyBeamEnergyCut = true;              // 3a: Coherent peak (8.2-8.8 GeV), 7-11.4 if FFStudy
	ApplyMinThetaCut = true;                // 3b: Minimum theta > 1.5 deg
	ApplyMaxThetaCut = false;               // 3c: Maximum theta cut (OFF - kills rho0 yield)
	Apply2DThetaAcceptanceMatchCut = false; // 3d: Sim/data acceptance match (OFF - kills rho0 yield)
	ApplyMomentumRangeCut = true;           // 3e: Valid FCAL correction range (0.45-7.92 GeV)
	ApplyInvariantMassCut = true;           // 3f: Invariant mass window (0.700-0.770 GeV for rho0)
	ApplyVertexZCut = true;                 // 3g: Vertex Z position (52-78 cm)
	ApplyKinFitCLCut = true;                // 3h: Kinematic fit confidence level > 1e-6

	OnlyBestChiSqComboInFlat = true;             // Only save combo with best chiSq to flat tree

	// =================================================================================
	// Additional test/development cuts (NOT used in official analysis)
	ApplyMinEoverPCut = false;              // Minimum E/P cut 
	ApplyMaxEoverPCut = false;              // Maximum E/P cut 
	ApplyMinFCALElasticityCut = false;      // Minimum FCAL elasticity
	ApplyMaxFCALElasticityCut = false;      // Maximum FCAL elasticity 
	ApplyMaxFCALDOCACut = false;            // Maximum FCAL DOCA cut 
	ApplyMaxTOFdEdxCut = false;             // Maximum TOF dE/dx 
	// Threshold values for test cuts (reference values from e+e- analysis)
	dMinEoverP = 0.4;                       // Minimum E/P (e+e- uses 0.3-0.7)
	dMaxEoverP = 1.2;                       // Maximum E/P 
	dMinFCALElasticity = 0.4;               // Minimum FCAL elasticity
	dMaxFCALElasticity = 1.2;               // Maximum FCAL elasticity (e+e- uses 1.2)
	dMaxFCALDOCA = 2.0;                     // Maximum FCAL DOCA in cm
	dMaxTOFdEdx = 2.0 * 0.00263;            // Maximum TOF dE/dx (2x mean from e+e-)
	// =================================================================================
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

	/************************** EXAMPLE USER INITIALIZATION: CUSTOM OUTPUT BRANCHES - MAIN TREE *************************/

	//EXAMPLE MAIN TREE CUSTOM BRANCHES (OUTPUT ROOT FILE NAME MUST FIRST BE GIVEN!!!! (ABOVE: TOP)):
	//The type for the branch must be included in the brackets
	//1st function argument is the name of the branch
	//2nd function argument is the name of the branch that contains the size of the array (for fundamentals only)
	/*
	dTreeInterface->Create_Branch_Fundamental<Int_t>("my_int"); //fundamental = char, int, float, double, etc.
	dTreeInterface->Create_Branch_FundamentalArray<Int_t>("my_int_array", "my_int");
	dTreeInterface->Create_Branch_FundamentalArray<Float_t>("my_combo_array", "NumCombos");
	dTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("my_p4");
	dTreeInterface->Create_Branch_ClonesArray<TLorentzVector>("my_p4_array");
	*/

	/************************** EXAMPLE USER INITIALIZATION: CUSTOM OUTPUT BRANCHES - FLAT TREE *************************/

	// RECOMMENDED: CREATE ACCIDENTAL WEIGHT BRANCH
	// dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("accidweight");

	//EXAMPLE FLAT TREE CUSTOM BRANCHES (OUTPUT ROOT FILE NAME MUST FIRST BE GIVEN!!!! (ABOVE: TOP)):
	//The type for the branch must be included in the brackets
	//1st function argument is the name of the branch
	//2nd function argument is the name of the branch that contains the size of the array (for fundamentals only)
	/*
	dFlatTreeInterface->Create_Branch_Fundamental<Int_t>("flat_my_int"); //fundamental = char, int, float, double, etc.
	dFlatTreeInterface->Create_Branch_FundamentalArray<Int_t>("flat_my_int_array", "flat_my_int");
	dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("flat_my_p4");
	dFlatTreeInterface->Create_Branch_ClonesArray<TLorentzVector>("flat_my_p4_array");
	*/

	/************************************* ADVANCED EXAMPLE: CHOOSE BRANCHES TO READ ************************************/

	//TO SAVE PROCESSING TIME
		//If you know you don't need all of the branches/data, but just a subset of it, you can speed things up
		//By default, for each event, the data is retrieved for all branches
		//If you know you only need data for some branches, you can skip grabbing data from the branches you don't need
		//Do this by doing something similar to the commented code below

	//dTreeInterface->Clear_GetEntryBranches(); //now get none
	//dTreeInterface->Register_GetEntryBranch("Proton__P4"); //manually set the branches you want

	/************************************** DETERMINE IF ANALYZING SIMULATED DATA *************************************/

	dIsMC = (dTreeInterface->Get_Branch("MCWeight") != NULL);

}

Bool_t DSelector_pippimmissp_psubMatrix::Process(Long64_t locEntry)
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
	//cout << "RUN " << Get_RunNumber() << ", EVENT " << Get_EventNumber() << endl;
	//TLorentzVector locProductionX4 = Get_X4_Production();

	/******************************************** GET POLARIZATION ORIENTATION ******************************************/

	//Only if the run number changes
	//RCDB environment must be setup in order for this to work! (Will return false otherwise)
	UInt_t locRunNumber = Get_RunNumber();
	if(locRunNumber != dPreviousRunNumber)
	{
		dIsPolarizedFlag = dAnalysisUtilities.Get_IsPolarizedBeam(locRunNumber, dIsPARAFlag);
		dPreviousRunNumber = locRunNumber;
	}

	/********************************************* SETUP UNIQUENESS TRACKING ********************************************/

	//ANALYSIS ACTIONS: Reset uniqueness tracking for each action
	//For any actions that you are executing manually, be sure to call Reset_NewEvent() on them here
	Reset_Actions_NewEvent();
	dAnalyzeCutActions->Reset_NewEvent(); // manual action, must call Reset_NewEvent()

	//PREVENT-DOUBLE COUNTING WHEN HISTOGRAMMING
		//Sometimes, some content is the exact same between one combo and the next
			//e.g. maybe two combos have different beam particles, but the same data for the final-state
		//When histogramming, you don't want to double-count when this happens: artificially inflates your signal (or background)
		//So, for each quantity you histogram, keep track of what particles you used (for a given combo)
		//Then for each combo, just compare to what you used before, and make sure it's unique

	//EXAMPLE 0: Event-specific info:
	Bool_t locUsedSoFar_Event = false; // Flag used to mark if the best chi-squared combo is filled in the histogram

	//EXAMPLE 1: Particle-specific info:
	set<Int_t> locUsedSoFar_BeamEnergy; //Int_t: Unique ID for beam particles. set: easy to use, fast to search. This container is used for the "hybrid" method dealing with combinatorics.

	//EXAMPLE 2: Combo-specific info:
		//In general: Could have multiple particles with the same PID: Use a set of Int_t's
		//In general: Multiple PIDs, so multiple sets: Contain within a map
		//Multiple combos: Contain maps within a set (easier, faster to search)
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_MissingMass;

	//INSERT USER ANALYSIS UNIQUENESS TRACKING HERE

	/**************************************** EXAMPLE: FILL CUSTOM OUTPUT BRANCHES **************************************/

	/*
	Int_t locMyInt = 7;
	dTreeInterface->Fill_Fundamental<Int_t>("my_int", locMyInt);

	TLorentzVector locMyP4(4.0, 3.0, 2.0, 1.0);
	dTreeInterface->Fill_TObject<TLorentzVector>("my_p4", locMyP4);

	for(int loc_i = 0; loc_i < locMyInt; ++loc_i)
		dTreeInterface->Fill_Fundamental<Int_t>("my_int_array", 3*loc_i, loc_i); //2nd argument = value, 3rd = array index
	*/

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
	// Sort the combos by ChiSq
	std::sort(loc_combos.begin(), loc_combos.end(), [](const std::pair<UInt_t, Double_t>& a, const std::pair<UInt_t, Double_t>& b) {
		return a.second < b.second;
	});

	// Track whether we've found the best in-time combo yet (for BestChiSq fills)
	bool foundBestInTimeCombo = false;
	//Loop over combos
	for(const auto& loc_combo : loc_combos)
	{
		UInt_t loc_i = loc_combo.first;
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
		Int_t locPiPlusTrackID = dPositiveWrapper->Get_TrackID();
		Int_t locPiMinusTrackID = dNegativeWrapper->Get_TrackID();

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

		TLorentzVector loc2TrackP4 = locPositiveP4 + locNegativeP4;
		TLorentzVector loc2TrackP4_Measured = locPositiveP4_Measured + locNegativeP4_Measured;

		double M2TrackKin = loc2TrackP4.M();  // Generic name for 2-track invariant mass

		/********************************************* GET COMBO RF TIMING INFO *****************************************/

		TLorentzVector locBeamX4_Measured = dComboBeamWrapper->Get_X4_Measured();
		TLorentzVector locBeam_X4_Measured = locBeamX4_Measured;
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

		/********************************************* COMBINE FOUR-MOMENTUM ********************************************/

		// DO YOUR STUFF HERE

		// Combine 4-vectors
		TLorentzVector locMissingP4_Measured = locBeamP4_Measured + dTargetP4;
		locMissingP4_Measured -= locPositiveP4_Measured + locNegativeP4_Measured;

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
		const Double_t tmva_FCAL_EoverPkin_plus = FCAL_EoverPkin_plus_hist;
		const Double_t tmva_FCAL_EoverPkin_minus = FCAL_EoverPkin_minus_hist;
		const Double_t tmva_FCAL_Saturation_plus = FCAL_Emax_plus;
		const Double_t tmva_FCAL_Saturation_minus = FCAL_Emax_minus;

        double FCAL_Asymmetry = fabs((FCAL_Energy_plus - FCAL_Energy_minus)/(FCAL_Energy_plus + FCAL_Energy_minus));
		// Geometric asymmetry from eq. 2.4: A_uv = |σ²_u - σ²_v| / |σ²_u + σ²_v|
		// SumU and SumV correspond to the second moments σ²_u and σ²_v
		double SumU_plus = dPositiveWrapper->Get_SumU_FCAL();
		double SumV_plus = dPositiveWrapper->Get_SumV_FCAL();
		double SumU_minus = dNegativeWrapper->Get_SumU_FCAL();
		double SumV_minus = dNegativeWrapper->Get_SumV_FCAL();
		double FCAL_UV_Asymmetry_plus = fabs(SumU_plus - SumV_plus) / fabs(SumU_plus + SumV_plus);
		double FCAL_UV_Asymmetry_minus = fabs(SumU_minus - SumV_minus) / fabs(SumU_minus + SumV_minus);

		// Calculate FCAL E/P variables for MLP/BDT inputs
		Double_t FCAL_EoverPkin_plus = FCAL_Energy_plus/locPositiveP4.Vect().Mag();
		Double_t FCAL_EoverPkin_minus = FCAL_Energy_minus/locNegativeP4.Vect().Mag();
		Double_t FCAL_EoverPmeas_plus = FCAL_Energy_plus/locPositiveP4_Measured.Vect().Mag();
		Double_t FCAL_EoverPmeas_minus = FCAL_Energy_minus/locNegativeP4_Measured.Vect().Mag();

		// dE/dx variables
		double TOF_dEdx_plus = dPositiveWrapper->Get_dEdx_TOF();
		double TOF_dEdx_minus = dNegativeWrapper->Get_dEdx_TOF();

		Int_t NumUnusedTracks = dComboWrapper->Get_NumUnusedTracks();
		double Energy_UnusedShowers = dComboWrapper->Get_Energy_UnusedShowers();

		/******************************************** EXECUTE ANALYSIS ACTIONS *******************************************/

		// Loop through the analysis actions, executing them in order for the active particle combo
		dAnalyzeCutActions->Perform_Action(); // Must be executed before Execute_Actions()
		if(!Execute_Actions()) //if the active combo fails a cut, IsComboCutFlag automatically set
			continue;

		//if you manually execute any actions, and it fails a cut, be sure to call:
			//dComboWrapper->Set_IsComboCut(true);



				/********************* CUTS **********************/
		cout << "DEBUG 7: Starting cuts section" << endl;

		// -------------------------------------------------------------
		// 0.) Preselection cut (optional, normally off):
		// -------------------------------------------------------------
		// Step 0: Measured E/p > 0.4 for both tracks
		if(ApplyPreselectionEoPCut) {
			if(FCAL_EoverPmeas_plus < 0.4 || FCAL_EoverPmeas_minus < 0.4) {
				dComboWrapper->Set_IsComboCut(true);
				continue;
			}
		}
		cout << "DEBUG 8: Passed preselection E/p cut" << endl;

		// -------------------------------------------------------------
		// 1.) Two track exclusive (not counting missing proton) cuts:
		// -------------------------------------------------------------
		// (a) No unused charged tracks
		if(ApplyNoUnusedTracksCut) {
			if(NumUnusedTracks > dMaxNumUnusedTracks) {
				dComboWrapper->Set_IsComboCut(true);
				continue;
			}
		}
		// (b) No unused neutral shower energy
		if(ApplyNoUnusedShowersCut) {
			if(Energy_UnusedShowers > 0.0){
				dComboWrapper->Set_IsComboCut(true);
				continue;
			}
		}

		// -----------------------------------------------------
		// 2.) Ensure tracks go into only the forward detectors:
		// -----------------------------------------------------
		// (a) FCAL E_1 and FCAL E_2 are greater than zero
		if(ApplyFCALEnergyNonZeroCut) {
			if (FCAL_Energy_plus <= 0.0 || FCAL_Energy_minus <= 0.0){
				dComboWrapper->Set_IsComboCut(true);
				continue;
			}
		}
		// (b) TOF dE/dx > 0 for both tracks
		if(ApplyTOFdEdxNonZeroCut) {
			if(TOF_dEdx_plus <= 0.0 || TOF_dEdx_minus <= 0.0){
				dComboWrapper->Set_IsComboCut(true);
				continue;
			}
		}

		// ------------------------------------------------------------------
		// 3.) Fiducial Cuts:
		// Turn these off if preparing a file to do systematics studies with.
		// ------------------------------------------------------------------

		// (a) Coherent peak cut (8.2 GeV < E_gamma < 8.8 GeV)
		if(ApplyBeamEnergyCut) {
			if( locBeamP4.E() < dMinBeamEnergy || locBeamP4.E() > dMaxBeamEnergy){
				dComboWrapper->Set_IsComboCut(true);
				continue;
			}
		}
		// (b) Minimum theta of tracks cut
		if(ApplyMinThetaCut) {
			if ( locPositiveP4.Theta()*TMath::RadToDeg() < 1.5 || locNegativeP4.Theta()*TMath::RadToDeg() < 1.5){
				dComboWrapper->Set_IsComboCut(true);
				continue;
			}
		}
		// (c) Maximum theta of tracks cut
		// DO NOT APPLY maximum theta of tracks cut for Rho0 file. It kills your yield.
		if(ApplyMaxThetaCut) {
			if ( locPositiveP4.Theta()*TMath::RadToDeg() > 7.5 || locNegativeP4.Theta()*TMath::RadToDeg() > 7.5){
				dComboWrapper->Set_IsComboCut(true);
				continue;
			}
		}
		// (d) Simulation/data acceptance match curve: theta_2 > 14 * (theta_1/(theta_1^(5/2) + 0.2)) deg
		// DO NOT APPLY this cut for Rho0 file--again, will kill your yield.
		if(Apply2DThetaAcceptanceMatchCut) {
			// Implementation would go here if needed
			// Currently just a placeholder flag
		}

		// (e) Valid range of momentum dependent non-linear FCAL corrections in the simulation, 0.45 GeV < |p| < 7.92 GeV
		if(ApplyMomentumRangeCut) {
			if ( locPositiveP4.Vect().Mag() < 0.45 || locPositiveP4.Vect().Mag() > 7.92){
				dComboWrapper->Set_IsComboCut(true);
				continue;
			}
			if ( locNegativeP4.Vect().Mag() < 0.45 || locNegativeP4.Vect().Mag()> 7.92){
				dComboWrapper->Set_IsComboCut(true);
				continue;
			}
		}
		// (f) Invariant mass cut (700 MeV < M_ππ < 770 MeV, to select rho0)
		if(ApplyInvariantMassCut) {
			if(M2TrackKin < 0.700 || M2TrackKin > 0.770){
				dComboWrapper->Set_IsComboCut(true);
				continue;
			}
		}
		// (g) Vertex position: 52 cm < Z < 78 cm
		if(ApplyVertexZCut) {
			if(locBeam_X4_Measured.Z() < 52 || locBeam_X4_Measured.Z() > 78){
				dComboWrapper->Set_IsComboCut(true);
				continue;
			}
		}
		// (h) Minimum chi squared confidence level: 1E-6
		if(ApplyKinFitCLCut) {
			if(dComboWrapper->Get_ConfidenceLevel_KinFit("") <= dMinKinFitCL) {
				dComboWrapper->Set_IsComboCut(true);
				continue;
			}
		}
		cout << "DEBUG 9: Passed all cuts" << endl;

		/****************************************** SECTION 4: TEST/DEVELOPMENT CUTS (NOT IN OFFICIAL ANALYSIS) **********************************************/
		// These cuts are provided for user convenience to quickly test various particle ID criteria.
		// They are typically disabled (set to false in Init()) and not part of the standard analysis.
		// Easy to copy/paste this entire section into other DSelectors for development work.

		// (a) Minimum E/P cut: FCAL energy divided by kinematic fit momentum
		if(ApplyMinEoverPCut) {
			// Use already-calculated E/P values
			if(FCAL_EoverPkin_plus < dMinEoverP || FCAL_EoverPkin_minus < dMinEoverP) {
				dComboWrapper->Set_IsComboCut(true);
				continue;
			}
		}

		// (b) Maximum E/P cut: FCAL energy divided by kinematic fit momentum
		if(ApplyMaxEoverPCut) {
			// Use already-calculated E/P values
			if(FCAL_EoverPkin_plus > dMaxEoverP || FCAL_EoverPkin_minus > dMaxEoverP) {
				dComboWrapper->Set_IsComboCut(true);
				continue;
			}
		}

		// (c) Minimum FCAL elasticity cut: (E_plus + E_minus) / E_beam
		if(ApplyMinFCALElasticityCut) {
			// Use already-calculated FCAL_Elasticity
			if(FCAL_Elasticity < dMinFCALElasticity) {
				dComboWrapper->Set_IsComboCut(true);
				continue;
			}
		}

		// (d) Maximum FCAL elasticity cut: (E_plus + E_minus) / E_beam
		if(ApplyMaxFCALElasticityCut) {
			// Use already-calculated FCAL_Elasticity
			if(FCAL_Elasticity > dMaxFCALElasticity) {
				dComboWrapper->Set_IsComboCut(true);
				continue;
			}
		}

		// (e) Maximum FCAL DOCA cut: distance of closest approach to shower center
		if(ApplyMaxFCALDOCACut) {
			// Use already-calculated DOCA values
			if(TrackFCAL_DOCA_plus > dMaxFCALDOCA || TrackFCAL_DOCA_minus > dMaxFCALDOCA) {
				dComboWrapper->Set_IsComboCut(true);
				continue;
			}
		}

		// (f) Maximum TOF dE/dx cut
		if(ApplyMaxTOFdEdxCut) {
			// Use already-calculated TOF dE/dx values
			if(TOF_dEdx_plus > dMaxTOFdEdx || TOF_dEdx_minus > dMaxTOFdEdx) {
				dComboWrapper->Set_IsComboCut(true);
				continue;
			}
		}


		// Declare variables for MLP/BDT responses (will be filled conditionally)
        Float_t mlp_response_minus = -999.0;  // Use sentinel value as default
        Float_t mlp_response_plus = -999.0;
		Float_t bdt_response_minus = -999.0;
		Float_t bdt_response_plus = -999.0;

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




		/**************************************** EXAMPLE: FILL CUSTOM OUTPUT BRANCHES **************************************/

		/*
		TLorentzVector locMyComboP4(8.0, 7.0, 6.0, 5.0);
		//for arrays below: 2nd argument is value, 3rd is array index
		//NOTE: By filling here, AFTER the cuts above, some indices won't be updated (and will be whatever they were from the last event)
			//So, when you draw the branch, be sure to cut on "IsComboCut" to avoid these.
		dTreeInterface->Fill_Fundamental<Float_t>("my_combo_array", -2*loc_i, loc_i);
		dTreeInterface->Fill_TObject<TLorentzVector>("my_p4_array", locMyComboP4, loc_i);
		*/

		/**************************************** EXAMPLE: BEST chi2 METHOD *****************************************/

		if(ThisComboIsBestChiSq) {
			dHist_BeamEnergy_BestChiSq->Fill(locBeamP4.E());

			dHist_InvMass_BestChiSq_total->Fill(M2TrackKin);
			dHist_MLP1_vs_MLP2_BestChiSq_total->Fill(mlp_response_minus, mlp_response_plus);
			dHist_BDT1_vs_BDT2_BestChiSq_total->Fill(bdt_response_minus, bdt_response_plus);


			if(FCAL_EoverPmeas_plus > 0.4 && FCAL_EoverPmeas_minus > 0.4) {
				dHist_InvMass_BestChiSq_EoPprecut->Fill(M2TrackKin);
				dHist_MLP1_vs_MLP2_BestChiSq_EoPprecut->Fill(mlp_response_minus, mlp_response_plus);
				dHist_BDT1_vs_BDT2_BestChiSq_EoPprecut->Fill(bdt_response_minus, bdt_response_plus);
				// Add more histograms as

			}
		}

        // Need to uncomment the section computing combo timing info before running this block of code
		//if(locUsedSoFar_Event == false)
		//{
			// Fill the histogram only when the beam bunch is in-time. 
			//if(!locRelBeamBucket)
			//{
			//	dHist_BeamEnergy_BestChiSq->Fill(locBeamP4.E());
			//	locUsedSoFar_Event = true;
			//}
		//}

		/**************************************** EXAMPLE: HISTOGRAM BEAM ENERGY *****************************************/

		//Histogram beam energy (if haven't already)
		if(locUsedSoFar_BeamEnergy.find(locBeamID) == locUsedSoFar_BeamEnergy.end())
		{
			dHist_BeamEnergy->Fill(locBeamP4.E()); // Fills in-time and out-of-time beam photon combos
			//dHist_BeamEnergy->Fill(locBeamP4.E(),locHistAccidWeightFactor); // Alternate version with accidental subtraction

			locUsedSoFar_BeamEnergy.insert(locBeamID);
		}

		/************************************ EXAMPLE: HISTOGRAM MISSING MASS SQUARED ************************************/

		//Missing Mass Squared
		double locMissingMassSquared = locMissingP4_Measured.M2();

		//Uniqueness tracking: Build the map of particles used for the missing mass
			//For beam: Don't want to group with final-state photons. Instead use "UnknownParticle" PID (not ideal, but it's easy).
		map<Particle_t, set<Int_t> > locUsedThisCombo_MissingMass;
		locUsedThisCombo_MissingMass[UnknownParticle].insert(locBeamID); //beam
		locUsedThisCombo_MissingMass[PiPlus].insert(locPiPlusTrackID);
		locUsedThisCombo_MissingMass[PiMinus].insert(locPiMinusTrackID);


		//compare to what's been used so far
		if(locUsedSoFar_MissingMass.find(locUsedThisCombo_MissingMass) == locUsedSoFar_MissingMass.end())
		{
			//unique missing mass combo: histogram it, and register this combo of particles
			dHist_MissingMassSquared->Fill(locMissingMassSquared); // Fills in-time and out-of-time beam photon combos
			//dHist_MissingMassSquared->Fill(locMissingMassSquared,locHistAccidWeightFactor); // Alternate version with accidental subtraction

			locUsedSoFar_MissingMass.insert(locUsedThisCombo_MissingMass);
		}

		//E.g. Cut
		//if((locMissingMassSquared < -0.04) || (locMissingMassSquared > 0.04))
		//{
		//	dComboWrapper->Set_IsComboCut(true);
		//	continue;
		//}

		/****************************************** FILL FLAT TREE (IF DESIRED) ******************************************/

		// RECOMMENDED: FILL ACCIDENTAL WEIGHT
		// dFlatTreeInterface->Fill_Fundamental<Double_t>("accidweight",locHistAccidWeightFactor);

		/*
		//FILL ANY CUSTOM BRANCHES FIRST!!
		Int_t locMyInt_Flat = 7;
		dFlatTreeInterface->Fill_Fundamental<Int_t>("flat_my_int", locMyInt_Flat);

		TLorentzVector locMyP4_Flat(4.0, 3.0, 2.0, 1.0);
		dFlatTreeInterface->Fill_TObject<TLorentzVector>("flat_my_p4", locMyP4_Flat);

		for(int loc_j = 0; loc_j < locMyInt_Flat; ++loc_j)
		{
			dFlatTreeInterface->Fill_Fundamental<Int_t>("flat_my_int_array", 3*loc_j, loc_j); //2nd argument = value, 3rd = array index
			TLorentzVector locMyComboP4_Flat(8.0, 7.0, 6.0, 5.0);
			dFlatTreeInterface->Fill_TObject<TLorentzVector>("flat_my_p4_array", locMyComboP4_Flat, loc_j);
		}
		*/

		//FILL FLAT TREE
		//Fill_FlatTree(); //for the active combo
	} // end of combo loop

	//FILL HISTOGRAMS: Num combos / events surviving actions
	Fill_NumCombosSurvivedHists();

	/******************************************* LOOP OVER THROWN DATA (OPTIONAL) ***************************************/
/*
	//Thrown beam: just use directly
	if(dThrownBeam != NULL)
		double locEnergy = dThrownBeam->Get_P4().E();

	//Loop over throwns
	for(UInt_t loc_i = 0; loc_i < Get_NumThrown(); ++loc_i)
	{
		//Set branch array indices corresponding to this particle
		dThrownWrapper->Set_ArrayIndex(loc_i);

		//Do stuff with the wrapper here ...
	}
*/
	/****************************************** LOOP OVER OTHER ARRAYS (OPTIONAL) ***************************************/
/*
	//Loop over beam particles (note, only those appearing in combos are present)
	for(UInt_t loc_i = 0; loc_i < Get_NumBeam(); ++loc_i)
	{
		//Set branch array indices corresponding to this particle
		dBeamWrapper->Set_ArrayIndex(loc_i);

		//Do stuff with the wrapper here ...
	}

	//Loop over charged track hypotheses
	for(UInt_t loc_i = 0; loc_i < Get_NumChargedHypos(); ++loc_i)
	{
		//Set branch array indices corresponding to this particle
		dChargedHypoWrapper->Set_ArrayIndex(loc_i);

		//Do stuff with the wrapper here ...
	}

	//Loop over neutral particle hypotheses
	for(UInt_t loc_i = 0; loc_i < Get_NumNeutralHypos(); ++loc_i)
	{
		//Set branch array indices corresponding to this particle
		dNeutralHypoWrapper->Set_ArrayIndex(loc_i);

		//Do stuff with the wrapper here ...
	}
*/

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
	if(!locIsEventCut && dOutputTreeFileName != "")
		Fill_OutputTree();
*/

	return kTRUE;
}

void DSelector_pippimmissp_psubMatrix::Finalize(void)
{
	//Save anything to output here that you do not want to be in the default DSelector output ROOT file.

	//Otherwise, don't do anything else (especially if you are using PROOF).
		//If you are using PROOF, this function is called on each thread,
		//so anything you do will not have the combined information from the various threads.
		//Besides, it is best-practice to do post-processing (e.g. fitting) separately, in case there is a problem.

	//DO YOUR STUFF HERE

	//CALL THIS LAST
	DSelector::Finalize(); //Saves results to the output file
}
