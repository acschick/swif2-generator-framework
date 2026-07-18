
void call_DSelector (TString file)
{
// issue the tree->Process, so that it can be run from the command line
//
cout << "call_eeDSelector2: file=" << file << endl;
gROOT->LoadMacro("$ROOT_ANALYSIS_HOME/scripts/Load_DSelector.C");

//epemmissprot__B2_U1_Tree->Process("DSelector_2eMissingProton.C+");
 epemmissprot__B2_Tree->Process("DSelector_2eMissingProton_Systematics.C+");
// gpb208_epemmisspb208__B2_Tree->Process("DSelector_2eMissingProton_Systematics.C+");
}
