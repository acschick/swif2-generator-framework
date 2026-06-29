
void call_DSelector (TString file)
{
// issue the tree->Process, so that it can be run from the command line
//
cout << "call_eeDSelector2: file=" << file << endl;
gROOT->LoadMacro("$ROOT_ANALYSIS_HOME/scripts/Load_DSelector.C");
//pippimmissprot__B4_Tree->Process("DSelector_2piMissingProton.C+"); //This is Rory's definitely bad.
//epemmissprot__U1_Tree->Process("DSelector_2eMissingProton.C+"); //This is Rory's potentially good

//epemmissprot__B2_U1_Tree->Process("DSelector_2eMissingProton.C+");

//pippimmissprot__B2_Tree->Process("DSelector_pippimmissprot.C+"); //for version 27 I found on cache.
//epemmissprot__B2_Tree->Process("DSelector_epemmisspTMVA.C+");
 epemmissprot__B2_Tree->Process("DSelector_2eMissingProton_Systematics.C+");
}
