
void call_DSelector (TString file)
{
// issue the tree->Process, so that it can be run from the command line
//
cout << "call_eeDSelector2: file=" << file << endl;
gROOT->LoadMacro("$ROOT_ANALYSIS_HOME/scripts/Load_DSelector.C");
pippimmissprot__B2_Tree->Process("DSelector_pippimmissp_psubMatrix.C+"); 
}
