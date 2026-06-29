#!/bin/csh -f                                                                               
set echo
                                                                                                     

root -b -q /volatile/halld/home/acschick/RunPeriod-2018-01/Analysis/gfdi2/tree_pippimmissprot__B2/042011/tree_pippimmissprot__B2_042011_000.root 'call_DSelector.C("DSelector_pippimmissp_psubMatrix.C+")' >! DSelector_cout.list


unset echo

    
                                                                                            
