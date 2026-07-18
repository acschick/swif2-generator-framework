#!/bin/csh -f                                                                                                                                                                                   
set echo
#                                                                                                                                                                                               
# streamline.csh                                                                                                                                                                                
# Andrew Schick. Bastardized from Elton Smith's streamline. June 29, 2018.                                                                                                                                     
# Streamline instructions to call a DSelector and to make plots of the results.

#root -b -q tree_epemmissprot__B2_041490_004.root 'call_DSelector.C("DSelector_2eMissingProton_Systematics.C+")' >! DSelector_cout.list
#root -b -q /volatile/halld/home/acschick/RBHG/FFS1/SIM/1801_0DEG_FF1_DBLRAD/Simulation/root/trees/tree_epemmissprot__B2_042265.root 'call_DSelector.C("DSelector_2eMissingProton_Systematics.C+")' >! DSelector_cout.list

#root -b -q /cache/halld/RunPeriod-2022-05/skims/epem/analysis_epemlskim/135DEG/tree_gpb208_epemmisspb208__B2/101549/tree_gpb208_epemmisspb208__B2_101549_000.root 'call_DSelector.C("DSelector_2eMissingProton_Systematics.C+")' >! DSelector_cout.list

#root -b -q /volatile/halld/home/acschick/RBHG/FFS1/SIM/mergedTest2.root 'call_DSelector.C("DSelector_2eMissingProton_Systematics.C+")' >! DSelector_cout.list
set i = 0

foreach file ( \
    /work/halld/home/acschick/channels/swif2-generator-framework/tree_1801_0DEG_FFN_DBLRAD.root \
    /work/halld/home/acschick/channels/swif2-generator-framework/tree_1801_135DEG_FFN_DBLRAD.root \
    /work/halld/home/acschick/channels/swif2-generator-framework/tree_1801_45DEG_FFN_DBLRAD.root \
    /work/halld/home/acschick/channels/swif2-generator-framework/tree_1801_90DEG_FFN_DBLRAD.root \
    /work/halld/home/acschick/channels/swif2-generator-framework/tree_1801_AMO_FFN_DBLRAD.root \
    /work/halld/home/acschick/channels/swif2-generator-framework/tree_1808_0DEG_FFN_DBLRAD.root \
    /work/halld/home/acschick/channels/swif2-generator-framework/tree_1808_135DEG_FFN_DBLRAD.root \
    /work/halld/home/acschick/channels/swif2-generator-framework/tree_1808_45DEG_FFN_DBLRAD.root \
    /work/halld/home/acschick/channels/swif2-generator-framework/tree_1808_90DEG_FFN_DBLRAD.root \
    /work/halld/home/acschick/channels/swif2-generator-framework/tree_1808_AMO_FFN_DBLRAD.root \
    )
    echo "Running iteration $i on file $file"
    root -b -q $file 'call_DSelector.C("DSelector_2eMissingProton_Systematics.C+")' >! DSelector_${i}.list
    @ i++
end
#foreach file ( \
#    /work/halld/home/acschick/channels/cpp_epemmisspb208/tree_gpb208_epemmisspb208__B2_135DEG_FFN_IRADOFF.root \
#    /work/halld/home/acschick/channels/cpp_epemmisspb208/tree_gpb208_epemmisspb208__B2_45DEG_FFN_IRADOFF.root \
#    )
#    echo "Running iteration $i on file $file"
#    root -b -q $file 'call_DSelector.C("DSelector_2eMissingProton_Systematics.C+")' >! DSelector_CPP_${i}.list
#    @ i++
#end

unset echo

    
                                                                                            
