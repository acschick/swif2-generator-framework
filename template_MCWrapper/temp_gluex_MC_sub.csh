#!/bin/tcsh
setenv APPTAINER_BIND "/work/halld,/work/osgpool,/volatile/halld/,$PWD,/cvmfs/singularity.opensciencegrid.org/"
$MCWRAPPER_CENTRAL/gluex_MC.py RBHGPATHTOMCCONFIG RBHGRUNRANGE RBHGNEVENTS per_file=RBHGPERFILE batch=RBHGBATCHMODE logdir=RBHGLOGDIR
 