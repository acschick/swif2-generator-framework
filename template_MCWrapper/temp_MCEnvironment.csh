setenv HOME /home/RBHGCUEUSERNAME
setenv CCDB_CONNECTION sqlite:////work/halld/ccdb_sqlite/21/ccdb.sqlite
setenv JANA_CALIB_URL $CCDB_CONNECTION
source /group/halld/Software/build_scripts/gluex_env_jlab.csh RBHGPATHTOXML
RBHGCPPGEOMETRY
RBHGCPPMCVARIATION
