import subprocess
import re
import struct
import os
import sys
import math
import difflib
import json
import hashlib
from datetime import datetime

# Parse command line arguments for interactive mode
INTERACTIVE_MODE = '--interactive' in sys.argv or '--local' in sys.argv
if INTERACTIVE_MODE:
    print("="*60)
    print("INTERACTIVE MODE: Jobs will run locally, not submitted to swif2")
    print("="*60)
    sys.argv = [arg for arg in sys.argv if arg not in ['--interactive', '--local']]  # Remove flag from argv

# Parse command line arguments for config file
CONFIG_FILE = "RBHG.config"  # Default
for i, arg in enumerate(sys.argv):
    if arg.startswith('--config='):
        CONFIG_FILE = arg.split('=', 1)[1]
        sys.argv.pop(i)
        break
    elif arg == '--config' and i + 1 < len(sys.argv):
        CONFIG_FILE = sys.argv[i + 1]
        sys.argv.pop(i)
        sys.argv.pop(i)  # Remove the next arg too
        break


#############################################################################
####################### Configuration Loading ##############################
#############################################################################

def load_config(config_file="RBHG.config"):
    """
    Load configuration from simple key-value format file.
    Format: KEY VALUE  (separated by whitespace)
    Lines starting with # are comments and ignored
    """
    config = {}
    
    # If config_file is an absolute path, use it directly
    # Otherwise, look for it relative to this script's directory
    if os.path.isabs(config_file):
        config_path = config_file
    else:
        config_path = os.path.join(os.path.dirname(__file__), config_file)
    
    if not os.path.exists(config_path):
        print(f"Warning: Config file {config_path} not found. Using defaults.")
        return {}
    
    with open(config_path, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            # Skip comments and empty lines
            if not line or line.startswith('#'):
                continue
            
            # Strip inline comments (everything after #)
            if '#' in line:
                line = line.split('#')[0].strip()
            
            parts = line.split()
            if len(parts) < 2:
                # Allow single-word lines with no value (treat as empty string)
                if len(parts) == 1:
                    config[parts[0]] = ""
                    continue
                print(f"Warning: Invalid config line {line_num}: {line}")
                continue
            
            key = parts[0]
            value = parts[1]
            
            # Convert boolean strings
            if value.lower() in ('true', 'false'):
                config[key] = value.lower() == 'true'
            # Convert numeric strings
            elif value.replace('.', '').replace('-', '').isdigit():
                if '.' in value:
                    config[key] = float(value)
                else:
                    config[key] = int(value)
            else:
                config[key] = value
    
    # Second pass: resolve variable substitutions like {USERNAME}
    for key, value in config.items():
        if isinstance(value, str) and '{' in value and '}' in value:
            # Simple variable substitution
            resolved_value = value
            for var_key, var_value in config.items():
                if isinstance(var_value, str):
                    resolved_value = resolved_value.replace(f'{{{var_key}}}', str(var_value))
            config[key] = resolved_value
    
    print(f"Loaded {len(config)} configuration parameters from {config_file}")
    return config

# Load configuration
config = load_config(CONFIG_FILE)

def get_config(key, default=None):
    """Get a configuration value with optional default."""
    return config.get(key, default)


def replace_path_placeholders(path_string, generator_type, username, framework_home):
    """
    Replace placeholders in path strings with actual values.
    
    Args:
        path_string: Path string potentially containing {GENERATOR_TYPE}, {USERNAME}, {FRAMEWORK_HOME}
        generator_type: Generator type (e.g., 'RBHG', 'SPIZG')
        username: Username for paths
        framework_home: Framework home directory path
        
    Returns:
        Path string with placeholders replaced
    """
    if not isinstance(path_string, str):
        return path_string
    
    replacements = {
        '{GENERATOR_TYPE}': generator_type,
        '{USERNAME}': username,
        '{FRAMEWORK_HOME}': framework_home
    }
    
    result = path_string
    for placeholder, value in replacements.items():
        result = result.replace(placeholder, value)
    
    return result


#############################################################################
############################ swif2_RBHG.py ##################################
#############################################################################

# Framework home directory (root of swif2-generator-framework)
# In the new framework structure, this should point to the swif2-generator-framework root
# Default: three directories up from this script (generators/RBHG/swif2_RBHG.py -> framework root)
FrameworkHomeDirectory = get_config("FRAMEWORK_HOME", os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

# Generator type (for directory paths)
GENERATOR_TYPE = get_config("GENERATOR_TYPE", "RBHG")

# Files in {FrameworkHomeDirectory}/generators/RBHG/template_RBHG/
###################################################
## a.) The template file version that will be edited and compiled into an .exe:
temp_fortran_file = get_config("FORTRAN_TEMPLATE", "temp_lepton_event_v113.f")
## b.) The bash script that constitutes the swif job:
if INTERACTIVE_MODE:
    RBHG_script = get_config("JOB_SCRIPT_INTERACTIVE", "swif2_RBHG_script_interactive.sh")
else:
    RBHG_script = get_config("JOB_SCRIPT", "swif2_RBHG_script.sh")
## c.) The python script that converts output vectors to hddm file (centralized)
ascii2hddm_script = os.path.join(FrameworkHomeDirectory, "hddm_scripts", "ascii2hddm.py")

# Environment file used by farm job:
ENVFILE = get_config("ENV_FILE", "/group/halld/www/halldweb/html/halld_versions/version_6.5.0.xml")
# can leave the above empty "", ver xml not particularly important for running a fortran exe
# but *it is* now that we are using hddm_s to make an hddm file.

# JLAB CUE username:
CUEusername = get_config("USERNAME", "acschick")

presimStudy = get_config("PRESIM_STUDY", False)                     # theta min changes from 0.75 deg to 1.5 deg
coherentPeakStudy = get_config("COHERENT_PEAK_STUDY", False)                # Sets E0 inside coherent peak range (8.2, 8.8) #edit this for CPP

# Event count matching
FIXED_EVENT_COUNT = get_config("FIXED_EVENT_COUNT", True)                # Override automatic event counts to use specific fixed values
TENX_FLAG = get_config("TENX_FLAG", False)                       # 10x the statistics of real data
RUN_PERIOD = str(get_config("RUN_PERIOD", ""))                     # valid: 1801, 1808, FULL2018, 2205, or ""
                                        # Comment out RUN_PERIOD line in config for empty preset
                                        # Add '_AMO' suffix to any preset to include AMO data
                                        # Add '_FFS' suffix for Form Factor Study (generates 4 datasets)
                                        # Suffixes can be combined: e.g. FULL2018_AMO_FFS
                                        # Examples: FULL2018_AMO, 1801_FFS, FULL2018_AMO_FFS
                                        # '_AMO': includes AMO data (0,90,135,45,AMO) vs polarized only (0,90,135,45)
                                        # '_FFS': generates qDATAq, BCFFN, BCFF1, SIM datasets for form factor analysis
                                        # automatically sets experiment to either GlueX/CPP

# Store original RUN_PERIOD state for later
RUN_PERIOD_SPECIFIED = bool(RUN_PERIOD.strip())

# Apply sensible defaults when RUN_PERIOD is blank (for cobrems lookup only, not preset behavior)
if not RUN_PERIOD_SPECIFIED:
    TENX_FLAG = False
    # Set default reference run period for cobrems/event lookup
    Experiment = get_config("EXPERIMENT", "GlueX")
    PolDeg = get_config("POL_DEG", "135")
    
    if Experiment.upper() in ["GLUEX", "GX"]:
        DEFAULT_RUN_PERIOD = "1808"  # Default to 2018-08 for GlueX (most common)
        print(f"No RUN_PERIOD specified - using {DEFAULT_RUN_PERIOD} as reference for cobrems/event lookup")
        print(f"  Will generate single job for POL_DEG={PolDeg} (not full {DEFAULT_RUN_PERIOD} preset)")
        print(f"  (GlueX has: 0, 45, 90, 135, AMO)")
    elif Experiment.upper() == "CPP":
        DEFAULT_RUN_PERIOD = "2205"  # Default to 2022-05 for CPP
        print(f"No RUN_PERIOD specified - using {DEFAULT_RUN_PERIOD} as reference for cobrems/event lookup")
        print(f"  Will generate single job for POL_DEG={PolDeg} (not full {DEFAULT_RUN_PERIOD} preset)")
        print(f"  (CPP has: 45, 135, AMO)")
        # Validate CPP polarization
        if str(PolDeg).upper() not in ["45", "135", "AMO"]:
            print(f"  WARNING: CPP only supports POL_DEG = 45, 135, or AMO")
            print(f"           Your setting POL_DEG={PolDeg} may not have cobrems/event data")
    else:
        print(f"Warning: Unknown experiment '{Experiment}', using 1808 as reference")
        DEFAULT_RUN_PERIOD = "1808"
else:
    DEFAULT_RUN_PERIOD = None  # Not needed when RUN_PERIOD is explicitly set

# Event count configuration
if FIXED_EVENT_COUNT:
    # Fixed event count mode: use these values regardless of RUN_PERIOD
    nevents_total = str(get_config("NEVENTS_TOTAL", 10000))
    nevents_perfile = str(get_config("NEVENTS_PERFILE", 5000))
else:
    # Auto event count mode: load from RunPeriods.json (requires RUN_PERIOD to be set)
    # Check if user put placeholder value indicating they know preset will be used
    config_nevents = get_config("NEVENTS_TOTAL", 10000)
    if isinstance(config_nevents, str) and config_nevents.upper() in ["USEPRESET", "USE_PRESET", "PRESET", "AUTO"]:
        # User explicitly indicated preset mode - use placeholder
        nevents_total = "0"  # Will be overwritten by RunPeriods.json
        print("NEVENTS_TOTAL set to use preset values from RunPeriods.json")
    else:
        # Use as fallback if RunPeriods.json lookup fails
        nevents_total = str(config_nevents)
    nevents_perfile = str(get_config("NEVENTS_PERFILE", 5000))



    
####################################
######## Generator Control #########
####################################
"""
Directory Structure:
  Regular studies: {outputTopDir}/{study_name}/{nametag}_{form_factor}/{run_period}_{polarization}_{extras}/
    Example: TS1/MyTest_FFN/1801_0DEG_DBLRAD/
  
  Form Factor Studies (_FFS): {outputTopDir}/{study_name}/{nametag}_{form_factor}/{run_period}_{polarization}_{form_factor}_{extras}/
    Example: FFS3/BC_FFN/1801_0DEG_FFN_DBLRAD/
             FFS3/BC_FF1/1801_0DEG_FF1_DBLRAD/
             FFS3/SIM_FF1/1801_0DEG_FF1_DBLRAD/
             FFS3/qDATAq_FFN/1801_0DEG_FFN_DBLRAD/
    
  Note: Consistent underscore naming (nametag_formfactor). FFS includes form factor at both levels.
inside t study number "x"-> TS(x)
Assuming full run period matching, 8 directories will be created for the two 2018 run periods,
and 4 polarization orientations. 
"""                                                                                                                                                     
#### Generator Control ####
study_name = get_config("STUDY_NAME", "JSONTESTING")              # Top level directory to contain modified configuration outputs
nametag = get_config("NAMETAG", "N2")          # e.g. BC (bin correction), SIM (simulation), qDATAq (psuedo-data)
nevents_perfile = str(get_config("NEVENTS_PERFILE", 5000))              

form_factor = get_config("FORM_FACTOR", "FFN")                     # Nuclear charge form factor. "FFN" for FF=dipole,  or "FF1" to set to 1.
lepton = get_config("LEPTON", "mu")                           # "ee" for electron, or "mu" for muon
# we currently have the cross section get overwritten to Heitler and CobremsDistribution to False if PolDeg = AMO.
# Should we set it up so that if PolDeg = 0,45,90,135 then xsctn = Berlin and Cobrems = True? Ask Rory.
BH_xsctn_formulation = get_config("BH_XSCTN_FORMULATION", "Berlin")         # set to "Berlin" (polarized) or "Heitler" (unpolarized)
CobremsDistribution = get_config("COBREMS_DISTRIBUTION", True)              # Use experimental Coherent Bremsstrahlung distribution (True) or 1/E distribution (False)
CobremsVarBin = get_config("COBREMS_VARBIN", False)                        # Use variable-width tagger binning for bremsstrahlung distribution
CobremsFileGlueX = get_config("COBREMS_FILE_GLUEX", "")                   # Override default GlueX cobrems file (blank = auto-select)
CobremsFileCPP = get_config("COBREMS_FILE_CPP", "")                       # Override default CPP cobrems file (blank = auto-select)
internal_radiation  = get_config("INTERNAL_RADIATION", True)              # Mo & Tsai, internally radiated brem. photon
single_radiation    = get_config("SINGLE_RADIATION", False)             # v113 feature: if True, only one lepton radiates; if False, both leptons radiate
hypgeom_radiation   = get_config("HYPGEOM_RADIATION", False)             # Future: hypergeometric radiation for muons (not implemented yet)
# automatically overwritten if RUN_PERIOD != "" 
PolDeg = str(get_config("POL_DEG", 45))                           # Set photon beam polarization orientation to '0', '90', '135', '45', or 'AMO'
Experiment = get_config("EXPERIMENT", "GlueX")                      # "GlueX" or "CPP"

#### MCWrapper / Simulation Control ####
MCWrapper_ReconEnv = get_config("MCWRAPPER_RECON_ENV_OVERRIDE", "")        # Override recon environment (blank = use RunPeriods.json default)
MCWrapper_AnalysisEnv = get_config("MCWRAPPER_ANALYSIS_ENV_OVERRIDE", "")  # Override analysis environment (blank = use RunPeriods.json default)
MCWrapper_BackgroundVersion = get_config("MCWRAPPER_BACKGROUND_VERSION", "") # Override background version (blank = use RunPeriods.json default)
MCWrapper_JanaConfig = get_config("MCWRAPPER_JANA_CONFIG_OVERRIDE", "")    # Override JANA config (blank = use RunPeriods.json default)
MCWrapper_RCDBQuery = get_config("MCWRAPPER_RCDB_QUERY", "")               # Override RCDB query (blank = auto-generate from RunPeriods.json)
MCWrapper_UseCharacteristicRun = get_config("MCWRAPPER_USE_CHARACTERISTIC_RUN", False) # Use characteristic run from RunPeriods.json

#### Histogram Control ####
hist_invariantMass    = get_config("HIST_INVARIANT_MASS", True)            # Histogram output control for invariant mass W
hist_energyFraction_x = get_config("HIST_ENERGY_FRACTION_X", True)            # Histogram output control for energy fraction x
hist_mandelstam_t     = get_config("HIST_MANDELSTAM_T", True)            # Histogram output control for Mandelstam t
hist_t_varBinWidth    = get_config("HIST_T_VAR_BIN_WIDTH", True)            # Same as above, but with GlueX resolution matched bin widths
hist_elasticity       = get_config("HIST_ELASTICITY", True)            # Histogram output control for elasticity
hist_missingmass      = get_config("HIST_MISSING_MASS", True)
hist_theta            = get_config("HIST_THETA", True)            # Histogram output control for lab theta
hist_phi_of_JT        = get_config("HIST_PHI_OF_JT", True)            # Histogram output control for polarization observable, phi of the transverse vector J
hist_Egamma           = get_config("HIST_EGAMMA", True)            # Histogram output control for the incident photon energy
####  Output Control  ####
ascii_vector_output   = get_config("ASCII_VECTOR_OUTPUT", True)            # Control for writing vector event file. Used as input to the GlueX simulation.
integral_xsection     = get_config("INTEGRAL_XSECTION", False)           # Control for outputting integrated cross sections 
verbose_output        = get_config("VERBOSE_OUTPUT", False)              # Control for test outputs to screen (max cross section convergence)

# Define GlueX and CPP Run Conditions                                                                                                                        
experiment_parameters = {
    "RBHG_ELO_GLUEX": str(get_config("GLUEX_PHOTON_ENERGY_MIN", 7.0)),            # GlueX: Minimum photon energy to generate
    "RBHG_EHI_GLUEX": str(get_config("GLUEX_PHOTON_ENERGY_MAX", 11.8)),           # GlueX: Maximum photon energy to generate
    "RBHG_ZTGT_GLUEX": "1",             # GlueX: Atomic number of Hydrogen (target)
    "RBHG_E0_GLUEX": str(get_config("GLUEX_ELECTRON_BEAM_ENERGY", 11.8)),            # GlueX: Electron Beam Energy (11 GeV for both GlueX and CPP)
    "RBHG_ECOHERENT_GLUEX": str(get_config("GLUEX_COHERENT_PEAK_EDGE", 8.8)),      # GlueX: Edge of the coherent peak
    "RBHG_POL_GLUEX": str(get_config("GLUEX_PHOTON_POLARIZATION", 1.0)),            # GlueX: Fraction of polarization of photon beam
    "RBHG_THETAMIN_GLUEX": str(get_config("GLUEX_LEPTON_THETA_MIN", 0.75)),      # GlueX: Minimum lepton theta to generate (degrees)
    "RBHG_THETAMAX_GLUEX": str(get_config("GLUEX_LEPTON_THETA_MAX", 14)),        # GlueX: Maximum lepton theta to generate (degrees) 13.12 to hit FCAL   (0.1 - 5.0 for ECAL)
    "RBHG_ELO_CPP": str(get_config("CPP_PHOTON_ENERGY_MIN", 4.5)),              # CPP  : Minimum photon energy to generate, currently set to lower bound of coherent peak  
    "RBHG_EHI_CPP": str(get_config("CPP_PHOTON_ENERGY_MAX", 6.0)),              # CPP  : "
    "RBHG_ZTGT_CPP": "82",              # CPP  : "
    "RBHG_E0_CPP": str(get_config("CPP_ELECTRON_BEAM_ENERGY", 11.0)),              # CPP  : "
    "RBHG_ECOHERENT_CPP": str(get_config("CPP_COHERENT_PEAK_EDGE", 5.7)),        # CPP  : "
    "RBHG_POL_CPP": str(get_config("CPP_PHOTON_POLARIZATION", 1.0)),              # CPP  : "
    "RBHG_THETAMIN_CPP": str(get_config("CPP_LEPTON_THETA_MIN", 0.3)),         # CPP  : "
    "RBHG_THETAMAX_CPP": str(get_config("CPP_LEPTON_THETA_MAX", 7))            # CPP  : " (was 5.3)
}


if coherentPeakStudy:
    experiment_parameters["RBHG_ELO_GLUEX"] = "8.2"
    experiment_parameters["RBHG_EHI_GLUEX"] = "8.8"

if presimStudy:
    experiment_parameters["RBHG_THETAMIN_GLUEX"] = "1.5"


####################################
####### swif2 job parameters #######
####################################
swif_addjob_dict_stub = {
    "-cores"           :   "1",                     # no advantage for multicore
    "-ram"             :   get_config("RAM_USAGE", "500MB"),                  # max 2GB or scicomp will spam you for inefficiency
    "-disk"            :   get_config("DISK_USAGE", "10GB"),                   # Space needed per job 
    "-time"            :   get_config("TIME_LIMIT", "60minutes"),                   # time limit per job
    "-account"         :   get_config("PROJECT", "halld"),                 
    "-partition"       :   get_config("TRACK", "production"),
    "-os"              :   "el9",
    "RBHG_bash_script" :  f"{RBHG_script}"
}


###########################  END OF USER INPUT ############################
##########################    DO NOT TOUCH    #############################

























def Fortran_TorF(pythonLogic):
    return ".true." if pythonLogic else ".false."

def create_rbhg_config(study_name, nametag, form_factor, lepton, BH_xsctn_formulation, 
                       Experiment, PolDeg, nevents_total, nevents_perfile, internal_radiation, 
                       CobremsDistribution, RUN_PERIOD, TENX_FLAG, fortran_version,
                       histogram_settings, output_settings, experiment_parameters, directories=None,
                       single_radiation=False, hypgeom_radiation=False, run_period=None, polarization=None):
    """Create RBHG configuration dictionary for JSON export"""
    
    # Extract version from fortran filename  
    version_match = re.search(r'v\d+', fortran_version) if fortran_version else None
    generator_version = version_match.group() if version_match else "unknown"
    study_type = determine_study_type(RUN_PERIOD, study_name)
    
    # Determine target based on experiment
    target = "p" if Experiment == "GlueX" else "pb208"
    
    # Load RunPeriods.json data for MCWrapper settings if run_period and polarization available
    mcwrapper_settings = {}
    if run_period and polarization and directories:
        rp_data = load_runperiods_data(run_period, polarization, FrameworkHomeDirectory)
        
        if rp_data:
            # Use explicit overrides from RBHG.config if set, otherwise use RunPeriods.json defaults
            mcwrapper_settings = {
                "recon_env": MCWrapper_ReconEnv if MCWrapper_ReconEnv.strip() else rp_data.get("recon_env", ""),
                "analysis_env": MCWrapper_AnalysisEnv if MCWrapper_AnalysisEnv.strip() else rp_data.get("analysis_env", ""),
                "background_version": MCWrapper_BackgroundVersion if MCWrapper_BackgroundVersion.strip() else rp_data.get("background_version", ""),
                "jana_config": MCWrapper_JanaConfig if MCWrapper_JanaConfig.strip() else rp_data.get("jana_epem_config", ""),
                "batch_system": rp_data.get("batch_system", "swif2"),
                "use_characteristic_run": MCWrapper_UseCharacteristicRun,
                "characteristic_run": rp_data.get("characteristic_run", "")
            }
            
            # Handle RCDB query - substitute {polDeg} placeholder if present
            if MCWrapper_RCDBQuery.strip():
                mcwrapper_settings["rcdb_query"] = MCWrapper_RCDBQuery
            elif rp_data.get("RCDB_query"):
                # Auto-generate from template in RunPeriods.json
                rcdb_template = rp_data["RCDB_query"]
                # Extract numeric polarization angle (remove "DEG" suffix)
                pol_numeric = str(polarization).replace("DEG", "").replace("DEGREES", "")
                if pol_numeric != "AMO":
                    rcdb_query = rcdb_template.replace("{polDeg}", pol_numeric)
                else:
                    rcdb_query = rcdb_template  # For AMO, leave as-is or handle differently
                mcwrapper_settings["rcdb_query"] = rcdb_query
            else:
                mcwrapper_settings["rcdb_query"] = ""
            
            # Add run_range from RunPeriods.json
            mcwrapper_settings["run_range"] = rp_data.get("run_range", "")
    
    # Extract experiment-specific parameters
    exp_suffix = "_GLUEX" if Experiment == "GlueX" else "_CPP"
    exp_params = {
        "photon_energy_min_GeV": float(experiment_parameters[f"RBHG_ELO{exp_suffix}"]),
        "photon_energy_max_GeV": float(experiment_parameters[f"RBHG_EHI{exp_suffix}"]),
        "target_atomic_number": int(experiment_parameters[f"RBHG_ZTGT{exp_suffix}"]),
        "electron_beam_energy_GeV": float(experiment_parameters[f"RBHG_E0{exp_suffix}"]),
        "coherent_peak_edge_GeV": float(experiment_parameters[f"RBHG_ECOHERENT{exp_suffix}"]),
        "polarization_fraction": float(experiment_parameters[f"RBHG_POL{exp_suffix}"]),
        "theta_min_deg": float(experiment_parameters[f"RBHG_THETAMIN{exp_suffix}"]),
        "theta_max_deg": float(experiment_parameters[f"RBHG_THETAMAX{exp_suffix}"])
    }
    
    # Determine run number if available
    run_number = None
    if RUN_PERIOD:
        try:
            DEFAULTS = {"CPP": "2022-05", "GlueX": "2018-08"}
            run_number = get_run_number(Experiment, RUN_PERIOD, PolDeg, default_periods=DEFAULTS)
        except:
            run_number = None
    
    # Build downstream output directory paths
    downstream_paths = {}
    if directories:
        # Create the study-specific directory pattern for this configuration
        radiation_mode = get_radiation_mode(internal_radiation, single_radiation, hypgeom_radiation)
        config_dir_pattern = f"{study_name}/{nametag}_{form_factor}/{PolDeg}DEG_{radiation_mode}"
        
        # Get the parent directory of the current generation output (where siblings would go)
        gen_output_parent = os.path.dirname(directories['output_directory'])
        
        # Process each output directory configuration
        mcwrapper_base = replace_path_placeholders(
            get_config("MCWRAPPER_OUTPUT_DIR_BASE", "/volatile/halld/home/{USERNAME}/{GENERATOR_TYPE}/"),
            GENERATOR_TYPE, CUEusername, FrameworkHomeDirectory
        )
        dselector_base = get_config("DSELECTOR_OUTPUT_DIR_BASE", "SIBLING")
        tmva_base = get_config("TMVA_OUTPUT_DIR_BASE", "SIBLING")
        
        # Resolve MCWrapper output directory
        if mcwrapper_base == "NESTED":
            mcwrapper_output = os.path.join(directories['output_directory'], "simulation")
        elif mcwrapper_base == "SIBLING":
            mcwrapper_output = os.path.join(gen_output_parent, "simulation")
        else:
            # Explicit base directory - build full structure
            mcwrapper_output = os.path.join(mcwrapper_base, config_dir_pattern, "simulation")
        
        # Resolve DSelector output directory
        if dselector_base == "NESTED":
            dselector_output = os.path.join(mcwrapper_output, "DSelector")
        elif dselector_base == "SIBLING":
            if mcwrapper_base == "NESTED":
                # If MCWrapper is nested, DSelector siblings go alongside generation output
                dselector_output = os.path.join(gen_output_parent, "DSelector")
            else:
                # DSelector siblings go alongside MCWrapper output
                dselector_output = os.path.join(os.path.dirname(mcwrapper_output), "DSelector")
        else:
            # Explicit base directory
            dselector_output = os.path.join(dselector_base, config_dir_pattern, "DSelector")
        
        # Resolve TMVA output directory 
        if tmva_base == "NESTED":
            tmva_output = os.path.join(dselector_output, "TMVA")
        elif tmva_base == "SIBLING":
            if dselector_base in ["NESTED", "SIBLING"]:
                # TMVA siblings go alongside DSelector
                tmva_output = os.path.join(os.path.dirname(dselector_output), "TMVA")
            else:
                # TMVA siblings go alongside DSelector output
                tmva_output = os.path.join(os.path.dirname(dselector_output), "TMVA")
        else:
            # Explicit base directory
            tmva_output = os.path.join(tmva_base, config_dir_pattern, "TMVA")
        
        downstream_paths = {
            "mcwrapper_simulation": {
                "output_directory": mcwrapper_output,
                "trees_directory": os.path.join(mcwrapper_output, "trees"),
                "hists_directory": os.path.join(mcwrapper_output, "hists"),
                "rest_directory": os.path.join(mcwrapper_output, "rest"),
                "input_hddm_directory": directories['vectorspath']  # Points to generation vectors/HDDM
            },
            "dselector_analysis": {
                "output_directory": dselector_output,
                "trees_directory": os.path.join(dselector_output, "trees"),
                "hists_directory": os.path.join(dselector_output, "hists"),
                "input_trees_directory": os.path.join(mcwrapper_output, "trees")  # Points to MCWrapper trees
            },
            "tmva_analysis": {
                "output_directory": tmva_output,
                "models_directory": os.path.join(tmva_output, "models"),
                "results_directory": os.path.join(tmva_output, "results"),
                "input_trees_directory": os.path.join(dselector_output, "trees")  # Points to DSelector trees
            }
        }
    
    config = {
        "rbhg_config": {
            "version": "1.0",
            "generation_info": {
                "study_name": study_name,
                "nametag": nametag,
                "study_type": study_type,
                "timestamp": datetime.now().isoformat(),
                "generator_version": generator_version,
                "script_version": "swif2_RBHG.py"
            },
            "physics_settings": {
                "experiment": Experiment,
                "lepton_type": lepton,
                "xsection_model": BH_xsctn_formulation,
                "form_factor": form_factor,
                "polarization_deg": PolDeg.replace("DEG", "") if "DEG" in str(PolDeg) else str(PolDeg),
                "run_period": RUN_PERIOD or None,
                "run_number": run_number,
                "target": target,
                "internal_radiation": internal_radiation,
                "single_radiation": single_radiation,
                "hypgeom_radiation": hypgeom_radiation,
                "radiation_mode": get_radiation_mode(internal_radiation, single_radiation, hypgeom_radiation),
                "coherent_brems": CobremsDistribution,
                "experiment_parameters": exp_params
            },
            "event_counts": {
                "total_events": int(nevents_total),
                "events_per_file": int(nevents_perfile),
                "total_jobs": math.ceil(int(nevents_total)/int(nevents_perfile)),
                "mode": "fixed_count" if FIXED_EVENT_COUNT else ("run_period" if RUN_PERIOD else "custom"),
                "tenx_flag": TENX_FLAG
            },
            "output_settings": {
                "histograms": histogram_settings,
                "ascii_output": output_settings.get("ascii_vector_output", True),
                "hddm_output": True,  # Since we're using ascii2hddm now
                "integral_xsection": output_settings.get("integral_xsection", False),
                "verbose_output": output_settings.get("verbose_output", True)
            },
            "directory_paths": {
                "base_paths": {
                    "framework_home": FrameworkHomeDirectory if directories else None,
                    "generator_type": GENERATOR_TYPE if directories else None,
                    "farm_out_base": f"/farm_out/{os.getenv('USER', 'acschick')}/{GENERATOR_TYPE}/{study_name}" if directories else None,
                    "output_base": os.path.join(FrameworkHomeDirectory, "output", GENERATOR_TYPE) if directories else None
                },
                "workflow_name": f"{study_name}_{nametag}_{form_factor}_{PolDeg}DEG_{lepton}" if directories else None,
                "rbhg_generation": {
                    "output_directory": directories['output_directory'] if directories else None,
                    "vectors_hddm_directory": directories['vectorspath'] if directories else None,  # Same dir for vectors and HDDM
                    "histograms_directory": directories['histspath'] if directories else None,
                    "fortran_files_directory": directories['submit_directory'] if directories else None,
                    "template_directory": directories['template_directory'] if directories else None,
                    "farm_logs": f"generation/{f'{study_name}_{nametag}_{form_factor}_{PolDeg}DEG_{lepton}'}" if directories else None
                },
                "downstream_analysis": downstream_paths if directories else {},
                "pipeline_farm_logs": {
                    "generation": f"generation/{f'{study_name}_{nametag}_{form_factor}_{PolDeg}DEG_{lepton}'}" if directories else None,
                    "simulation": f"simulation/{f'{study_name}_{nametag}_{form_factor}_{PolDeg}DEG_{lepton}'}" if directories else None,
                    "dselector": f"DSelector/{f'{study_name}_{nametag}_{form_factor}_{PolDeg}DEG_{lepton}'}" if directories else None,
                    "tmva": f"TMVA/{f'{study_name}_{nametag}_{form_factor}_{PolDeg}DEG_{lepton}'}" if directories else None,
                    "analysis_hists": f"analysis_hists/{f'{study_name}_{nametag}_{form_factor}_{PolDeg}DEG_{lepton}'}" if directories else None
                },
                "swif2_workflow_names": {
                    "generation": f"gen_{study_name}_{nametag}_{form_factor}_{run_period}_{PolDeg}DEG_{lepton}" if directories else None,
                    "simulation": f"sim_{study_name}_{nametag}_{form_factor}_{run_period}_{PolDeg}DEG_{lepton}" if directories else None,
                    "dselector": f"dsel_{study_name}_{nametag}_{form_factor}_{run_period}_{PolDeg}DEG_{lepton}" if directories else None,
                    "tmva": f"tmva_{study_name}_{nametag}_{form_factor}_{run_period}_{PolDeg}DEG_{lepton}" if directories else None,
                    "analysis_hists": f"hist_{study_name}_{nametag}_{form_factor}_{run_period}_{PolDeg}DEG_{lepton}" if directories else None
                },
                "templates": {
                    "rbhg_template_directory": directories['template_directory'] if directories else None,
                    "mcwrapper_template_directory": os.path.join(os.path.dirname(directories['template_directory']), "template_MCWrapper") if directories else None
                },
                "scripts": {
                    "ascii2hddm_script": ascii2hddm_script
                }
            },
            "downstream_configs": {
                "mc_simulation": {
                    "events_per_mc_file": int(nevents_perfile),
                    "target": target,
                    "input_hddm_pattern": "vectors_*.hddm",  # Pattern for input HDDM files
                    "expected_input_files": math.ceil(int(nevents_total)/int(nevents_perfile))  # Number of HDDM files expected
                }
            },
            "mcwrapper_settings": mcwrapper_settings  # Add MCWrapper configuration from RunPeriods.json
        }
    }
    return config

def write_rbhg_config(config, output_directory):
    """Write RBHG configuration to JSON file"""
    config_file = os.path.join(output_directory, "rbhg_config.json")
    with open(config_file, 'w') as f:
        json.dump(config, f, indent=2)
    print(f"Configuration saved to: {config_file}")
    return config_file

def determine_study_type(RUN_PERIOD, study_name):
    """
    Determine study type based on run period preset
    
    Logic:
    - _FFS suffix → "Form Factor Study" (can be with or without _AMO)
    - Specific polarization presets (1801, 1808, 1701, FULL2018) WITHOUT _AMO → "Polarization Study"
    - Anything with _AMO (except _FFS) → use study_name (not a standard polarization study)
    - Empty/null preset → use study_name as study_type
    - Other cases → use study_name as study_type
    
    Future: Completely different study types (ECAL, etc.) will use different 
    configuration systems, not RUN_PERIOD extensions.
    """
    if not RUN_PERIOD:
        return study_name
    
    if "_FFS" in RUN_PERIOD:
        return "Form Factor Study"
    
    # If it has _AMO but not _FFS, it's not a standard polarization study
    if "_AMO" in RUN_PERIOD:
        return study_name
    
    # Check for specific polarization study presets (without _AMO)
    polarization_presets = ["1801", "1808", "1701", "FULL2018"]
    
    if RUN_PERIOD in polarization_presets:
        return "Polarization Study"
    
    # Default: use study name
    return study_name

def create_ffs_master_config(study_name, RUN_PERIOD, dataset_configs, base_output_dir):
    """Create master configuration file for Form Factor Study"""
    
    # Parse the run period to understand what was requested
    include_amo = "_AMO" in RUN_PERIOD
    base_period = RUN_PERIOD.replace("_AMO", "").replace("_FFS", "")
    study_type = determine_study_type(RUN_PERIOD, study_name)
    
    master_config = {
        "ffs_master_config": {
            "version": "1.0",
            "study_info": {
                "study_name": study_name,
                "study_type": study_type,
                "preset": RUN_PERIOD if RUN_PERIOD else None,
                "timestamp": datetime.now().isoformat(),
                "generator_script": os.path.abspath(__file__)
            },
            "datasets": {
                "qDATAq_FFN": {
                    "description": "Pseudo-data mimicking real experimental data",
                    "form_factor": "FFN",
                    "statistics_multiplier": 1,
                    "purpose": "Reference dataset for comparison",
                    "config_file": "qDATAq_FFN/run_period_master_config.json"
                },
                "BC_FFN": {
                    "description": "Bin correction with dipole form factor",
                    "form_factor": "FFN", 
                    "statistics_multiplier": 10,
                    "purpose": "Account for q² slope across bins with dipole FF",
                    "config_file": "BC_FFN/run_period_master_config.json"
                },
                "BC_FF1": {
                    "description": "Bin correction with flat form factor",
                    "form_factor": "FF1",
                    "statistics_multiplier": 10, 
                    "purpose": "Account for q² slope across bins with flat FF",
                    "config_file": "BC_FF1/run_period_master_config.json"
                },
                "SIM_FF1": {
                    "description": "Additional simulation with flat form factor",
                    "form_factor": "FF1",
                    "statistics_multiplier": 10,
                    "purpose": "Independent simulation check for systematic studies",
                    "config_file": "SIM_FF1/run_period_master_config.json"
                }
            },
            "dataset_summary": {
                "total_datasets": len(dataset_configs),
                "total_directories": sum(len(config.get("rbhg_config", {}).get("event_counts", {})) 
                                      for config in dataset_configs.values()),
                "dataset_names": list(dataset_configs.keys()),
                "analysis_workflow": [
                    "1. Event generation (swif2_RBHG.py/lepton_event_vXXX.f)", 
                    "2. MC simulation (prepareSimulation.py/MCWrapper)",
                    "3. DSelector analysis", 
                    "4. TMVA application",
                    "5. Final histogram production",
                    "6. Form factor extraction"
                ]
            },
            "batch_processing": {
                "prepSim_command": f"python prepareSimulation.py {study_name}/ffs_master_config.json",
                "individual_configs": ["qDATAq_FFN/run_period_master_config.json",
                                     "BC_FFN/run_period_master_config.json",
                                     "BC_FF1/run_period_master_config.json",
                                     "SIM_FF1/run_period_master_config.json"]
            }
        }
    }
    
    return master_config

def write_ffs_master_config(master_config, study_name, base_output_dir):
    """Write master FFS configuration to top-level directory"""
    study_output_dir = os.path.join(base_output_dir, study_name)
    master_config_file = os.path.join(study_output_dir, "ffs_master_config.json") 
    
    with open(master_config_file, 'w') as f:
        json.dump(master_config, f, indent=2)
    
    print(f"Master FFS configuration saved to: {master_config_file}")
    return master_config_file

def create_run_period_master_config(study_name, dataset_name, form_factor, RUN_PERIOD, polarization_configs, base_output_dir):
    """Create Level 2 master configuration for a run period (coordinates multiple polarizations)"""
    
    # Parse the run period to understand what was requested
    include_amo = "_AMO" in RUN_PERIOD
    base_period = RUN_PERIOD.replace("_AMO", "").replace("_FFS", "")
    study_type = determine_study_type(RUN_PERIOD, study_name)
    
    # Count total directories and extract run period info
    total_directories = sum(len(configs) for configs in polarization_configs.values())
    run_periods = list(polarization_configs.keys())
    
    master_config = {
        "run_period_master_config": {
            "version": "1.0",
            "study_info": {
                "study_name": study_name,
                "dataset_name": dataset_name,
                "form_factor": form_factor,
                "study_type": study_type,
                "preset": RUN_PERIOD if RUN_PERIOD else None,
                "timestamp": datetime.now().isoformat(),
                "generator_script": os.path.abspath(__file__)
            },
            "run_periods": {
                rp: {
                    "polarizations": list(configs.keys()),
                    "total_directories": len(configs),
                    "includes_amo": "AMO" in configs
                } for rp, configs in polarization_configs.items()
            },
            "coordination_summary": {
                "total_run_periods": len(run_periods),
                "total_polarizations": sum(len(configs) for configs in polarization_configs.values()),
                "total_directories": total_directories,
                "directories": [f"{rp}_{pol}_{study_name}_{dataset_name}_{form_factor}_{config['rbhg_config']['physics_settings']['radiation_mode']}" 
                              for rp, configs_by_pol in polarization_configs.items() 
                              for pol, config in configs_by_pol.items()],
                "analysis_workflow": [
                    "1. Event generation (swif2_RBHG.py)",
                    "2. MC simulation (prepareSimulation.py/MCWrapper)", 
                    "3. DSelector analysis",
                    "4. TMVA application", 
                    "5. Final histogram production"
                ]
            },
            "batch_processing": {
                "prepSim_command": f"python prepareSimulation.py --run-period-config {study_name}/{dataset_name}_{form_factor}/run_period_master_config.json",
                "individual_configs": [f"{rp}_{pol}_*/rbhg_config.json" 
                                     for rp, configs in polarization_configs.items() 
                                     for pol in configs.keys()]
            }
        }
    }
    
    return master_config

def write_run_period_master_config(master_config, study_name, dataset_name, form_factor, base_output_dir):
    """Write Level 2 run period master configuration"""
    dataset_output_dir = os.path.join(base_output_dir, study_name, f"{dataset_name}_{form_factor}")
    master_config_file = os.path.join(dataset_output_dir, "run_period_master_config.json")
    
    # Create directory if it doesn't exist
    os.makedirs(dataset_output_dir, exist_ok=True)
    
    with open(master_config_file, 'w') as f:
        json.dump(master_config, f, indent=2)
    
    print(f"Run period master configuration saved to: {master_config_file}")
    return master_config_file

def load_runperiods_data(run_period, polarization, framework_home_dir):
    """
    Load data from RunPeriods.json for a specific run period and polarization.
    Returns a dict with polarization-specific and run-period-level data, or None if not found.
    """
    try:
        import json
        runperiods_path = os.path.join(framework_home_dir, 'RunPeriods.json')
        if not os.path.exists(runperiods_path):
            print(f"Warning: RunPeriods.json not found at {runperiods_path}")
            return None
        
        with open(runperiods_path, 'r') as f:
            runperiods_data = json.load(f)
        
        # Parse run_period (strip FFS/AMO suffixes, convert FULL2018 -> 1808)
        base_period = run_period.replace("_AMO", "").replace("_FFS", "")
        if base_period == "FULL2018":
            base_period = "1808"  # Use 1808 as representative for FULL2018
        
        # Normalize polarization format
        pol_key = str(polarization).upper().replace("DEGREES", "DEG").replace("DEGREE", "DEG")
        if not pol_key.endswith("DEG") and pol_key not in ["AMO"]:
            pol_key = pol_key + "DEG"
        
        # Get run period data
        if base_period not in runperiods_data:
            print(f"Warning: Run period {base_period} not found in RunPeriods.json")
            return None
        
        period_data = runperiods_data[base_period]
        
        # Get polarization-specific data
        pol_data = {}
        if "Polarizations" in period_data and pol_key in period_data["Polarizations"]:
            pol_data = period_data["Polarizations"][pol_key]
        else:
            print(f"Warning: Polarization {pol_key} not found in RunPeriods.json for run period {base_period}")
        
        # Return combined data (polarization-specific overrides run-period defaults)
        result = {
            "cobrems_distribution": pol_data.get("cobrems_distribution", ""),
            "epem_event_count": pol_data.get("epem_event_count", ""),
            "characteristic_run": pol_data.get("characteristic_run", ""),
            "total_triggers": pol_data.get("total_triggers", ""),
            "recon_env": period_data.get("recon_env", ""),
            "analysis_env": period_data.get("analysis_env", ""),
            "background_version": period_data.get("background_version", ""),
            "jana_epem_config": period_data.get("jana_epem_config", ""),
            "RCDB_query": period_data.get("RCDB_query", ""),
            "batch_system": period_data.get("batch_system", ""),
            "run_range": period_data.get("run_range", "")
        }
        
        return result
        
    except Exception as e:
        print(f"Warning: Error loading RunPeriods.json: {e}")
        return None

# Accounting for acceptance, cuts, and other attrition. Determined empirically.

# Reminder: the nevents for an orientation/run period in this dictionary don't need to be strings
# because only the nevent_perfile is used in the regex replace for the individual fortran files.
# These numbers are only used to compute the total number of jobs through nevents/nevents_perfile
nevents = {
    "1801": {
        "0DEG":   6600000,
        "90DEG":  6300000,
        "135DEG": 6600000,
        "45DEG":  6800000,
        "AMO":    6600000
    },
    "1808": {
        "0DEG":   4000000,
        "90DEG":  3600000,
        "135DEG": 3100000,
        "45DEG":  3000000,
        "AMO":    3600000
    },
    "2205": {
        "135DEG": 40000000,
        "45DEG":  40000000,
        "AMO":     40000000  # TODO: Update with actual AMO event count when available
    }
}
# e.g. events_1801_45DEG = nevents["1801"]["45DEG"]

# Load event counts from RunPeriods.json when FIXED_EVENT_COUNT is False
if not FIXED_EVENT_COUNT and RUN_PERIOD:
    print("FIXED_EVENT_COUNT is False - loading event counts from RunPeriods.json...")
    try:
        # Load event counts for each run period and polarization
        for run_period_key in list(nevents.keys()):
            for pol_key in list(nevents[run_period_key].keys()):
                # Create full run_period string with potential suffixes
                full_run_period = RUN_PERIOD if RUN_PERIOD else run_period_key
                
                rp_data = load_runperiods_data(run_period_key, pol_key, FrameworkHomeDirectory)
                if rp_data and rp_data.get("epem_event_count"):
                    event_count = rp_data["epem_event_count"]
                    # Convert to int (might be string in JSON)
                    try:
                        nevents[run_period_key][pol_key] = int(event_count)
                        print(f"  {run_period_key}/{pol_key}: {nevents[run_period_key][pol_key]:,} events (from RunPeriods.json)")
                    except (ValueError, TypeError) as e:
                        print(f"  Warning: Could not parse event count for {run_period_key}/{pol_key}: {event_count}")
    except Exception as e:
        print(f"Warning: Could not load event counts from RunPeriods.json: {e}")
        print("Using hardcoded event counts as fallback.")

# overwrite the number of events if using fixed event count:
if FIXED_EVENT_COUNT:
    # Use the same NEVENTS_TOTAL value for all run periods and polarizations
    for run_period in nevents:
        for polarization in nevents[run_period]:
            nevents[run_period][polarization] = int(nevents_total)
if TENX_FLAG:
    for run_period in nevents:
        for polarization in nevents[run_period]:
            nevents[run_period][polarization] *= 10
    # Also 10x the events per file
    nevents_perfile = str(int(nevents_perfile) * 10)
    print(f"TENX_FLAG is True: Multiplied total events and events per file by 10")
    print(f"  New nevents_perfile = {nevents_perfile}")





# Create dict for loop in main()
def get_run_configurations(RUN_PERIOD):
    """
    Parse RUN_PERIOD to extract base period and suffixes.
    Supports: _AMO (include AMO data), _FFS (Form Factor Study - 4 datasets)
    Examples: FULL2018, FULL2018_AMO, FULL2018_FFS, FULL2018_AMO_FFS
    """
    # Parse suffixes
    include_amo = "_AMO" in RUN_PERIOD
    include_ffs = "_FFS" in RUN_PERIOD
    
    # Extract base period by removing all suffixes
    base_period = RUN_PERIOD.replace("_AMO", "").replace("_FFS", "")
    
    # Get base run configurations
    if base_period == "FULL2018":
        # FULL2018: Both 1801 and 1808 run periods
        base_config = {}
        for run_period in ("1801", "1808"):
            if include_amo:
                base_config[run_period] = nevents[run_period]  # includes AMO
            else:
                base_config[run_period] = {k: v for k, v in nevents[run_period].items() if k != "AMO"}  # excludes AMO
    else:
        # Single run period (1801, 1808, 2205, etc.)
        if base_period in nevents:
            if include_amo:
                base_config = {base_period: nevents[base_period]}  # includes AMO if it exists
            else:
                base_config = {base_period: {k: v for k, v in nevents[base_period].items() if k != "AMO"}}  # excludes AMO
        else:
            return {}  # unknown run period
    
    # If FFS (Form Factor Study) is requested, return special structure
    if include_ffs:
        return {"FFS": base_config}  # Mark as Form Factor Study
    else:
        return base_config



def get_run_number(EXP, RUNPERIOD, POL, default_periods=None) -> str:
    """
    Robust lookup with normalization and optional defaults.
      - EXP: case-insensitive ("gluex", "gx", "cpp")
      - RUNPERIOD: accepts "", "2018-01", "201801", "1801", "2018/01"
                   If empty and default_periods provided, uses that.
      - POL: accepts 0/45/90/135, "0", "0deg", "0°", "AMO" (any case)
    """

    RunNumbers = {
        "GlueX": {
            "2017-01": {"0DEG": "30327", "90DEG": "30891", "135DEG": "30480", "45DEG": "30484", "AMO": "30888"},
            "2018-01": {"0DEG": "42559", "90DEG": "41541", "135DEG": "41172", "45DEG": "41425", "AMO": "42554"},
            "2018-08": {"0DEG": "50771", "90DEG": "50699", "135DEG": "51233", "45DEG": "51056", "AMO": "51029"},
        },
        "CPP": {
            "2022-05": {"135DEG": "101582", "45DEG": "101586", "AMO": "101615"}
        }
    }

    # ---------- Normalizers (no free vars) ----------
    def norm_expt(x: str) -> str:
        s = str(x).strip().lower()
        if s in {"gluex", "gx"}: return "GlueX"
        if s == "cpp":           return "CPP"
        choices = ["GlueX", "CPP"]
        sugg = difflib.get_close_matches(s, [c.lower() for c in choices], n=1)
        hint = f" (did you mean '{choices[[c.lower() for c in choices].index(sugg[0])]}')?" if sugg else ""
        raise ValueError(f"Unknown experiment '{x}'. Choices: {choices}{hint}")

    def period_choices(expt: str):
        return sorted(RunNumbers[expt].keys())

    def norm_period(expt: str, p, default_periods) -> str:
        s = "" if p is None else str(p).strip()
        # If empty, try default
        if not s:
            if default_periods and expt in default_periods:
                s = default_periods[expt]
            else:
                raise ValueError(f"Missing run period for {expt}. Choices: {period_choices(expt)}")

        # Remove separators for flexible parsing
        s_flat = s.replace(" ", "").replace("_", "").replace("/", "").replace("-", "")

        # YYYYMM -> YYYY-MM
        if re.fullmatch(r"\d{6}", s_flat):
            cand = f"{s_flat[:4]}-{s_flat[4:]}"
            if cand in RunNumbers[expt]: return cand

        # YYMM -> 20YY-MM (e.g., 1801 -> 2018-01, 2205 -> 2022-05)
        if re.fullmatch(r"\d{4}", s_flat):
            cand = f"20{s_flat[:2]}-{s_flat[2:]}"
            if cand in RunNumbers[expt]: return cand

        # YYYY or YYYY-M(M)
        m = re.fullmatch(r"(\d{4})(?:[-]?(\d{1,2}))?$", s.replace("/", "-"))
        if m:
            year, month = m.groups()
            if month is None:
                matches = [rp for rp in RunNumbers[expt] if rp.startswith(year + "-")]
                if len(matches) == 1:
                    return matches[0]
                raise ValueError(f"Ambiguous period '{s}' for {expt}. Available for {year}: {sorted(matches)}")
            cand = f"{year}-{int(month):02d}"
            if cand in RunNumbers[expt]: return cand

        # Already canonical?
        if s in RunNumbers[expt]: return s

        sugg = difflib.get_close_matches(s, period_choices(expt), n=1)
        hint = f" (did you mean '{sugg[0]}'?)" if sugg else ""
        raise ValueError(f"Unknown run period '{s}' for {expt}. Choices: {period_choices(expt)}{hint}")

    def norm_pol(expt: str, period: str, pol) -> str:
        # numeric like 0, 45, 90, 135
        if isinstance(pol, (int, float)) and float(pol).is_integer():
            key = f"{int(pol)}DEG"
            if key in RunNumbers[expt][period]: return key
        s = str(pol).strip().upper().replace("°", "DEG")
        s = s.replace("DEGREES", "DEG").replace("DEGREE", "DEG")
        if s in {"0", "45", "90", "135"}: s = f"{int(s)}DEG"
        if s.endswith("DEG") and s[:-3].isdigit(): s = f"{int(s[:-3])}DEG"
        if s == "AMO": s = "AMO"
        if s in RunNumbers[expt][period]: return s
        choices = sorted(RunNumbers[expt][period].keys())
        sugg = difflib.get_close_matches(s, choices, n=1)
        hint = f" (did you mean '{sugg[0]}'?)" if sugg else ""
        raise ValueError(f"Unknown polarization '{pol}' for {expt} {period}. Choices: {choices}{hint}")

    # ---------- Apply in safe order ----------
    expt = norm_expt(EXP)
    period = norm_period(expt, RUNPERIOD, default_periods or {})
    pol_key = norm_pol(expt, period, POL)
    return RunNumbers[expt][period][pol_key]






##############################################
########### DIRECTORIES AND NAMES ############
##############################################
def create_ifnot_directory(Dir):
    if not os.path.exists(Dir):
        os.makedirs(Dir)

def create_directory_plusone(Dir):
    i_outdir = 0
    outdir = Dir
    while os.path.exists(outdir):
        i_outdir += 1
        outdir = f"{Dir}_{i_outdir}"
    os.makedirs(outdir)
    return outdir

def check_directory(Dir):
    if not os.path.exists(Dir):
        print(f"{Dir} does not exist. Exiting..")
        sys.exit(1)

def get_radiation_mode(internal_radiation, single_radiation, hypgeom_radiation):
    """
    Determine the radiation mode string based on radiation settings
    
    Args:
        internal_radiation (bool): Mo & Tsai internal radiation
        single_radiation (bool): Single lepton radiation (v113 feature)
        hypgeom_radiation (bool): Hypergeometric radiation (future)
    
    Returns:
        str: Radiation mode string (IRADOFF, ONERAD, DBLRAD, HYPGEOM)
    """
    if hypgeom_radiation:
        return "HYPGEOM"
    elif not internal_radiation:
        return "IRADOFF"
    elif single_radiation:
        return "ONERAD"
    else:
        return "DBLRAD"

def abbreviate_number(s):
    num = float(s)
    # Check for divisibility by 1k, 1mil, etc., and format accordingly
    if num < 1_000:
        return s
    elif num % 1_000 == 0 and num < 1_000_000:
        return f"{num / 1_000:.0f}k"
    elif num % 1_000_000 == 0 and num < 1_000_000_000:
        return f"{num / 1_000_000:.0f}mil"
    elif num % 1_000_000_000 == 0 and num < 1_000_000_000_000:
        return f"{num / 1_000_000_000:.0f}bil"
    elif num % 1_000_000_000_000 == 0:
        return f"{num / 1_000_000_000_000:.0f}tril"
    else:
        return s  # Return the original string if not evenly divisible


def make_generation_dirname(study_name, nametag, form_factor, PolDeg, run_period, lepton, BH_xsctn_formulation, nevents_total, 
                            RUNPERIOD, TENXFLAG, Experiment, RBHGfortranVer, internal_radiation, single_radiation=False, hypgeom_radiation=False):
    # Handle the TENXFLAG for the '_10X_' part
    tenx = "_10X" if TENXFLAG else ""
    # Determine radiation mode
    radiation_mode = get_radiation_mode(internal_radiation, single_radiation, hypgeom_radiation)
    intrad = f"_{radiation_mode}"

    # Save original PolDeg for directory naming (before Heitler blanks it)
    PolDeg_for_dirname = PolDeg
    
    if BH_xsctn_formulation.lower() == "heitler":
        # For Heitler (unpolarized), blank PolDeg for physics but keep original for directory name
        PolDeg = ""
    elif "DEG" not in PolDeg.upper(): # elifs only proceed if first 'if' is not true
        PolDeg += "DEG"
        PolDeg_for_dirname = PolDeg  # Update dirname version too

    # Determine if this is a Form Factor Study
    is_form_factor_study = RUNPERIOD and "_FFS" in RUNPERIOD
    
    dirParts = []
    dirParts.append(study_name)

    # For organized structure, individual configs should go inside dataset folders
    if nametag:
        # Always use underscore between nametag and form_factor (consistent naming)
        dataset_folder = f"{nametag}_{form_factor}"
        dirParts.append(dataset_folder)
    else:
        dirParts.append(form_factor)

    if is_form_factor_study:
        # Form Factor Study: Keep intermediate level for dataset distinction
        # Build final directory for Form Factor Study
        final_dir_parts = []
    else:
        # Regular study: Build individual directory name without dataset prefix
        final_dir_parts = []

    # Build the final directory name
    if is_form_factor_study:
        # Form Factor Study: Include form factor at run_period level for clarity
        if not RUNPERIOD:
            final_dir_parts.append(abbreviate_number(nevents_total) + "Events")
            final_dir_parts.append(f"{PolDeg_for_dirname}_{form_factor}{intrad}")
        else:
            if not RUNPERIOD.upper() in ["FULL2018", "FULL2018_AMO"]:
                final_dir_parts.append(f"{run_period}_{PolDeg_for_dirname}_{form_factor}{intrad}")
            else:
                final_dir_parts.append(f"{run_period}_{PolDeg_for_dirname}_{form_factor}{intrad}")
    else:
        # Regular study: individual directory name without redundant dataset prefix
        if not RUNPERIOD:
            final_dir_parts.append(abbreviate_number(nevents_total) + "Events")
        else:
            if not RUNPERIOD.upper() in ["FULL2018", "FULL2018_AMO"]:
                final_dir_parts.append(run_period)
            else:
                final_dir_parts.append(run_period)
        
        # Add polarization and radiation mode
        if PolDeg_for_dirname:
            final_dir_parts.append(f"{PolDeg_for_dirname}{intrad}")
        else:
            final_dir_parts.append(intrad.lstrip("_"))  # Remove leading underscore if no polarization
        
        # Add technical details for regular studies
        final_dir_parts.extend([Experiment, lepton, BH_xsctn_formulation, RBHGfortranVer])
    
    # Build final directory name and complete path
    final_dir_name = '_'.join(filter(None, final_dir_parts))
    dirParts.append(final_dir_name)
    directory_name = '/'.join(dirParts)

    return directory_name


def make_all_gen_directories(FrameworkHomeDirectory, gen_dir_name, CUEusername, logicals, run_period=None, polarization=None):
    outputDirectoryTop = os.path.join(FrameworkHomeDirectory, "output")
    template_directory = os.path.join(FrameworkHomeDirectory, "generators", "RBHG", "template_RBHG")
    check_directory(template_directory)
    # Create output directory and generator-specific subdirectory if they don't exist
    create_ifnot_directory(outputDirectoryTop)
    outputDirectoryGenerator = os.path.join(outputDirectoryTop, GENERATOR_TYPE)
    create_ifnot_directory(outputDirectoryGenerator)
    gen_output_directory = os.path.join(outputDirectoryGenerator, gen_dir_name)
    # Extract study name from gen_dir_name for farm_out structure
    study_name_from_dir = gen_dir_name.split('/')[0]  # First part of path is study name
    workflow_name_base = '_'.join(gen_dir_name.split('/')[1:])  # Rest becomes workflow name
    # Extract study name from gen_dir_name for farm_out structure
    study_name_from_dir = gen_dir_name.split('/')[0]  # First part of path is study name
    workflow_name_base = '_'.join(gen_dir_name.split('/')[1:])  # Rest becomes workflow name
    farm_out_directory = os.path.join(f"/farm_out/{CUEusername}/{GENERATOR_TYPE}/{study_name_from_dir}/generation", workflow_name_base)
    
    out_dir = create_directory_plusone(gen_output_directory)
    sub_dir = os.path.join(f"{out_dir}", "FortranFiles")
    create_ifnot_directory(sub_dir)
    farm_out_dir = create_directory_plusone(farm_out_directory)
    
    # Dictionaries
    dir_name = gen_output_directory    
    Directories = {
        'logpath': farm_out_dir,
        'vectorspath': os.path.join(out_dir, 'vectors'),
        'histspath': os.path.join(out_dir, 'hists'),
    }


    RBHGOUTPUTDICT = {
        "FRAMEWORKHOMEDIRECTORY": FrameworkHomeDirectory,
        "RBHGOUTPUTDIRTOP": outputDirectoryTop,
        "RBHGVECTORSPATH": Directories["vectorspath"].replace(outputDirectoryTop, ''),        # fortran has line length limit--had to break up paths to make it fit (hence replace)
        "RBHGHISTSPATH": Directories["histspath"].replace(outputDirectoryTop, '')             # fixed length limit problem again.
    }

    # Check if RBHGOUTPUT_EVENT is ".false."                                                                                                                 
    if logicals["RBHGOUTPUT_EVENT"] == ".false.":
        del Directories['vectorspath']
        RBHGOUTPUTDICT['RBHGVECTORSPATH'] = ""
        print("No output 3-vector file requested.")

    # Check if all logicals containing 'HIST' are ".false."                                                                                                  
    if all(logicals[key] == ".false." for key in logicals if 'HIST' in key):
        del Directories['histspath']
        RBHGOUTPUTDICT['RBHGHISTSPATH'] = ""
        print("No histograms requested.")

    for Dir in Directories.values():
        create_ifnot_directory(Dir)

    Directories2 = {
        'submit_directory': sub_dir,
        'output_directory': out_dir,
        'template_directory': template_directory
    }

    Directories.update(Directories2)
    
    # Add cobrems file override processing - copy file to submit directory for easy access
    cobrems_override = ""
    cobrems_source_path = ""
    
    if CobremsFileGlueX.strip():
        # User explicitly specified a file in RBHG.config
        cobrems_source_path = CobremsFileGlueX if os.path.isabs(CobremsFileGlueX) else os.path.join(FrameworkHomeDirectory, CobremsFileGlueX)
    elif CobremsFileCPP.strip():
        # User explicitly specified a file in RBHG.config
        cobrems_source_path = CobremsFileCPP if os.path.isabs(CobremsFileCPP) else os.path.join(FrameworkHomeDirectory, CobremsFileCPP)
    elif run_period or DEFAULT_RUN_PERIOD:
        # Use explicit run_period or DEFAULT_RUN_PERIOD for cobrems lookup
        lookup_period = run_period if run_period else DEFAULT_RUN_PERIOD
        # No override in config - look up from RunPeriods.json
        try:
            import json
            runperiods_path = os.path.join(FrameworkHomeDirectory, 'RunPeriods.json')
            if os.path.exists(runperiods_path):
                with open(runperiods_path, 'r') as f:
                    runperiods_data = json.load(f)
                
                # Parse run_period (strip FFS/AMO suffixes)
                base_period = lookup_period.replace("_AMO", "").replace("_FFS", "").replace("FULL2018", "1808")
                
                # Normalize polarization format
                pol_key = str(polarization).upper().replace("DEGREES", "DEG").replace("DEGREE", "DEG")
                if not pol_key.endswith("DEG") and pol_key not in ["AMO"]:
                    pol_key = pol_key + "DEG"
                
                # Look up cobrems file path
                if base_period in runperiods_data:
                    period_data = runperiods_data[base_period]
                    if "Polarizations" in period_data and pol_key in period_data["Polarizations"]:
                        cobrems_path = period_data["Polarizations"][pol_key].get("cobrems_distribution", "")
                        if cobrems_path:
                            cobrems_source_path = cobrems_path
                            print(f"Using cobrems file from RunPeriods.json ({base_period}/{pol_key}): {cobrems_source_path}")
        except Exception as e:
            print(f"Warning: Could not lookup cobrems from RunPeriods.json: {e}")
    elif current_CobremsDistribution:
        # User wants cobrems but didn't specify file and no run_period - warn them
        print(f"WARNING: COBREMS_DISTRIBUTION=True but no cobrems file specified and RUN_PERIOD is blank!")
        print(f"  Either:")
        print(f"    1. Set COBREMS_FILE_GLUEX or COBREMS_FILE_CPP in RBHG.config")
        print(f"    2. Set RUN_PERIOD to auto-select from RunPeriods.json")
        print(f"    3. Set COBREMS_DISTRIBUTION=False to use 1/E distribution")
        print(f"  Falling back to 1/E distribution for now.")
        # Don't set a cobrems file - Fortran will use 1/E
    
    # If we have a cobrems file, copy it to CobremsDistributions subdirectory in template_RBHG
    if cobrems_source_path and os.path.exists(cobrems_source_path):
        cobrems_filename = os.path.basename(cobrems_source_path)
        # Copy to template_RBHG/CobremsDistributions/ directory
        cobrems_dest_dir = os.path.join(template_directory, "CobremsDistributions")
        os.makedirs(cobrems_dest_dir, exist_ok=True)
        cobrems_dest_path = os.path.join(cobrems_dest_dir, cobrems_filename)
        subprocess.run(["cp", cobrems_source_path, cobrems_dest_path])
        # Pass just the filename - Fortran will prepend templateDir/CobremsDistributions/
        cobrems_override = cobrems_filename
        print(f"Copied cobrems file to: {cobrems_dest_path}")
    elif cobrems_source_path:
        print(f"Warning: Cobrems file not found: {cobrems_source_path}")
        cobrems_override = ""
    
    # Add cobrems file override to replacement dictionary
    RBHGOUTPUTDICT["RBHGCOBREMS_FILE_OVERRIDE"] = cobrems_override
    
    logicals.update(RBHGOUTPUTDICT)

    return Directories, logicals


### At this point, we have the top level directories all created-- stderr and stdout can all be in one directory, same with hists and vectors.
## maybe just have fortran exe files?

def MakeConfigReport(study_name, nametag, form_factor, PolDeg, run_period, lepton, BH_xsctn_formulation, nevents_total):
    
    
    return 

### Fortran File
def fortran_basename(fortran_file):
    # Using regex to extract full name for renaming template file, and creating .exe's
    match = re.search(r"temp_(.+)\.f", fortran_file)
    if match:
        fortran_file_basename = match.group(1)  # removing .f
        extracted_string = "RBHG version = " + fortran_file_basename
        #print(extracted_string)  # prints "lepton_event_v108"
    else:
        print("Couldn't find RBHG version. Perhaps you named the template file wrong?")

    return fortran_file_basename

        
def determine_RBHG_fortran_ver(fortran_file):
    match = re.search(r'v\d+', fortran_file)
    if match:
        VER = match.group()
    else:
        VER = None

    return VER


def fortranfile_cp_and_regex(seed, job_i, neventsperfile, dir_dict, temp_fortran_file, logicals, experiment_parameters):    
    currentSeed = str(seed)
    sjob_i = str(job_i)
    FortranFile_i = os.path.join(f"{dir_dict['submit_directory']}", fortran_basename(temp_fortran_file) + '_' + str(sjob_i) +".f")
    command = ["cp", f"{dir_dict['template_directory']}/{temp_fortran_file}", f"{FortranFile_i}"]
    subprocess.run(command)

    logicals.update(experiment_parameters)
    logicals.update({
        "RBHGJOBNUMBER": sjob_i,
        "RBHGSEEDVALUE": currentSeed,
        "RBHGNEVENTS": neventsperfile,
    })

    # to prevent RBHGHIST_T from replacing RBHGHIST_T_VARIABLE and RBHGHIST_THETA.
    length_sorted_logicals = sorted(logicals.items(), key=lambda x: len(x[0]), reverse=True)  
    
    # use 'sed' to replace REGEXs in copied template file with values from dictionary                                                              
    for placeholder, replacement in length_sorted_logicals:
        command = ['sed', '-i', f's|{placeholder}|{replacement}|g', f'{FortranFile_i}']
        subprocess.run(command)
            
    compiled_fortran_exe_path = os.path.join(f"{dir_dict['submit_directory']}", fortran_basename(temp_fortran_file) + f"_{sjob_i}.exe")
    command = ['gfortran', '-std=legacy', '-ffixed-line-length-132', '-fd-lines-as-code', f"{FortranFile_i}", '-o', compiled_fortran_exe_path]
    subprocess.run(command)

    return compiled_fortran_exe_path


######################################
########## swif2 functions ###########
######################################

def capture_swif_stdout(command):
    try:
        # If command is a list, join it into a single string
        if isinstance(command, list):
            command = ' '.join(command)
        result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, shell=True)

        if result.returncode != 0:
            print("Error running command:", command, "\nError Message:", result.stderr)
            return None

        return result.stdout
    except Exception as e:
        print("An error occurred while running command:", command, "\nError:", str(e))
        return None

def parse_workflow_names(swif_output):
    # Split the output into lines
    lines = swif_output.strip().split('\n')
    # Skip the header line and possibly the footer or empty lines
    workflow_lines = lines[1:]
    # Extract the workflow names
    workflow_names = []
    for line in workflow_lines:
        columns = line.split('|')
        if columns and len(columns) > 1:
            workflow_name = columns[0].strip()
            workflow_names.append(workflow_name)

    return workflow_names

def workflow_exists(workflow_name):
    swif2_list_output = capture_swif_stdout(['swif2', 'list'])
    workflows = parse_workflow_names(swif2_list_output)
    return workflow_name in workflows

def create_swif2_workflow(workflow_name):
    original_name = workflow_name
    n = 1
    while workflow_exists(workflow_name):
        workflow_name = f"{original_name}_{n}"
        n += 1
    
    subprocess.run(["swif2", "create", "-workflow", workflow_name])
    #print(f"Workflow created with name: {workflow_name}")
    return workflow_name

def safe_workflow_name(name, maxlen=48):
    """Ensure workflow name does not exceed SWIF2 limits; append short hash if truncated."""
    if len(name) <= maxlen:
        return name
    short_hash = hashlib.md5(name.encode('utf-8')).hexdigest()[:6]
    return f"{name[:maxlen-7]}_{short_hash}"


def make_workflow_name(study_name, nametag, form_factor, run_period, PolDeg, lepton, internal_radiation, single_radiation=False, hypgeom_radiation=False):
    parts = []
    radiation_mode = get_radiation_mode(internal_radiation, single_radiation, hypgeom_radiation)
    intrad = radiation_mode
    parts.append(study_name)
    
    # Corrected to use form_factor
    if nametag.upper() == "BC" or nametag.upper() == "SIM":
        parts.append(f"{nametag}{form_factor}")
    else:
        parts.append(f"{nametag}_{form_factor}")

    # Use extend() to add multiple elements at once
    parts.extend([run_period, PolDeg, intrad, lepton])
    
    workflow_name = '_'.join(filter(None, parts))
    return safe_workflow_name(workflow_name)

def swif2_add_job_with_args(stub_dict, workflow_name, ENVFILE, leptonexe, directories, job_i, seed, nevents, lepton, experiment, ascii2hddm_script, hddm_run_number):
    """Add swif2 job that passes seed, job_i, and nevents as command-line arguments to the executable"""
    exp_2_target = {
        "GlueX": "p",
        "CPP": "pb208"
        }
    lepton_2_finalstate = {
        "ee": "epem",
        "mu": "mupmum"
        }
    target = exp_2_target[experiment]
    particle = lepton_2_finalstate[lepton]

    job_name_i = f"{workflow_name}_{job_i}"
    RBHG_bash_script = stub_dict["RBHG_bash_script"]
    RBHG_bash_script_path = os.path.join(directories['template_directory'], RBHG_bash_script)
    stub_dict['RBHG_bash_script'] = RBHG_bash_script_path

    add_job_dict = {
        "-workflow": workflow_name,
        "-name": job_name_i,
        "-stdout": os.path.join(directories['logpath'], f"stdout_{job_name_i}.out"),
        "-stderr": os.path.join(directories['logpath'], f"stderr_{job_name_i}.err")
    }
    add_job_dict.update(stub_dict)

    command_list = ["swif2", "add-job"]

    # Append other options and their values
    for key, value in add_job_dict.items():
        if key.startswith("-"):
            command_list.extend([key, value])
        else:
            # If the key doesn't start with "-", only append the value
            command_list.append(value)

    # Arguments for bash script - now includes seed, job_i, nevents as exe arguments
    command_list.append(ENVFILE)
    # Build Fortran command as single string EXACTLY like hddm_command
    fortran_command = f"{leptonexe} {seed} {job_i} {nevents}"
    command_list.append(fortran_command)
    command_list.append(directories['vectorspath']) # cd to vectors path in bash script
    # Build hddm command as single string (SAME format as fortran_command)
    hddm_command = (
        f"python {ascii2hddm_script} {particle} {target} "
        f"{os.path.join(directories['vectorspath'], f'vectors_{job_i}.txt')} "
        f"{os.path.join(directories['vectorspath'], f'vectors_{job_i}.hddm')} "
        f"--run {hddm_run_number} "
        f"--vertex '0 0 0 0'"
    )
    command_list.append(hddm_command)
    # Execute the command
    subprocess.run(command_list)


def run_job_locally_with_args(ENVFILE, leptonexe, directories, job_i, seed, nevents, lepton, experiment, ascii2hddm_script, hddm_run_number, RBHG_script):
    """Run job locally - directly execute Fortran and HDDM conversion without shell wrapper"""
    exp_2_target = {
        "GlueX": "p",
        "CPP": "pb208"
    }
    lepton_2_finalstate = {
        "ee": "epem",
        "mu": "mupmum"
    }
    target = exp_2_target[experiment]
    particle = lepton_2_finalstate[lepton]
    
    print(f"Running: {leptonexe} {seed} {job_i} {nevents}")
    
    # Run Fortran executable directly
    fortran_command = [leptonexe, str(seed), str(job_i), str(nevents)]
    fortran_result = subprocess.run(fortran_command, capture_output=False, cwd=directories['vectorspath'])
    
    if fortran_result.returncode != 0:
        print(f"ERROR: Fortran executable failed with exit code {fortran_result.returncode}")
        print(f"  ✗ Job {job_i} failed with return code {fortran_result.returncode}")
        return fortran_result.returncode
    
    # Run HDDM conversion
    hddm_command = [
        'python', ascii2hddm_script, particle, target,
        os.path.join(directories['vectorspath'], f'vectors_{job_i}.txt'),
        os.path.join(directories['vectorspath'], f'vectors_{job_i}.hddm'),
        '--run', str(hddm_run_number),
        '--vertex', '0', '0', '0', '0'
    ]
    hddm_result = subprocess.run(hddm_command, capture_output=False)
    
    if hddm_result.returncode != 0:
        print(f"ERROR: HDDM conversion failed with exit code {hddm_result.returncode}")
        print(f"  ✗ Job {job_i} failed with return code {hddm_result.returncode}")
        return hddm_result.returncode
    
    print(f"  ✓ Job {job_i} completed successfully")
    return 0


def swif2_add_job(stub_dict, workflow_name, ENVFILE, leptonexe, directories, job_i, lepton, experiment, ascii2hddm_script, hddm_run_number):
    exp_2_target = {
        "GlueX": "p",
        "CPP": "pb208"
        }
    lepton_2_finalstate = {
        "ee": "epem",
        "mu": "mupmum"
        }
    target = exp_2_target[experiment]
    particle = lepton_2_finalstate[lepton]

    job_name_i = f"{workflow_name}_{job_i}"
    RBHG_bash_script = stub_dict["RBHG_bash_script"]
    RBHG_bash_script_path = os.path.join(directories['template_directory'], RBHG_bash_script)
    stub_dict['RBHG_bash_script'] = RBHG_bash_script_path

    add_job_dict = {
        "-workflow": workflow_name,
        "-name": job_name_i,
        "-stdout": os.path.join(directories['logpath'], f"stdout_{job_name_i}.out"),
        "-stderr": os.path.join(directories['logpath'], f"stderr_{job_name_i}.err")
    }
    add_job_dict.update(stub_dict)

    command_list = ["swif2", "add-job"]

    # Append other options and their values
    for key, value in add_job_dict.items():
        if key.startswith("-"):
            command_list.extend([key, value])
        else:
            # If the key doesn't start with "-", only append the value
            command_list.append(value)

    # Arguments for bash script
    command_list.append(ENVFILE)
    command_list.append(leptonexe)
    command_list.append(directories['vectorspath']) # cd to vectors path in bash script
    hddm_command = (
        f"python {ascii2hddm_script} {particle} {target} "
        f"{os.path.join(directories['vectorspath'], f'vectors_{job_i}.txt')} "
        f"{os.path.join(directories['vectorspath'], f'vectors_{job_i}.hddm')} "
        f"--run {hddm_run_number} "
        f"--vertex '0 0 0 0'"
    )
    command_list.append(hddm_command)
    # Execute the command
    subprocess.run(command_list)


def run_job_locally(ENVFILE, leptonexe, directories, job_i, lepton, experiment, ascii2hddm_script, hddm_run_number):
    """
    Run a job locally without swif2 submission.
    Executes the same bash script that swif2 would run on the farm.
    """
    exp_2_target = {
        "GlueX": "p",
        "CPP": "pb208"
    }
    lepton_2_finalstate = {
        "ee": "epem",
        "mu": "mupmum"
    }
    target = exp_2_target[experiment]
    particle = lepton_2_finalstate[lepton]
    
    # Construct the hddm conversion command
    hddm_command = (
        f"python {ascii2hddm_script} {particle} {target} "
        f"{os.path.join(directories['vectorspath'], f'vectors_{job_i}.txt')} "
        f"{os.path.join(directories['vectorspath'], f'vectors_{job_i}.hddm')} "
        f"--run {hddm_run_number} "
        f"--vertex '0 0 0 0'"
    )
    
    # Path to the bash script
    bash_script = os.path.join(directories['template_directory'], RBHG_script)
    
    print(f"Running job {job_i} locally...")
    print(f"  Fortran exe: {leptonexe}")
    print(f"  Output dir: {directories['vectorspath']}")
    
    # Run the bash script with the same arguments as swif2 would
    command = [
        'bash', bash_script,
        ENVFILE,
        leptonexe,
        directories['vectorspath'],
        'python', ascii2hddm_script, particle, target,
        os.path.join(directories['vectorspath'], f'vectors_{job_i}.txt'),
        os.path.join(directories['vectorspath'], f'vectors_{job_i}.hddm'),
        '--run', str(hddm_run_number),
        '--vertex', '0 0 0 0'
    ]
    
    result = subprocess.run(command, capture_output=False)
    
    if result.returncode == 0:
        print(f"  ✓ Job {job_i} completed successfully")
    else:
        print(f"  ✗ Job {job_i} failed with return code {result.returncode}")
    
    return result.returncode


def append_unique_seed(seed_list):
    while True:
        # Generate 4 random bytes from OS-specific randomness source
        random_bytes = os.urandom(4)

        # Convert these bytes into an unsigned integer and ensure it's within the 32-bit signed integer range
        seed = struct.unpack("<I", random_bytes)[0] % 2147483648

        # Check if the generated seed is unique
        if seed not in seed_list:
            seed_list.append(seed)
            break  # Exit the loop once a unique seed is found

    return seed, seed_list  # Return the unique seed and the updated list



########################## MAIN #############################

RBHGfortranVer = determine_RBHG_fortran_ver(temp_fortran_file)
print(f"Fortran version = {RBHGfortranVer}")

list_of_seeds = []

if RUN_PERIOD:
    run_configurations = get_run_configurations(RUN_PERIOD)
else: 
    run_configurations = {"": {PolDeg+"DEG": nevents_total}}

# Determine what type of study this is
is_ffs = "FFS" in run_configurations
is_preset = RUN_PERIOD != ""  # Any non-empty run period is a preset
needs_run_period_master = is_preset  # Level 2 master for any preset
needs_ffs_master = is_ffs  # Level 3 master for FFS only

if is_ffs:
    # FFS mode: Generate 4 datasets (qDATAq, BC_FFN, BC_FF1, SIM_FF1)
    # Form factor designation moved to run_period level for cleaner naming
    ffs_datasets = [
        {"nametag": "qDATAq", "form_factor": "FFN", "tenx_flag": False},
        {"nametag": "BC", "form_factor": "FFN", "tenx_flag": True}, 
        {"nametag": "BC", "form_factor": "FF1", "tenx_flag": True},
        {"nametag": "SIM", "form_factor": "FF1", "tenx_flag": True}
    ]
    actual_run_configs = run_configurations["FFS"]  # Get the actual run/polarization data
    ffs_dataset_configs = {}  # Collect configs for Level 3 master
    run_period_master_configs = {}  # Collect configs for Level 2 masters (one per dataset)
else:
    # Normal mode: Single dataset with current settings
    ffs_datasets = [{"nametag": nametag, "form_factor": form_factor, "tenx_flag": TENX_FLAG}]
    actual_run_configs = run_configurations
    ffs_dataset_configs = None
    run_period_master_configs = {}  # Collect configs for Level 2 master (if preset)

# Loop through FFS datasets (or single dataset for normal mode)
for dataset in ffs_datasets:
    current_nametag = dataset["nametag"]
    current_form_factor = dataset["form_factor"] 
    current_tenx_flag = dataset["tenx_flag"]
    
    print(f"\n{'='*60}")
    print(f"Processing dataset: {current_nametag} (form_factor={current_form_factor}, 10x={current_tenx_flag})")
    print(f"{'='*60}")
    
    # Loop through run periods and polarizations
    for run_period, polarizations in actual_run_configs.items():
        for polarization, events in polarizations.items():
            print(f"Run Period: {run_period}, Polarization: {polarization}, Events: {events}")
            if RUN_PERIOD: 
                PolDeg = polarization.replace("DEG", "")
                print(polarization)
                nevents_total = events
                
                # Apply 10x multiplier if this dataset requires it
                if current_tenx_flag:
                    nevents_total = events * 10
                    print(f"10x flag enabled: {events} -> {nevents_total} events")

            # Automatic configuration for AMO (Amorphous/unpolarized) data
            current_BH_xsctn = BH_xsctn_formulation
            current_CobremsDistribution = CobremsDistribution
            if PolDeg == "AMO":
                current_BH_xsctn = "Heitler"    # AMO data should use unpolarized Heitler formulation
                current_CobremsDistribution = False  # AMO data should use 1/E distribution
                print("AMO polarization detected: automatically setting cross-section to Heitler and CobremsDistribution to False")

            logicals = {
            "RBHGHIST_W": Fortran_TorF(hist_invariantMass),
            "RBHGHIST_X": Fortran_TorF(hist_energyFraction_x),
            "RBHGHIST_T": Fortran_TorF(hist_mandelstam_t),
            "RBHGHIST_T_VARIABLE": Fortran_TorF(hist_t_varBinWidth),
            "RBHGHIST_ELASTICITY": Fortran_TorF(hist_elasticity),
            "RBHGHIST_MM": Fortran_TorF(hist_missingmass),
            "RBHGHIST_THETA": Fortran_TorF(hist_theta),
            "RBHGHIST_PHI_JT": Fortran_TorF(hist_phi_of_JT),
            "RBHGHIST_EGAMMA": Fortran_TorF(hist_Egamma),
            "RBHGOUTPUT_EVENT": Fortran_TorF(ascii_vector_output),
            "RBHGINTEGRAL_XSCTN": Fortran_TorF(integral_xsection),  
            "RBHGVERBOSE_OUTPUT": Fortran_TorF(verbose_output),     
            ####
            "RBHGBAKMAEV": Fortran_TorF(current_BH_xsctn == "Bakmaev"),  # B.H. cross section - Bakmaev's formulation. Do not use, unless you have a good reason to.
            "RBHGHEITLER": Fortran_TorF(current_BH_xsctn == "Heitler"),  # B.H. cross section - Unpolarized Heitler formulation          
            "RBHGBERLIN": Fortran_TorF(current_BH_xsctn == "Berlin"),    # B.H. cross section - Berlin polarized formulation. Should be selected as default.
            "RBHGPARA": Fortran_TorF(PolDeg in ["0", "45"] and PolDeg != "AMO"),  # Linear polarization at 0 (135) deg. False for 90 (45) and AMO
            "RBHGZERO_NINETY": Fortran_TorF(PolDeg in ["0", "90"] and PolDeg != "AMO"),  # Polarization angles of 0/90. False for 135/45 and AMO
            "RBHGCOBREMS": Fortran_TorF(current_CobremsDistribution),                # Use Coherent Bremsstrahlung distribution (true) or 1/E (false)
            "RBHGCOBREMS_VARBIN": Fortran_TorF(CobremsVarBin),                  # Use variable-width tagger binning for bremsstrahlung
            "RBHGGLUEX": Fortran_TorF(Experiment == "GlueX"),                # Use GlueX run conditions (true) or CPP run conditions (false)
            "RBHGELECTRON": Fortran_TorF(lepton == "ee"),                    # e+e- control. set .false. to generate BH muons
            "RBHGRADIATION": Fortran_TorF(internal_radiation),               # internal radiation control
            "RBHGSINGLE_RAD": Fortran_TorF(single_radiation),                # v113: single lepton radiation control
            "RBHGHYPGEOM": Fortran_TorF(hypgeom_radiation),                  # Future: hypergeometric radiation control
            "RBHGNUC_FF": Fortran_TorF(current_form_factor == "FFN")                 # True for FF = dipole form factor (FFN). False for FF = 1 (FF1).
        }

            job_i = 0 
            total_jobs = math.ceil(int(nevents_total)/int(nevents_perfile))

            # Create a name that accurately describes the generator's configuration for output directories
            gen_dir_name = make_generation_dirname(study_name, current_nametag, current_form_factor, polarization, run_period, lepton, current_BH_xsctn, nevents_total,
                                                  RUN_PERIOD, current_tenx_flag, Experiment, RBHGfortranVer, internal_radiation, single_radiation, hypgeom_radiation)
            print(f"Generation directory name = {gen_dir_name}")

            DEFAULTS = {"CPP": "2022-05", "GlueX": "2018-08"}
            hddm_run_number = get_run_number(Experiment, run_period, PolDeg, default_periods=DEFAULTS)
            
            # Create all relevant directories 
            alldirs, updatedlogicals = make_all_gen_directories(FrameworkHomeDirectory, gen_dir_name, CUEusername, logicals, run_period, polarization)

            # Create and save JSON configuration
            histogram_settings = {
                "invariant_mass": hist_invariantMass,
                "energy_fraction": hist_energyFraction_x,
                "mandelstam_t": hist_mandelstam_t,
                "t_variable_bins": hist_t_varBinWidth,
                "elasticity": hist_elasticity,
                "missing_mass": hist_missingmass,
                "theta": hist_theta,
                "phi_of_JT": hist_phi_of_JT,
                "egamma": hist_Egamma
            }
            
            output_settings = {
                "ascii_vector_output": ascii_vector_output,
                "integral_xsection": integral_xsection,
                "verbose_output": verbose_output
            }
            
            config = create_rbhg_config(
                study_name, current_nametag, current_form_factor, lepton, current_BH_xsctn,
                Experiment, PolDeg, nevents_total, nevents_perfile, internal_radiation,
                current_CobremsDistribution, RUN_PERIOD, current_tenx_flag, RBHGfortranVer,
                histogram_settings, output_settings, experiment_parameters, alldirs,
                single_radiation, hypgeom_radiation, run_period=run_period, polarization=polarization
            )
            
            config_file = write_rbhg_config(config, alldirs['output_directory'])
            
            # Collect config for Level 2 master (run period coordination) if needed
            if needs_run_period_master:
                dataset_key = f"{current_nametag}_{current_form_factor}"
                if dataset_key not in run_period_master_configs:
                    run_period_master_configs[dataset_key] = {}
                if run_period not in run_period_master_configs[dataset_key]:
                    run_period_master_configs[dataset_key][run_period] = {}
                run_period_master_configs[dataset_key][run_period][polarization] = config
            
            # Collect config for Level 3 master (FFS coordination) if needed
            if is_ffs:
                ffs_dataset_configs[current_nametag] = config

            workflow_name_base = make_workflow_name(study_name, current_nametag, current_form_factor, run_period, polarization, lepton, internal_radiation, single_radiation, hypgeom_radiation)
            workflow_name = f"gen_{workflow_name_base}"  # Add generation prefix
            
            if not INTERACTIVE_MODE:
                workflow_name = create_swif2_workflow(workflow_name)  # Capture the actual workflow name used

            # Generate job-specific parameters FIRST (we need them for compilation defaults and JSON)
            print(f"\n{'='*60}")
            print(f"Generating job parameters for {total_jobs} jobs...")
            print(f"{'='*60}")
            job_params = []
            job_i = 0
            
            while job_i < total_jobs:
                # handle the case where a user has requested a total number events not equally divisible by nevents_perfile:
                current_nevents_perfile = nevents_perfile
                if (job_i + 1) == total_jobs:
                    if (int(nevents_total) % int(nevents_perfile) != 0):
                        current_nevents_perfile = str(int(nevents_total) % int(nevents_perfile))
                
                # Generate unique initial seed
                latest_seed, updated_list_of_seeds = append_unique_seed(list_of_seeds)
                list_of_seeds = updated_list_of_seeds
                
                job_params.append({
                    'job_i': job_i,
                    'seed': latest_seed,
                    'nevents': current_nevents_perfile
                })
                job_i += 1
            
            print(f"Job parameters generated!\n")
            
            # Add job parameters to config for audit trail
            config['rbhg_config']['job_parameters'] = job_params
            
            # Update the config file with job parameters
            config_file = write_rbhg_config(config, alldirs['output_directory'])
            print(f"Job parameters saved to: {config_file}\n")

            # Compile ONCE - create single executable for this workflow
            print(f"{'='*60}")
            print(f"Compiling Fortran executable (once for all {total_jobs} jobs)...")
            print(f"{'='*60}")
            
            # Create one Fortran file with all workflow settings
            FortranFile = os.path.join(alldirs['submit_directory'], fortran_basename(temp_fortran_file) + ".f")
            command = ["cp", f"{alldirs['template_directory']}/{temp_fortran_file}", FortranFile]
            subprocess.run(command)
            
            # Do regex replacements for all non-job-specific settings
            updatedlogicals.update(experiment_parameters)
            # Use FIRST JOB's actual values as compile-time defaults (clear, not confusing dummy values)
            # These will be overridden by command-line arguments at runtime for jobs 1-N
            first_job = job_params[0]
            updatedlogicals['RBHGJOBNUMBER'] = str(first_job['job_i'])
            updatedlogicals['RBHGSEEDVALUE'] = str(first_job['seed'])
            updatedlogicals['RBHGNEVENTS'] = str(first_job['nevents'])
            
            length_sorted_logicals = sorted(updatedlogicals.items(), key=lambda x: len(x[0]), reverse=True)
            
            for placeholder, replacement in length_sorted_logicals:
                command = ['sed', '-i', f's|{placeholder}|{replacement}|g', FortranFile]
                subprocess.run(command)
            
            # Compile the single executable
            fortranexe = os.path.join(alldirs['submit_directory'], fortran_basename(temp_fortran_file) + ".exe")
            command = ['gfortran', '-std=legacy', '-ffixed-line-length-none', '-fd-lines-as-code', FortranFile, '-o', fortranexe]
            compile_result = subprocess.run(command, capture_output=True, text=True)
            
            if compile_result.returncode != 0:
                print(f"ERROR: Compilation failed!")
                print(compile_result.stderr)
                sys.exit(1)
            
            print(f"Compilation complete! Executable: {fortranexe}\n")
            
            # Submit all jobs using the single executable with different arguments
            if INTERACTIVE_MODE:
                print(f"{'='*60}")
                print(f"Running {total_jobs} jobs locally...")
                print(f"{'='*60}")
                for params in job_params:
                    print(f"Running job {params['job_i']}/{total_jobs}...")
                    run_job_locally_with_args(ENVFILE, fortranexe, alldirs, params['job_i'], params['seed'], 
                                             params['nevents'], lepton, Experiment, ascii2hddm_script, hddm_run_number, RBHG_script)
            else:
                print(f"{'='*60}")
                print(f"Submitting {total_jobs} jobs to swif2...")
                print(f"{'='*60}")
                for params in job_params:
                    if (params['job_i'] % 100 == 0) or (params['job_i'] == total_jobs - 1):
                        print(f"Submitting job {params['job_i']}/{total_jobs}...")
                    swif2_add_job_with_args(swif_addjob_dict_stub, workflow_name, ENVFILE, fortranexe, alldirs, 
                                           params['job_i'], params['seed'], params['nevents'], lepton, Experiment, 
                                           ascii2hddm_script, hddm_run_number)
                
                print(f"\nAll {total_jobs} jobs submitted to workflow '{workflow_name}'!")


# Create Level 2 master configurations (run period coordination) if needed
outputDirectoryTop = os.path.join(FrameworkHomeDirectory, "output", GENERATOR_TYPE)
run_period_master_files = []

if needs_run_period_master and run_period_master_configs:
    print(f"\n{'='*60}")
    print("Creating run period master configurations...")
    print(f"{'='*60}")
    
    for dataset_key, run_period_data in run_period_master_configs.items():
        dataset_name, form_factor = dataset_key.split('_', 1)
        
        print(f"Creating master config for {dataset_name}_{form_factor}...")
        master_config = create_run_period_master_config(
            study_name, dataset_name, form_factor, RUN_PERIOD, 
            run_period_data, outputDirectoryTop
        )
        master_file = write_run_period_master_config(
            master_config, study_name, dataset_name, form_factor, outputDirectoryTop
        )
        run_period_master_files.append(master_file)

# Create Level 3 master configuration (FFS coordination) if needed
if is_ffs and ffs_dataset_configs:
    print(f"\n{'='*60}")
    print("Creating Form Factor Study master configuration...")
    print(f"{'='*60}")
    
    master_config = create_ffs_master_config(study_name, RUN_PERIOD, ffs_dataset_configs, outputDirectoryTop)
    ffs_master_file = write_ffs_master_config(master_config, study_name, outputDirectoryTop)
    
    print(f"\nForm Factor Study setup complete!")
    print(f"Level 3 Master (FFS): {ffs_master_file}")
    print(f"Level 2 Masters (Run Periods): {len(run_period_master_files)} files")
    for file in run_period_master_files:
        print(f"  {file}")
    print(f"\nTo run prepareSimulation.py on all datasets:")
    print(f"  python prepareSimulation.py --master-config {study_name}/ffs_master_config.json")

elif needs_run_period_master and run_period_master_files:
    print(f"\nRun period coordination setup complete!")
    print(f"Level 2 Master configs: {len(run_period_master_files)} files")
    for file in run_period_master_files:
        print(f"  {file}")
    print(f"\nTo run prepareSimulation.py on this study:")
    for file in run_period_master_files:
        rel_path = file.replace(outputDirectoryTop + '/', '')
        print(f"  python prepareSimulation.py --run-period-config {rel_path}")





