#gxenv /group/halld/www/halldweb/html/halld_versions/version_root_python3.6.xml
import os
import shutil
import glob
import numpy as np
import ROOT
from array import array

def find_prefixes(hists_path):
    files = os.listdir(hists_path)
    # Exclude files that end with '_merged.txt' in case you need to run over the directory twice.
    prefixes = set(os.path.splitext(name)[0].rsplit('_', 1)[0] for name in files if name.startswith('array') and not name.endswith('_merged.txt'))
    return list(prefixes)

def merge_histograms(base_directory, cleanup_mode=None, cleanup_folder=None, output_filename='merged_histograms.root'):
    hists_path = os.path.join(base_directory, 'hists')
    prefixes = find_prefixes(hists_path)

    # Create or overwrite the ROOT file
    root_file = ROOT.TFile(os.path.join(hists_path, output_filename), 'RECREATE')

    for prefix in prefixes:
        if "array_tvar" in prefix:
            continue

        files = glob.glob(os.path.join(hists_path, '{}_*.txt'.format(prefix)))
        if not files:
            print("No files found for prefix: {}".format(prefix))
            continue

        # Initialize bins and counts from the first file
        first_file_data = np.loadtxt(files[0], usecols=(0, 1))
        bins = first_file_data[:, 0]
        counts = first_file_data[:, 1]

        # Add counts from the remaining files
        for f in files[1:]:
            data = np.loadtxt(f, usecols=(0, 1))
            counts += data[:, 1]

        # Create the bin edges
        bin_edges = np.concatenate(([bins[0] - (bins[1] - bins[0]) / 2], bins[:-1] + np.diff(bins) / 2, [bins[-1] + (bins[-1] - bins[-2]) / 2]))
        bin_edges = array('d', bin_edges)  # Convert to array for ROOT

        # Create and fill the histogram
        hist = ROOT.TH1D(prefix, prefix, len(bin_edges) - 1, bin_edges)
        for i, count in enumerate(counts, 1):
            hist.SetBinContent(i, count)

        # Write the histogram to the ROOT file
        hist.Write()


        # Write merged histogram data to a text file
        errors = np.sqrt(counts)
        merged_data = np.column_stack((bins, counts, errors))
        merged_file_path = os.path.join(hists_path, '{}_merged.txt'.format(prefix))
        np.savetxt(merged_file_path, merged_data, fmt='%f\t%f\t%f')

    
    # Create variable-width histogram
    var_bin_edges = array('d', [
        0.0000, 0.0002, 0.0006, 0.0013, 0.002, 0.003, 0.0042, 0.0055, 0.007, 0.009, 0.011,
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
        0.365, 0.367, 0.369, 0.371, 0.373, 0.375, 0.377, 0.379, 0.381, 0.383, 0.385, 0.387, 0.389, 0.391])
    
    var_width_hist = ROOT.TH1D("array_tvar", ";-t (GeV/c)^2", len(var_bin_edges) - 1, var_bin_edges)

    # Initialize arrays to store bin centers and widths
    bin_centers = []
    bin_widths = []


    tvar_files = glob.glob(os.path.join(hists_path, 'array_tvar_*.txt'))
    for file_path in tvar_files:
        data = np.loadtxt(file_path)
        if len(bin_centers) == 0:  # If empty, initialize bin_centers and bin_widths
            bin_centers = data[:, 0]
            bin_widths = data[:, 1]

        counts = data[:, 2]
        for i, count in enumerate(counts, 1):
            var_width_hist.AddBinContent(i, count)

    var_width_hist.Write()

    # Write merged histogram data for variable width bins to a text file
    merged_counts = np.array([var_width_hist.GetBinContent(i+1) for i in range(var_width_hist.GetNbinsX())])
    merged_data = np.column_stack((bin_centers, bin_widths, merged_counts))
    merged_file_path = os.path.join(hists_path, 'array_tvar_merged.txt')
    np.savetxt(merged_file_path, merged_data, fmt='%f\t%f\t%f')

    # Close the ROOT file
    root_file.Close()

    if cleanup_mode in ['move', 'delete']:
        for prefix in prefixes:
            files = glob.glob(os.path.join(hists_path, '{}_*.txt'.format(prefix)))
            for file_path in files:
                if file_path.endswith('_merged.txt') or file_path.endswith('.root'):
                    continue
                if cleanup_mode == 'move':
                    if not cleanup_folder:
                        cleanup_folder = os.path.join(hists_path, 'archived_histograms')
                    if not os.path.exists(cleanup_folder):
                        os.makedirs(cleanup_folder)
                    shutil.move(file_path, os.path.join(cleanup_folder, os.path.basename(file_path)))
                elif cleanup_mode == 'delete':
                    os.remove(file_path)

# End definition of functions
##################################################################################################################

# Example usage (only runs when script is executed directly, not when imported)
if __name__ == "__main__":
    #directories = [
    #"/work/halld/home/acschick/channels/epemmissprot/newRBHG/output/TS6_11E_qDATAq_1801_MissingEvts_FFN_0DEG_epem_GlueX_v108/",
    #"/work/halld/home/acschick/channels/epemmissprot/newRBHG/output/TS6_11E_qDATAq_1801_MissingEvts_FFN_45DEG_epem_GlueX_v108/",
    #"/work/halld/home/acschick/channels/epemmissprot/newRBHG/output/TS6_11E_qDATAq_1808_MissingEvts_FFN_45DEG_epem_GlueX_v108/",
    #]


    directories = [
    "/work/halld/home/acschick/channels/epemmissprot/newRBHG/output/TS6_11E_SIMFF1_FF1_FULL2018_10x/TS6_11E_SIMFF1_1801_10x_FF1_0DEG_epem_GlueX_v108",
    "/work/halld/home/acschick/channels/epemmissprot/newRBHG/output/TS6_11E_SIMFF1_FF1_FULL2018_10x/TS6_11E_SIMFF1_1801_10x_FF1_135DEG_epem_GlueX_v108",
    "/work/halld/home/acschick/channels/epemmissprot/newRBHG/output/TS6_11E_SIMFF1_FF1_FULL2018_10x/TS6_11E_SIMFF1_1801_10x_FF1_45DEG_epem_GlueX_v108",
    "/work/halld/home/acschick/channels/epemmissprot/newRBHG/output/TS6_11E_SIMFF1_FF1_FULL2018_10x/TS6_11E_SIMFF1_1801_10x_FF1_90DEG_epem_GlueX_v108",
    "/work/halld/home/acschick/channels/epemmissprot/newRBHG/output/TS6_11E_SIMFF1_FF1_FULL2018_10x/TS6_11E_SIMFF1_1808_10x_FF1_0DEG_epem_GlueX_v108",
    "/work/halld/home/acschick/channels/epemmissprot/newRBHG/output/TS6_11E_SIMFF1_FF1_FULL2018_10x/TS6_11E_SIMFF1_1808_10x_FF1_135DEG_epem_GlueX_v108",
    "/work/halld/home/acschick/channels/epemmissprot/newRBHG/output/TS6_11E_SIMFF1_FF1_FULL2018_10x/TS6_11E_SIMFF1_1808_10x_FF1_45DEG_epem_GlueX_v108",
    "/work/halld/home/acschick/channels/epemmissprot/newRBHG/output/TS6_11E_SIMFF1_FF1_FULL2018_10x/TS6_11E_SIMFF1_1808_10x_FF1_90DEG_epem_GlueX_v108"
    ]



    for directory in directories:
        merge_histograms(directory, cleanup_mode='move')


    # to keep the small text files unmerged, and move them to another folder away from the merged use:
    # merge_histograms(base_directory, cleanup_mode='move', cleanup_folder='/path/to/archive/folder')
    # Or to delete:
    # merge_histograms(base_directory, cleanup_mode='delete')
