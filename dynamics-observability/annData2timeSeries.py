"""
This script converts datasets from an AnnData object to a format compatible with building models of time series data. AnnData from a .h5ad file are converted to a .pkl file that will contain the state type or coordinates space (raw data, geneformer embedinggs, etc.), the state names, and the time series data organized as states (coordinates) by timepoint (hours) by replicates (experiment replicates).

Arguments (three arguments):
    1. output directory: chatstatus['output-directory']
    2. output file: <name of output file>
    3. input file: <file created in previous step>

Usage:
    Command Line:
    ```
    python <path/to/script/>annData2timeSeries.py <output path> <output file> <input file>
    ```
                                                    |              |           |
                                                Argument 1     Argument 2   Argument 3
    BRAD Line:
    ```
    subprocess.call([sys.executable, '<path/to/script/>/annData2timeSeries.py', chatstatus['output-directory'], <output file>, <input file>])
    ```

**OUTPUT FILE NAME INSTRUCTIONS**
1. Output path should be chatstatus['output-directory']
2. Output file name should be `<descriptive name>.pkl`
"""
import os  # All BRAD tools
import sys # All BRAD tools

import anndata as an
import numpy as np
import pickle


def main():
    outputPath = sys.argv[1] # chatstatus['output-directory']
    outputFile = sys.argv[2] # output file name
    dataset    = sys.argv[3] # input file

    print('******************************')
    print('  Convert to Gene Coordinates ')
    print('******************************')
    
    print('Dataset=', end='')
    print(dataset)

    ad = an.read(dataset)
    
    # Build the object in which to save the trajectories. This will be state variables by time by trajectory
    n = ad.shape[1]
    T = ad.obs['order'].nunique()
    R = ad.obs['replicate'].nunique()
    X = np.zeros((n, T, R))

    # List possible timepoints & replicates
    timepoints = ad.obs['order'].unique()
    replicates = ad.obs['replicate'].unique()

    # Loop over dataframe
    for i in range(ad.obs.shape[0]):
        order = ad.obs['order'][i]
        replicate = ad.obs['replicate'][i]
    
        time = np.where(timepoints == order)[0][0]
        rep = np.where(replicates == replicate)[0][0]
    
        X[:, time, rep] = ad.X[i, :]
    states = list(ad.var.index)
    metadata = ad.var
    
    output_file = os.path.join(outputPath, outputFile)
    with open(output_file, 'wb') as f:
        pickle.dump({'X': X, 'states': states, 'metadata': metadata}, f)
    
    print(f"Data saved to {output_file}")


if __name__ == "__main__":
    main()
