"""
This script is a file chooser to load a gene expression dataset from the lab's turbo partition to give access to the data to BRAD.

Arguments (three arguments):
    1. output directory: chatstatus['output-directory']
    2. output file: <name of output file>
    3. dataset: this can be any of the following:
        - "2015": a bulk RNAseq time series dataset of synchronized Fibroblast proliferation
        - "2018": a bulk RNAseq time series dataset of Weintraubs Myogenic reprogramming experiment

Usage:
    Command Line:
    ```
    python <path/to/script/>datasetSelector.py <output path> <output file> <gene name>
    ```
                                                               |              |           |
                                                          Argument 1     Argument 2   Argument 3
    BRAD Line:
    ```
    subprocess.call([sys.executable, '<path/to/script/>/datasetSelector.py', chatstatus['output-directory'], <output file>, <gene name>])
    ```

**OUTPUT FILE NAME INSTRUCTIONS**
1. Output path should be chatstatus['output-directory']
2. Output file name should be `S0-<data set name>.h5ad`
"""
import os
import sys
import pandas as pd
import anndata as an

def main():
    outputPath = sys.argv[1] # chatstatus['output-directory']
    outputFile = sys.argv[2] # output file name
    dataset    = sys.argv[3]

    outputFile = os.path.join(outputPath, outputFile)
    
    print('******************************')
    print('       Dataset Selector       ')
    print('******************************')

    out_path = "/nfs/turbo/umms-indikar/shared/projects/geneformer/data/rajapakse_lab_data.h5ad"
    ad = an.read(out_path)
    
    if dataset == '2015':
        ds = 'chen_2015'
    elif dataset == '2018':
        ds = 'liu_2018'
    else:
        print('Dataset is invalid!')

    adDs = ad[ad.obs['dataset'] == ds]
    adDs.write(outputFile)
    

if __name__ == "__main__":
    main()
