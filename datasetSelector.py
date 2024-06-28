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
2. Output file name should be `S0-<data set name>.csv`
"""
import os
import sys
import pandas as pd

def main():
    outputPath = sys.argv[1] # chatstatus['output-directory']
    outputFile = sys.argv[2] # output file name
    dataset    = sys.argv[3]

    outputFile = os.path.join(outputPath, outputFile)
    
    print('******************************')
    print('Dataset Selector')
    print('******************************')

    if dataset == '2015':
        filepath = '/nfs/turbo/umms-indikar/Joshua/bioObsv/notebooks/obsvArticle/dataProcessing/2015_avg.csv'
    elif dataset == '2018':
        filepath = '/nfs/turbo/umms-indikar/Joshua/bioObsv/notebooks/obsvArticle/dataProcessing/2018_avg.csv'
    else:
        print('Dataset is invalid!')

    df = pd.read_csv(filepath)
    df.to_csv(outputFile, index=False)

if __name__ == "__main__":
    main()
