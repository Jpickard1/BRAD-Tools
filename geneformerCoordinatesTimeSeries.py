"""
This script organizes datasets into a format to work with time series gene expression using individual geneformer embeddings as the state space. From a .h5ad file, a .pkl file is constructed here that will contain the state type being geneformer, the state names being geneformer coordinates, and the time series data organized as states (geneformer coordinates) by timepoint (hours) by replicates (experiment replicates).

Arguments (three arguments):
    1. output directory: chatstatus['output-directory']
    2. output file: <name of output file>
    3. input file: <file created in previous step>

Usage:
    Command Line:
    ```
    python <path/to/script/>datasetSelector.py <output path> <output file> <input file>
    ```
                                                    |              |           |
                                                Argument 1     Argument 2   Argument 3
    BRAD Line:
    ```
    subprocess.call([sys.executable, '<path/to/script/>/geneformerCoordinateTimeSeries.py', chatstatus['output-directory'], <output file>, <input file>])
    ```

**OUTPUT FILE NAME INSTRUCTIONS**
1. Output path should be chatstatus['output-directory']
2. Output file name should be `S1-<data set name>.h5ad`
"""