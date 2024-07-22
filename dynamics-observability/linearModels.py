"""
This file builds linear and linear time variant (LTI and LTV) models of gene expressions.

Arguments (four required arguments, one optional argument):
    1. output path: chatstatus['output-directory']
    2. output file: <name of output file>
    3. input file: <name of the input file with data to build models from>
    4. model parameters:
        - type: <DMD or DGC>
        - DMD Rank = 0 through 20. This is only used for DMD. Do not use this with DGC.

Dynamic Mode Decomposition (DMD) facilitates model reduction. Lower ranks enhance model reduction, while higher ranks improve model accuracy. However, excessively high ranks increase the risk of overfitting.

Data Guided Control (DGC) uses linear time variant models to fit the data. This requires averaging over multiple replicates of the experiment.

Usage:
    Command Line:
    ```
    python <path/to/script/>geneformerEmbedding.py <output path> <output file> <input file> <model type> <DMD rank>
    ```
                                                       |              |            |             |            |
                                                   Argument 1     Argument 2   Argument 3    Argument 4   Argument 5
    BRAD Line:
    ```
    subprocess.call([sys.executable, '<path/to/script/>/geneformerEmbedding.py', chatstatus['output-directory'], <output file>, <input file>])
    ```

Output file name should be `S3-Model.pkl`

This script uses Geneformer to embed gene expression data.
"""