"""
Identify highly variable genes in single-cell RNA-seq data and save the results

BRAD's Command to Run this file:
subprocess.run([sys.executable, '<path/to/script/>/highly_variable_genes.py', chatstatus['output-directory'], <input file>, <output file>, <min mean>, <max mean>, <min dispersion>], capture_output=True, text=True)

Arguments:
  1. output directory: chatstatus['output-directory']
  2. input file: <.h5ad file that can be read by scanpy>
  3. output file: <name of output file>
  4. min mean: <minimum mean for identifying highly variable genes>
  5. max mean: <maximum mean for identifying highly variable genes>
  6. min dispersion: <minimum dispersion for identifying highly variable genes>

**OUTPUT FILE NAME INSTRUCTIONS**
1. Output path should be chatstatus['output-directory']
2. Output file name should be `<descriptive name>.h5ad`
"""

# Auth: Joshua Pickard
#       jpic@umich.edu
# Date: July 22, 2024

import os
import sys
import scanpy as sc

def main():
    outputPath = sys.argv[1]  # chatstatus['output-directory']
    inputFile = sys.argv[2]   # input file
    outputFile = sys.argv[3]  # output file name
    min_mean = float(sys.argv[4])  # minimum mean
    max_mean = float(sys.argv[5])  # maximum mean
    min_disp = float(sys.argv[6])  # minimum dispersion

    # Load the input file
    adata = sc.read_h5ad(inputFile)
    
    # Identify highly variable genes
    sc.pp.highly_variable_genes(adata, min_mean=min_mean, max_mean=max_mean, min_disp=min_disp)
    
    # Save the processed AnnData object to a file
    output_path_full = os.path.join(outputPath, outputFile)
    adata.write(output_path_full)
    print(f'Data written to {output_path_full}')

