"""
DO NOT RUN THIS CODE.

"""

# ts_embedding_1.py
import os
import sys
import scanpy as sc
import anndata as ad

sys.path.append('../anndata2embedding')
from embed import *

def main():
    job_number = int(sys.argv[1])
    print('Start')
    print('tGPT embeddings')
    print(job_number)
    
    DATAPATH = "/nfs/turbo/umms-indikar/shared/projects/DGC/data/tabula_sapiens/extract/"
    OUTPUTPATH = '/nfs/turbo/umms-indikar/shared/projects/foundation_models/experiments/tabulaSapiens'
    h5ad_files = []
    for root, dirs, files in os.walk(DATAPATH):
        for file in files:
            if file.endswith('.h5ad'):
                h5ad_files.append(os.path.join(root, file))

    file = h5ad_files[job_number]
    print(f"Reading file: {file}")
    adata = sc.read_h5ad(file)
    print(f"adata shape: {adata.shape}")
    adata.var['ensembl_id'] = adata.var['ensemblid']
    adata.var['ensembl_id'] = adata.var['ensembl_id'].str.split('.').str[0]

    embedAdTgpt = embed(adata, 'tGPT',
                      output_directory = OUTPUTPATH,
                      output_file      = os.path.splitext(os.path.basename(file))[0]
                     )
    embedAdTgpt.write(os.path.join(OUTPUTPATH, os.path.splitext(os.path.basename(file))[0] + '_Tgpt.h5ad'))
    print('Job well done')
    sys.exit(0)


if __name__ == "__main__":
    main()
