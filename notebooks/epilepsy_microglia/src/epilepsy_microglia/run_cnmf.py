import sys
from cnmf import cNMF

worker_i = int(sys.argv[1])
total_workers = int(sys.argv[2])

cnmf_obj = cNMF(
    output_dir="/data/SJLEE/Epilepsy_Microglia/notebooks/epilepsy_microglia/data/cnmf_results/snRNA_microglia_run3_filtered",
    name="snRNA_microglia"
)

cnmf_obj.factorize(
    worker_i=worker_i,
    total_workers=total_workers
)