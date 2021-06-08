#!/usr/bin/python3
## file: all_jobs.py


# Load Python
module load python/3.6

# Preprocess
python3 ../../rnaSeq/code/preprocess/submit_preprocess_jobs.py

# Cluster
# OLD: python3 rnaSeq/code/cluster/sc3/submit_sc3_jobs.py
# OLD: python3 rnaSeq/code/cluster/seurat/submit_seurat_jobs.py

python3 rnaSeq/code/cluster/submit_cluster_jobs.py

# Visualize
python3 rnaSeq/code/visualize/submit_visualize_jobs.py

# Compare
python3 rnaSeq/code/compare/submit_compare_jobs.py

print("Done submitting all jobs!")
