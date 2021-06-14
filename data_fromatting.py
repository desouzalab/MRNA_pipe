import subprocess
import sys
def install(package):
    subprocess.check_call([sys.executable, "-m", "pip", "install", package])

install("pandas")
install("numpy")
import pandas as pd
import numpy as np




p1 = "/home/emiliano/projects/def-cdesouza/Lab/data/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt"

p2 = "/home/emiliano/projects/def-cdesouza/Lab/data/GSE97930_CerebellarHem_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt"

p3 = "/home/emiliano/projects/def-cdesouza/Lab/data/GSE97930_VisualCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt"


data=pd.read_csv(p1,delimiter = "\t")
data2=pd.read_csv(p2,delimiter = "\t")
data3=pd.read_csv(p3,delimiter = "\t")

full=pd.concat([data,data2,data3],axis = 1)
end_p = "/home/emiliano/projects/def-cdesouza/Lab/data/raw"
full = full.fillna(0)
full.to_csv(end_p)