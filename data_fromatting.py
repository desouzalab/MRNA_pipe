#!/usr/bin/python3
## file: data_formatting.py
import pandas as pd
import numpy as np
p1 = input("path to first file")
p2 = input("path to second file")
p3 = input("path to third file")

data=pd.read_csv(p1,delimiter = "\t")
data2=pd.read_csv(p2,delimiter = "\t")
data3=pd.read_csv(p3,delimiter = "\t")

full=pd.concat([data,data2,data3],axis = 1)
end_p = input("path to output")
full = full.fillna(0)
full.to_csv(end_p)