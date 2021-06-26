import requests
import re
import pandas as pd
import numpy as np
from bs4 import BeautifulSoup

############################## WAS USED TO SCRAPE AND FORMAT DATA FOR CELL CLUSTERS ##############################

column_list = []
target_list = []
not_found = []
for i in df["data.GSM.ID"]:
    try:
        URL = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={0}'.format(i)
        page = requests.get(URL)

        soup = BeautifulSoup(page.content, 'html.parser')
        web= str(page.content) 
        target = re.search("(?<=processed data column:).*",web)
        target_s=str(target.group(0) )
        column=target_s.split("<")[0]
        column_list.append(column)
        target_list.append(i)
     
    except AttributeError:
        not_found.append(i)
        
matches=pd.DataFrame({"GSM_ID":target_list,"label":column_list})
matches.to_csv(r"C:\Users\Emiliano\Documents\1_Lab_DIR_2020\data\labels_GSE98816.csv",index=False)
data = pd.read_csv(r"C:\Users\Emiliano\Documents\1_Lab_DIR_2020\data\GSE98816.csv",index_col = 0)