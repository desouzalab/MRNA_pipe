import requests

URL = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2612539'
page = requests.get(URL)
print(page)
