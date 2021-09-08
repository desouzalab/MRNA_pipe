v# MRNA_pipe

## This Repository belongs to the RNA Sequence Pipeline for DR.Camila de Souza's Laboratory 

<br>

## Installation
To make use of this Pipeline and it's functionalities please install the latest version of [Anaconda](https://www.anaconda.com/)

Set up Anaconda Enviroment 
`conda install -c conda-forge mamba`


`mamba create -c conda-forge -c bioconda -n rnaPipe snakemake `

`conda activate rnaPipe`

Install Backspin 

`conda install -c bioconda backspinpy`

Install [Giniclust](https://github.com/lanjiangboston/GiniClust) and edit the path in /rnaPipe/giniclust.R to the corresponding Giniclust path 


## Execution 

run `cd /Snakemake/<PROCESS TO RUN>`

Then execute the desired process with 

`snakemake --configfile config.yaml -c1 `
