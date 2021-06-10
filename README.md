# MRNA_pipe

## This Repository belongs to the RNA Sequence Pipeline for DR.Camila de Souza's Laboratory 

### Pipeline Structure 
<br>

![Alt text](https://github.com/desouzalab/MRNA_pipe/blob/Master/docs/Pipeline_Diagram.png "Title")

Task can be executed independantlly with 4 main tasks:

* Preprocess
* Clustering
* Visualization
* Comparison

When processing new data pipeline should be executed in sequential order, but when tuning of a specitfic task is desired without rerunning other procedures, specific tasks can be executed in isolation.
<br>
### Repo Structure

<br>

![Alt text](docs\Repo_structure.png "Title")

<br>

The Repo follows a Classic version control structure 

<br>

Feature branches stem from Master and are merged back to master when all testing is complete.
This allows for each featur eto be tested independantly of each other in order to not change any existing output.

<br>

### Execution 
<br>

Pipeline is Executed from Jenkins to have a history log of executed procedures:

<br>

![Alt text](docs\rna_jenkins.PNG "Title")

This task calls [main.sh](main.sh) using ssh to connect to sharcnet.

It is parameterized in order to be able to execute various tasks when desired as well as to specify desired method being worked on.
