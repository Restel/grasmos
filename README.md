The following folder contains source code of GRASMOS: Graph Signage Model Selection for Gene Regulatory Networks:

[1] /data/

Syntetic and real datasets (when downloaded through data_load.R function) used for analysis

[2] data_load.R and data_clean.R scripts

[3] grasmos_SC_TC_explore.py 

Function to run the exploration of likelihoods of theta belonging to SC, TC cases on the specified input

[4] Functions to run the coarse and fine-grid exploration of theta belonging to BNC cases on Regulon and SubtiWiki datasets

### Usage

#### Real GRN data loading and cleaning
```
>>> Rscript data_load.R # load RegulonDB and SubtiWiki
trying URL 'https://regulondb.ccg.unam.mx/menu/download/datasets/files/NetWorkTFGene.txt'
Content type 'text/plain' length 300283 bytes (293 KB)
==================================================
downloaded 293 KB

trying URL 'http://subtiwiki.uni-goettingen.de/v4/regulation/exporter'
downloaded 308 KB
>>> Rscript data_clean.R 
```

#### Create tested theta parameters (BNC case)

```
>>> python3 create_params.py # saves all used parameters for into a local folder parameters/
```

#### BNC likelihood estimation

IF LOCAL COMPUTATION

```
>>> /bin/bash color_models_local.sh # a function to launch the estimation of ALL BNC theta parameters from /parameter/ folder locally on PC
```

IF PARALLEL COMPUTATION
```
>>> /bin/bash color_models_RC.sh 
```
