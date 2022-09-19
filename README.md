# GRASMOS: Graph Signage Model Selection for Gene Regulatory Networks

[1] /data/

Real and syntetic datasets used for analysis

[2] /data_source/

Containts cleaning scripts for real datasets, and their original uncleaned sources

[3] grasmos_SC_TC_explore.py - Function to run the exploration of likelihoods of theta belonging to SC, TC cases on the specified input

[4] functions to run the coarse and fine-grid exploration of theta belonging to BNC cases on Regulon and SubtiWiki datasets

## Usage 
Create parameters
```
python3 create_params.py # saves all used parameters for into a local folder parameters/
```  
**If local computation**

```
/bin/bash color_models_local.sh  # a function to launch the estimation of ALL parameters from /parameter/ folder locally on PC
```
**If parallel computation**

```
/bin/bash color_models_RC.sh # the syntax is specific to our Research computing cluster, so you might need to change that
```
