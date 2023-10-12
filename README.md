# polyploid_expansion_load
SLiM simulator for polyploid expansion load study. All sims were run using SLiM v 4.0.1 https://messerlab.org/slim/

Once cloned, the scripts polyploid_expansion.slim and polyploid_single_population.slim scripts should run with the following commands as an example (if slim is installed and added to your path'. variables can be changed to alter parameters of the script. 

```
slim -d g_size=999999 -d K=100 -d "r=log(2)" -d mig_rate=0.005 -d u_del=2.5e-8 -d u_ben=2.5e-9 -d b_s=0.000 -d d_s=-0.005 -d rho=2.5e-8 -d "inheritance='auto'" -d "dom_pattern='recessive'" polyploid_expansion.slim 
```
or 

```
slim -d g_size=999999 -d K=500 -d "r=log(2)" -d mig_rate=0.005 -d u_del=2.5e-8 -d u_ben=2.5e-9 -d b_s=0.000 -d d_s=-0.005 -d rho=2.5e-8 -d "inheritance='auto'" -d "dom_pattern='recessive'" polyploid_single_population.slim  
```

For the diploidization scripts, if using the dominant meiotic model run the following script:

```
slim -d "out_dir='output/'" -d g_size=999999 -d K=100 -d "r=log(2)" -d dip_lambda=2 -d rho=2.5e-8 -d remove_dip_muts=1 -d mig_rate=0.05 -d u_del=2.5e-8 -d u_ben=2.5e-9 -d u_dip=2.5e-10 -d b_s=0.00 -d d_s=-0.005 -d "dom_pattern='auto'" -d "dip_model='dom'" polyploid_diploidization_dominant_meiotic.slim
```

and for the pairing efficiency model:

```
slim -d "out_dir='output/'" -d g_size=999999 -d K=100 -d "r=log(2)" -d dip_lambda=100 -d rho=2.5e-8 -d remove_dip_muts=1 -d meiotic_fitness=1 -d mig_rate=0.05 -d u_del=2.5e-8 -d u_ben=2.5e-9 -d u_dip=2.5e-10 -d b_s=0.00 -d d_s=-0.005 -d "dom_pattern='auto'" -d "dip_model='diff'" polyploid_diploidization_pairing_efficiency.slim
```
for diploidization, the extra parameter 'remove_dip_muts' governs whether diploidization mutations are removed prior to expansion (1 is yes, 0 no), and for the pairing efficiency model, 'meiotic fitness' regulates if pairing associate fitness effects are included (PAFCs; 1 for yes, 0 for no) 
