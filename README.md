# The genetic consequences of expansion and its influence on diploidization in polyploids
SLiM simulator for polyploid expansion load study. All sims were run using SLiM v 4.0.1 https://messerlab.org/slim/

Once cloned, the scripts polyploid_expansion.slim and polyploid_single_population.slim scripts should run with the following commands as an example (if slim is installed and added to your path). variables can be changed to alter parameters of the script. 

For keyword parameters, the following are accepted:

| inheritance         |                                                                  |
|:-------------------|:----------------------------------------------------------------------------|
| __auto__          |  autotetraploid                                       |
| __allo__               |  allotetraploid         |
| __diploid__                |  diploid    |

| dom_pattern         |  h_vector                                                                |
|:-------------------|:----------------------------------------------------------------------------|
| __additive__          |  beneficial: (0.0, 0.25, 0.5, 0.75, 1.0), deleterious: (0.0, 0.25, 0.5, 0.75, 1.0)                                      |
| __recessive__               |  beneficial: (0.0, 0.0, 0.0, 0.0, 1.0), deleterious: beneficial: (0.0, 0.0, 0.0, 0.0, 1.0)|
| __bd_dr__                |  beneficial: (0.0, 1.0, 1.0, 1.0, 1.0), deleterious (0.0, 0.0, 0.0, 0.0, 1.0)    |
| __br_dd__          |  beneficial: (0.0, 0.0, 0.0, 0.0, 1.0), deleterious (0.0, 1.0, 1.0, 1.0, 1.0)                                      |
| __duplex__               |  beneficial: (0.0, 1.0, 1.0, 1.0, 1.0), deleterious (0.0, 0.0, 0.0, 1.0, 1.0)         |
| __DFE__               |  beneficial: (0.0, 1.0, 1.0, 1.0, 1.0), deleterious (estimated h-s relationship)         |

NOTE: for the paper, bd_dr is used for the recessive model as deleterious are recessive and beneficial are dominant


| s_dist         |                                                                  |
|:-------------------|:----------------------------------------------------------------------------|
| __fixed__          |  both beneficial and deleterious set to b_s and d_s, respectively                                       |
| __exp__               |  deleterious drawn from an exponential distribution with mean d_s, beneficial fixed at b_s         |
| __gamma__                |  deleterious drawn from a gamma distribution at mean -0.001472, beneficial fixed at b_s   |


To run each of the scripts, the following will work as examples 
```
slim -d "out_dir='out'" -d g_size=999999 -d K=100 -d "r=log(2)" -d mig_rate=0.05 -d u_del=2.5e-8 -d u_ben=2.5e-9 -d b_s=0.005 -d d_s=-0.005 -d rho=2.5e-8 -d "inheritance='auto'" -d "dom_pattern='bd_dr'" -d "s_dist='gamma'" polyploid_expansion.slim 

```

For the diploidization scripts, if using the dominant meiotic model run the following script:

```
slim -d "out_dir='out'" -d g_size=999999 -d K=100 -d "r=log(2)" -d dip_lambda=2 -d rho=2.5e-8 -d remove_dip_muts=1 -d mig_rate=0.05 -d u_del=2.5e-8 -d u_ben=2.5e-9 -d u_dip=2.5e-10 -d b_s=0.05 -d d_s=-0.005 -d "dom_pattern='auto'" -d "dip_model='dom'" polyploid_diploidization_dominant_meiotic.slim
```

and for the pairing efficiency model:

```
slim -d "out_dir='out'" -d g_size=999999 -d K=100 -d "r=log(2)" -d dip_lambda=100 -d rho=2.5e-8 -d remove_dip_muts=1 -d meiotic_fitness=1 -d mig_rate=0.05 -d u_del=2.5e-8 -d u_ben=2.5e-9 -d u_dip=2e-3 -d b_s=0.005 -d d_s=-0.005 -d "dom_pattern='bd_dr'" -d "dip_model='diff'" -d pe_inflection=85 -d pe_slope=1 -d "s_dist='gamma'" polyploid_diploidization_pairing_efficiency.slim
```
for diploidization, the extra parameter 'remove_dip_muts' governs whether diploidization mutations are removed prior to expansion (1 is yes, 0 no), and for the pairing efficiency model, 'meiotic fitness' regulates if pairing associate fitness effects are included (PAFCs; 1 for yes, 0 for no), pe_inflection is the inflection point of the model, and pe_slope is the slope parameter of the model.


For a single population, an example script has been added to the slim_extras repository https://github.com/MesserLab/SLiM-Extras/blob/master/models/autotetraploid.slim