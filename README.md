# SME
spatial meta-ecosystem study: code for reproducing the results

For running the scripts, you need to create a /results/images folder in your directory, and also have /data/parameter_set.csv

The main script is (surprisingly)

### SME_main_script

it generates all files for generating the figures of the manuscript. For replicating the manuscript, run it three times with parameter configurations 6,7, and 11.

Then, the figures are generated with

### SME_plot_results

and the one related to network structure with

### SME_network_metrics_effect_type

the sensitivity analyses of the appendix need additional simulations. Check the main script adapted to it, and the script for generating its results:

### SME_main_script_sensitivity
### SME_plot_sensitivity_results