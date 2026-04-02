# NIC EEG Clean Reproduction Repository

This repository contains the final clean scripts for the NIC EEG manuscript mainline analysis.

The repository is designed for clean reproduction of the manuscript-linked workflow only.  
Original scripts, manuscript files, review materials, and local result outputs are not included.

## Included contents

- final clean scripts
- usage instructions

## Not included

- original scripts
- manuscript files
- reviewer files
- local output/result files

## Canonical EEG network output folders

Under the network output directory, use the following four folder names:

- `NeuroEPO_post`
- `NeuroEPO_pre`
- `Placebo_post`
- `Placebo_pre`

## Mainline analysis chain

EEG network outputs  
→ signal activity  
→ topz nodes  
→ topz rates  
→ GLMM  
→ mediation  
→ figures

## Scripts

### EEG network
- `EEG_network_NeuroEPO_post_clean_main_english_prompt.R`
- `EEG_network_NeuroEPO_pre_clean_main_english_prompt.R`
- `EEG_network_Placebo_post_clean_main_english_prompt.R`
- `EEG_network_Placebo_pre_clean_main_english_prompt.R`

### TopZ + GLMM
- `05_topz_glmm_main_clean_final_new_EEG_paths.R`

### Mediation
- `hub_mediation_20230401_clean_TOPZ_only.R`

### Figures
- `Figure_2_pathfixed_current_pc_v3.R`
- `Figure_3A_pathfixed_current_pc.R`
- `Figure_3B_pathfixed_current_pc.R`
- `Figure_4_pathfixed_current_pc.R`

## Recommended running order

1. Run the four EEG network scripts.
2. Run `05_topz_glmm_main_clean_final_new_EEG_paths.R`.
3. Run `hub_mediation_20230401_clean_TOPZ_only.R`.
4. Run the figure scripts.

## Notes

- This repository keeps the manuscript mainline only.
- TopZ definition, TopZ rate calculation, and GLMM settings were preserved for reproducibility.
- Figure 1 is not included in the final clean mainline repository because Panel A depends on the classic hub branch.
