# NIC EEG Clean Reproduction Repository

This repository contains the final clean scripts for the NIC EEG manuscript mainline analysis.

It is designed for clean reproduction of the manuscript-linked workflow using the final clean scripts only.

## Repository structure

- `scripts_clean/` final clean scripts for the manuscript mainline analysis
- `README.md` usage instructions

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

## Scripts included

### EEG network
- `EEG_network_NeuroEPO_post_clean_main_english_prompt.R`
- `EEG_network_NeuroEPO_pre_clean_main_english_prompt.R`
- `EEG_network_Placebo_post_clean_main_english_prompt.R`
- `EEG_network_Placebo_pre_clean_main_english_prompt.R`

### TopZ + GLMM
- `topz_glmm_main_clean_final_new_EEG_paths.R`

### Mediation
- `hub_mediation_20230401_clean_TOPZ_only.R`

### Figures
- `Figure_1.R`
- `Figure_2.R`
- `Figure_3A .R`
- `Figure_3B .R`
- `Figure_4.R`

## Recommended running order

1. Run the four EEG network scripts.
2. Run `topz_glmm_main_clean_final_new_EEG_paths.R`.
3. Run `hub_mediation_20230401_clean_TOPZ_only.R`.
4. Run the figure scripts.

## Notes

- This repository keeps the manuscript mainline only.
- The clean pipeline is intended to preserve the manuscript-linked analysis logic.
- EEG network output folders should follow the canonical folder names listed above.
