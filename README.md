# NIC EEG Clean Reproduction Repository (Shell)

This is an empty repository skeleton for the NIC EEG clean project.

It contains only the directory structure needed for the final GitHub repository.
No scripts, data, results, or figures are included in this shell.

## Suggested workflow
1. Put raw input files into `01_raw_data/`.
2. Put auxiliary files into `02_auxiliary/`.
3. Put original scripts into `03_original_scripts/`.
4. Put path-fixed scripts into `04_pathfixed_scripts/`.
5. Put final manuscript mainline scripts into `scripts_clean/`.
6. Put generated outputs into `05_results/`.

## Canonical EEG network output folders
Under `05_results/network/`, use these four folder names:
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
