# Mithra: A Cross-Platform Toolkit for Intracranial Electrode Visualization

**Mithra** is a unified Python and MATLAB toolkit for visualizing intracranial electrophysiological recordings, including ECoG and sEEG. It enables subject-specific visualization, standardized-space mapping, and systematic within- and across-subject analysis.

## Overview

Intracranial recordings provide high spatiotemporal resolution but require careful handling of anatomy, coordinate systems, and cross-subject comparisons. Existing tools often address only isolated components of this workflow. Mithra provides an integrated solution:

- Visualize electrodes on each participant’s pial surface
- Preserve precise subject-specific anatomy
- Map electrodes to common average spaces (MNI, FreeSurfer average)
- Load FreeSurfer annotations for anatomical context
- Project electrodes to the cortical surface
- Generate Gaussian heatmaps of spatial activation

For scientific details, refer to the preprint:
**bioRxiv:** [Mithra](https://www.biorxiv.org/content/10.1101/2025.11.20.689539v1)

### Data Formats
- FreeSurfer surfaces (`.pial`, `.inflated`, `fsaverage` equivalents, `.mat`)
- Electrode coordinate files (CSV, MATLAB structs, JSON, mat )
- FreeSurfer annotation files (`.annot`)

### Citation

```bib
@article {Mithra_2025,
	author = {Khalilian-Gourtani, Amirhossein and Esmaeili, Yasamin and Michalak, Andrew J and Flinker, Adeen},
	title = {Mithra: An Open-Source and Cross-Platform Visualization Toolbox for Human Intracranial Recordings},
	elocation-id = {2025.11.20.689539},
	year = {2025},
	doi = {10.1101/2025.11.20.689539},
	URL = {https://www.biorxiv.org/content/early/2025/11/21/2025.11.20.689539},
	eprint = {https://www.biorxiv.org/content/early/2025/11/21/2025.11.20.689539.full.pdf},
	journal = {bioRxiv}
}
```

