## Single-cell translatome and RNA sequencing (sctRNA-seq) analysis repository
Repository for reanalyzing and reproducing single-cell tRNA-seq analyses

### Background and gap filling framework
<img width="1259" height="314" alt="Screenshot 2026-11-05 at 9 05 48 PM" src="https://github.com/user-attachments/assets/a7f0f6b1-a4d9-44ad-a70f-5b6d2d375336" />

<br>
<br>

### Data analysis
Analyses from all main figures are coordinated in the following folders:
```
Figure 1 | 4cellline_noBetaine,4cellline_Betaine
Figure 2 | perturb_pilot
Figure 3 | perturb_protein-coding
Figure 4 | cd34_differentiation
Figure 5 | Tcell_stim
Figure 6 | perturb_treatment
```

<br>

### Setup
Code and data files are organized in the following this path structure:
```
file_structure/
├── git_repo/
│   ├── global_functions
│   │   ├── tRNA_import.R
│   │   └── tRNA_valid.R
│   └── project
│       ├── code
│       │   └── setup.R
│       └── output
└── trna_large_data
    ├── giftwrap
    │   ├── Sample_A
    │   │   ├── barcodes.tsv.gz
    │   │   ├── manifest.tsv
    │   │   └── probe_reads.tsv.gz
    │   └── Sample_B
    │       ├── barcodes.tsv.gz
    │       ├── manifest.tsv
    │       └── probe_reads.tsv.gz
    ├── mRNA
    │   ├── Sample_A.wta.h5
    │   └── Sample_B.wta.h5
    └── output
```
