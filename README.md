# sctrna-reproducibility
Repository for reanalyzing and reproducing single-cell tRNA-seq analyses


file_structure/
├── git_repo/
│   ├── global_functions
│   │         ├── tRNA_import.R
│   │         └── tRNA_valid.R
│   └── project
│                ├── code
│                │        └── setup.R
│                └── output
└── trna_large_data
                   ├── giftwrap
                   │        ├── Sample_A
                   │        │        ├── barcodes.tsv.gz
                   │        │        ├── manifest.tsv
                   │        │        └── probe_reads.tsv.gz
                   │        └── Sample_B
                   │                   ├── barcodes.tsv.gz
                   │                   ├── manifest.tsv
                   │                   └── probe_reads.tsv.gz
                   ├── mRNA
                   │        ├── Sample_A.wta.h5
                   │        └── Sample_B.wta.h5
                   └── output

                   
