# hiv_polyA_struct_phylo_analysis

### Goal: to analyze the structural conservation of the 5' polyA hairpin the HIV genome.

`insert figure here (i.e. graphical abstract)`

#### Approach
There are three main components of the pipeline and each is seperated into a seperate directory in the repository. Each component is described further. 

- ##### polyA_detection
- ##### make_tree-post_detection_analysis
- ##### concensus_structure-post_detection_analysis

---
### polyA_detection
Goal: Parse out the putative polyA hairpin from a set of target HIV genomes using de novo discovery method described below. 

`insert graphic describing discovery pipeline`

Requirements: `viennaRNA`, `clustalo`

Input file (parameters are explained in '<< >>' - do not include brackets in file when running):
```json
{
    "paths": {
        "full_genomes_path": ""<<Path to a JSON file with the target genomes. See example for input format.>>,
        "rna_fold_path": ""<<Path to viennaRNA executable (RNAfold)>>,
        "clustal_path": "" <<Path to clustalo executable>>
    },
    "initial_seg": {
        "ext_range": 75 <<# of basepairs to extend up and down stream from the AAUAAA signal to pull the initital fragment to pull the polyA from>>
    },
    "minimization": {
        "min_len": 35, <<Minimum length of fragment to be conisdered a potential putative polyA hairpin>>
        "max_len": 50, <<Maximum length of fragment to be conisdered a potential putative polyA hairpin>>
        "max_unbp_ratio": 0.45, <<Ratio representing the max percentage of unbase paired residues>>
        "hwd_limit": 0.25 <<A restraint on how far from the center of the hairpin the AAUAAA can be for a potential putative polyA hairpin.>>
    },
    "output": {
        "folder": "../output/", << >>
        "polyA_output": "polyA.fasta" << >>
    },
    "strains_force_include":["KM390026", "K03455", "AB286862", "U51188", "X04415", "MH705144", "AB485658", "MH705163", "KC156211"] <<Strains that will be forced into the final output - regardless of whether or not they are represenative of their subtype.>>
}

```