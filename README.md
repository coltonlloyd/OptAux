# OptAux

Files and scripts needed to reproduce the figures in "Model-driven
design and evolution of non-trivial synthetic syntrophic pairs"

## Docker
Project can be installed using docker by running from the docker folder
```
docker build -t optaux:everything .
```

this will install everything needed to reproduce most of the results outlined
below including the community iJL1678b-ME model.

## Reproducing Figures/Tables
### OptAux simulations
Run `python make_optaux_supplement.py`

This runs OptAux algorithm in aerobic glucose minimal media conditions for
all carbon containing exchange reactions in iJO1366. It will by default
run OptAux for 4 `trace_metabolite_thresholds` (0, .01, .1, and 2) and
output the results in `supplement_1_optaux_solutions.xls`

### Relative abundance approximations
By running `python output_relative_abundance_results.py` this will:

1. Check if `abundance_by_characteristic_[pair_name].csv`
and `abundance_by_coverage_[pair_name].csv` are located in
`[optaux]/scripts_figures_and_tables/relative_abundance/tables`. If
not, these CSVs will be created using the information in the read Bowtie2
alignment BAM files for all sequencing samples (not provided for now)
and breseq generated mutation calls found in
`[optaux]/optaux/resources/resequencing_data/AUX [pair_name] Mutation Table.csv`

2. Output bar charts of the characteristic mutation/alignment coverage
base approximations of relative strain abundances in
`[optaux]/scripts_figures_and_tables/relative_abundance/figures` as well
as a comparison of the predictions based on each method.

## Software
The following software and versions were used for publication:

- Python 2.7/3.6
- An MILP solver. We recommend Gurobi.
- [cobrame](https:/github.com/sbrg/cobrame) v0.0.7
- [ecolime](https:/github.com/sbrg/ecolime) v0.0.7
- [solvempy](https:/github.com/sbrg/solvemepy) v1.0.1
    - Including the qMINOS solver
- pysam v0.14.1
- openpyxl v2.3.2
- pandas v0.22.0
- matplotlib v2.0.2
- numpy v1.14.2
- scipy v0.19.0
- Biopython v1.66