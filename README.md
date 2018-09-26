# OptAux

Files and scripts needed to reproduce the figures in "Model-driven
design and evolution of non-trivial synthetic syntrophic pairs"

## Docker
Project can be installed using docker by running from the docker folder
```
docker build -t optaux:everything .
```

this will install everything needed to reproduce most of the results outlined
below including the community iJL1678b ME-model.

## Reproducing Figures/Tables
### OptAux simulations
Run `python [optaux]/scripts_figures_and_tables/make_optaux_supplement.py`

This runs OptAux algorithm in aerobic glucose minimal media conditions for
all carbon containing exchange reactions in iJO1366. It will by default
run OptAux for 4 `trace_metabolite_thresholds` (0, 0.01, 0.1, and 2) and
output the results in `supplement_1_optaux_solutions.xls`

### Relative abundance approximations
Running `python output_relative_abundance_results.py` will:

1. Check if `abundance_by_characteristic_[pair_name].csv`
and `abundance_by_coverage_[pair_name].csv` are located in
`[optaux]/scripts_figures_and_tables/relative_abundance/tables`. If
not, these CSVs will be created using the information in the read Bowtie2
alignment BAM files for all sequencing samples (not provided for now)
and breseq generated mutation calls found in
`[optaux]/optaux/resources/resequencing_data/AUX [pair_name] Mutation Table.csv`
      -  Some values in the Mutation Tables for the characteristic strain
      mutations are missing because breseq
      rounds high/low mutation frequencies to 100% or 0%, respectively.
      These values are filled in with their frequencies in
      `relative_abundance.py`

2. Output bar charts of the characteristic mutation/alignment coverage
base approximations of relative strain abundances in
`[optaux]/scripts_figures_and_tables/relative_abundance/figures` as well
as a comparison of the predictions based on each method.


### Duplications
Read coverages are plotted using the alignment files produced from breseq 
(not provided for now). They coverage of each base pair, however, is included 
as a compressed dictionary and must first be unpacked by running `python unpack_read_coverage_dict`

Running `python output_duplications.py` will:

1. Check if `[optaux]/scripts_figures_and_tables/duplications/[pair_name]_coverage_dict.json`
exists. If not it will get this dictionary using the bam files output from breseq.

2. Find the genes with >80% of their base pairs above the 1.25x fit mean cutoff.
These are compiled and output in 
`[optaux]/scripts_figures_and_tables/duplications/duplicated_genes`

3. Output the plots used to create Figure 7 and the supplementary figures in
`[optaux]/scripts_figures_and_tables/duplications/`

### Community ME-model Sims and Plotting
To run the community ME-model simulations the, community ME-model must first be
constructed (if the docker image is used, the community model is already built)

1. Build iJL1768b-ME by running `python [ecolime]/ecolime/build_me_model.py`
2. Copy the model to optaux with:
 ```cp [ecolime]/ecolime/me_models/iJL1678b.pickle [optaux]/optaux/resources/```
3. Build iJL1678b-community with: `python [optuax]/optaux/me_community/make_me_communityl.py`

4. All of the simulations needed to reproduce Figure 8 can then be ran with:
`python run_community_me_sims.py`

   - This code uses python's multiprocessesing functionality with 2 processes 
    by default. To speed up these simulations, more processes can be used. A 
    single simulation typically requires 4-8 gb of RAM to solve with 
    qMINOS.
   - This will output the results into `[optaux]/scripts_figures_and_tables/community_me_sims`.
Alternatively, tar.gz files are containing the output of these simulations are 
already included in this package. **Note:** If using these files instead of
running new simulations you must run `python unpack_community_me_sims.py` to unpack
the tar.gz files in order to plot.

To recreate Figure 8 and the supplementary community ME-figures:

5. To plot community growth rates for varying strain abundances run: 
```python [optaux]/scripts_figures_and_tables/output_computed_community_growth_rates.py```

6. To plot metabolite cross-feeding for varying strain abundances run: 
```python [optaux]/scripts_figures_and_tables/output_computed_metabolite_crossfeeding.py```


## Software
The following software and versions were used for publication:

- Python 3.6
- An MILP solver. The scripts here use Gurobi by default.
- COBRApy v0.5.11
- [COBRAme](https:/github.com/sbrg/cobrame) v0.0.9
- [ECOLIme](https:/github.com/sbrg/ecolime) v0.0.9
- [solvemepy](https:/github.com/sbrg/solvemepy) v1.0.1
    - Including the qMINOS solver
- pysam v0.14.1
- openpyxl v2.3.2
- pandas v0.22.0
- matplotlib v2.0.2
- numpy v1.14.2
- scipy v0.19.0
- Biopython v1.66
- Seaborn v0.7.1