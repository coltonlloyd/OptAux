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
run OptAux for 4 `competing_metabolite_uptake_thresholds` (0, 0.01, 0.1, and 2) and
output the results in `supplement_1_optaux_solutions.xls` as well as
intermediate results as `optaux_intermediate_trace_[threhold_value].xls`

To process the results into an MSE summary spreadsheet and the set of
EBC designs with the smallest number of knockouts, run
`output_optaux_summaries.xls`

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
(Not provided due to filesize limits on github. They can be obtained by
contacting cjlloyd@ucsd.edu for access to the alignment files or accessing
the raw reads hosted on the Sequence Read Archive under accession no. SRP161177).

Running `python output_duplications.py` will:

1. Check if `[optaux]/scripts_figures_and_tables/duplications/[pair_name]_coverage_dict.json`
exists. If not it will produce this dictionary using the bam files output from breseq.

2. Find the genes with >80% of their base pairs above the 1.25x fit mean cutoff.
These are compiled and output in 
`[optaux]/scripts_figures_and_tables/duplications/duplicated_genes`

3. Output the plots used to create Figure 7 and the supplementary figures in
`[optaux]/scripts_figures_and_tables/duplications/`

### Community Modeling Sims and Plotting
The iJL1678b ME-model (constructed using COBRAme/ECOLIme v0.0.9) is provided
as json files with two different k<sub>eff</sub> parameter sets:

  - `iJL1678b.json`: The model with default parameters

  - `iJL1678b_null_keffs.json`: The model with k<sub>eff</sub>s values set to
those obtained from [Dividi et. al.](http://www.pnas.org/content/113/12/3401).
Metabolic k<sub>eff</sub>s not included in this dataset were imputed
with the median value, 6.2 s<sup>-1</sup>. Non-metabolic reactions outside the scope
of this dataset were set to 65 s<sup>-1</sup> as in the default model.

The key results from the study can be reproduced with the following
1. Build iJL1678b-community models with: `python [optuax]/optaux/me_community/make_me_communityl.py`.
   - This will make community models based on the two ME-models described
   as well as community model with all k<sub>eff</sub>s set to 65 s<sup>-1</sup>.

2. All of the simulations needed to reproduce Figure 8 can then be ran with:
`python run_community_me_sims.py`

   - This code uses python's multiprocessesing functionality with 2 processes 
    by default. To speed up these simulations, more processes can be used. A 
    single simulation typically requires 4-8 gb of RAM to solve with 
    qMINOS.
   - This will output the results into `[optaux]/scripts_figures_and_tables/community_sims_output_[null/65/default]_keffs`.
Alternatively, tar.gz files are containing the output of these simulations are 
already included in this package. **Note:** If using these files instead of
running new simulations you must run `python unpack_community_me_sims.py` to unpack
the tar.gz files in order to plot.

To recreate Figure 8, part of Figure 7, and the supplementary community ME-figures:

3. To plot community growth rates for varying strain abundances run:
```python [optaux]/scripts_figures_and_tables/output_computed_community_growth_rates.py```

4. To plot metabolite cross-feeding for varying strain abundances run:
```python [optaux]/scripts_figures_and_tables/output_computed_metabolite_crossfeeding.py```

5. To plot community growth rates for substrate limited ME-model and M-model run:
``python [optaux]/scripts_figures_and_tables/output_glucose_limited_me_m_comparison.py``

6. To output the supplement community growth comparison between steadycom and jointfba M-model:
``python [optaux]/scripts_figures_and_tables/output_steadycom_jointfba_comparison.py``

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