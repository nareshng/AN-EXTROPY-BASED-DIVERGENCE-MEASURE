NONPARAMETRIC INFERENCE FOR AN EXTROPY-BASED DIVERGENCE MEASURE
CODE AND DATA SUPPLEMENT
==============================================================

Authors
-------
Naresh Garg, Isha Dewan and Sudheesh Kumar Kattumannil

Repository and archived version
--------------------------------
Repository: https://github.com/nareshng/AN-EXTROPY-BASED-DIVERGENCE-MEASURE
Final commit: AUTHOR ACTION REQUIRED
Archived release/DOI: AUTHOR ACTION REQUIRED

Scope
-----
This supplement contains the code and data needed to reproduce Tables 1-14
and Figures 1-4 in the manuscript. It also contains intermediate simulation
results, run logs, random-number seeds, checksums and software-environment
records from the final manuscript run.

Reproduction status
-------------------
AUTHOR ACTION REQUIRED: Replace this paragraph after the full run and the
independent clean-machine verification have passed. State the date, final
commit hash, tester, operating system, R version and whether every printed
number and figure was reproduced.

No manual source-code editing is required. All settings are supplied by the
master runner or command-line arguments. Quick-mode output is only a software
check and must not be used in the manuscript.

Required software
-----------------
Reference environment:

- Operating system: AUTHOR ACTION REQUIRED
- R 4.6.1
- emplik 1.3-3
- survival 3.5-8
- TH.data 1.1-2
- jpeg 0.1-10
- ggplot2 3.4.4
- survminer 0.4.9

The complete dependency snapshot is stored in renv.lock. The exact loaded
packages, locale and platform for the final run are recorded in
results/<final-run>/environment/sessionInfo.txt and environment.txt.

Create the environment with the instructions in ENVIRONMENT_SETUP.txt. The
analysis scripts load packages but do not install or update them.

Folder structure
----------------
.
|-- README.txt
|-- renv.lock
|-- ENVIRONMENT_SETUP.txt
|-- run_all.R
|-- capture_environment.R
|-- make_manuscript_tables.R
|-- verify_submission.R
|-- DATA_PROVENANCE.txt
|-- CLEAN_MACHINE_TEST_RECORD.txt
|-- MONTE_CARLO_STABILITY_RECORD.txt
|-- INDEPENDENT_SEED_RUNS.txt
|-- MANUSCRIPT_CODE_CROSSWALK.csv
|-- REPRODUCIBILITY_CHECKLIST_STATUS.csv
|-- CHECKLIST_SUMMARY.txt
|-- MANUSCRIPT_REVISIONS_REQUIRED.txt
|-- Uncensored Simulations/
|   |-- run_tables_1_6.R
|   |-- validate_tables_1_6.R
|   |-- test_tables_1_6_determinism.R
|   `-- simulation source files for Tables 1-6
|-- Censored Simulations/
|   |-- run_tables_7_10.R
|   |-- Point_estimation_right_censoring.R
|   |-- Confidence_Intervals_right_censoring.R
|   |-- Aggregate_right_censoring_results.R
|   `-- shared functions, diagnostics and atomic I/O
|-- Real Data Analysis/
|   |-- Censored_Real_Data_Analysis.R
|   |-- Image_based_Real_Data_Analysis.R
|   |-- shared analysis and figure functions
|   `-- Images/
|       |-- image_manifest.csv
|       `-- nine analysed MRI images
`-- results/
    `-- final full-run code, intermediate results and paper outputs

One-command reproduction
------------------------
Run from any directory. The master script resolves all paths relative to its
own location.

Quick software check:

  Rscript --vanilla run_all.R --mode=quick --cores=2 --seed-offset=0

Full manuscript run:

  Rscript --vanilla run_all.R --mode=full --cores=4 --seed-offset=0

On Windows, use --cores=1. The master runner creates a new timestamped output
directory and never overwrites an existing run.

Approximate full-run time
-------------------------
- Tables 1-6: approximately 1-2 CPU-hours in total.
- Tables 7-8: approximately 5 minutes.
- Tables 9-10: approximately 2-3 CPU-hours.
- Table 11 and Figures 1-3: approximately 1 minute.
- Tables 12-14 and Figure 4: a few seconds.

Actual hardware and elapsed times for the final run:
AUTHOR ACTION REQUIRED

Table and figure map
--------------------
- Tables 1-6: Uncensored Simulations/run_tables_1_6.R
- Tables 7-10: Censored Simulations/run_tables_7_10.R
- Table 11 and Figures 1-3:
  Real Data Analysis/Censored_Real_Data_Analysis.R
- Tables 12-14 and Figure 4:
  Real Data Analysis/Image_based_Real_Data_Analysis.R

The master workflow calls make_manuscript_tables.R after all analyses. Exact
manuscript-formatted tables are collected in paper_outputs/ as
Table_1_manuscript.csv through Table_14_manuscript.csv. Figure panels are
copied to the same folder using manuscript figure numbers.

Random-number generation
------------------------
All simulations use explicit seeds and an explicit R RNG configuration.
Tables 1-6 use the table-specific seeds recorded in their run manifest.
Tables 7-8 use base seed 5401. Tables 9-10 use base seed 5402. Table 11 uses
base seed 2026. Cell-level seeds and source checksums are saved automatically.

The final supplement must also contain the independent-seed Monte Carlo
stability report demonstrating that conclusions and method rankings are not
materially changed by Monte Carlo error.

Intermediate results and spot checks
------------------------------------
The long censored simulations save every completed cell as an atomic
checkpoint. The final supplement includes these checkpoints, replication-level
results, task manifests and logs so an editor can rerun a selected setting
without repeating the entire study. Instructions for selecting task IDs are in
Censored Simulations/README.md.

Data sources
------------
1. Veteran lung cancer data: survival::veteran.
2. NCCTG lung cancer data: survival::lung.
3. GBSG2 breast cancer data: TH.data::GBSG2.
4. Nine MRI images: selected from the Kaggle source documented in
   Real Data Analysis/Images/image_manifest.csv.

The manuscript run requires TH.data::GBSG2. It must not silently substitute
survival::gbsg because the row ordering changes the seeded bootstrap result.
Complete provenance, documentation links, licensing information and image
selection details are recorded in DATA_PROVENANCE.txt and image_manifest.csv.

Verification
------------
After a full run, execute:

  Rscript --vanilla verify_submission.R --run-root=<FULL_RUN_DIRECTORY>

The command checks required sources, completed provenance fields, environment
records, manuscript-formatted outputs and checksums. A co-author must then
repeat the full workflow on a clean machine and complete
CLEAN_MACHINE_TEST_RECORD.txt.

Known inferential scope
-----------------------
The kernel estimator is used for point estimation only. The right-censored
analysis targets divergence over a fixed follow-up interval. Percentile
bootstrap intervals under censoring are diagnostic comparators and are not
claimed to be uniformly valid near a degenerate null or under sparse tail
follow-up.
