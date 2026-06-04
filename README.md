# scrna-coarsening-sims

Simulation and real-data analysis code for the paper:

> Alexander Towell. "Coarsening-at-random conditions for scRNA-seq zero
> inflation: a reliability-theoretic perspective." Preprint, Zenodo,
> DOI: 10.5281/zenodo.20414734.

The paper recasts the single-cell RNA-seq dropout / zero-inflation controversy
as a masked-data coarsening-at-random identifiability problem. This repository
contains every script behind the paper's numerical claims: the core simulation
library, the per-study drivers, and the real-data analyses on Tabula Muris and
the Klein et al. (2015) inDrops data.

## Layout

- `sim.R`: core simulation library (data-generating processes, estimators,
  diagnostics). Sourced by the driver scripts.
- `run.R`, `run_pb.R`, `run_v3*.R` ... `run_v11*.R`: per-study drivers. Each
  writes a `results_*.rds` and a `log_*.txt`. The `v*` numbering matches the
  studies described in the paper.
- `run_v8_tabula_muris.R`, `run_v10_tabula_muris_cross_tissue.R`: real-data
  application to Tabula Muris FACS Smart-seq2 (spleen, marrow, liver).
- `run_v11_klein_decent.R`, `run_v11b_klein_decent_v2.R`,
  `run_v11c_klein_compare.R`: DECENT head-to-head on Klein et al. (2015)
  inDrops UMI data (GEO GSE65525).
- `make_figures.R`: regenerates the paper figures from the `results_*.rds`.
- `results_*.md`: human-readable write-ups of each study.
- `results_*.rds`: derived results (small; kept for reproducibility).

## Data (not included; fetch from canonical sources)

The raw datasets are public and are NOT redistributed here (see `.gitignore`).
To reproduce the real-data analyses, fetch them into the indicated directories:

- Tabula Muris FACS (Smart-seq2) counts and annotations, into `tabula_muris/`:
  figshare DOI 10.6084/m9.figshare.5715040 (The Tabula Muris Consortium, 2018).
  The drivers expect `tabula_muris/FACS/<Tissue>-counts.csv` and
  `tabula_muris/annotations_facs.csv`.
- Klein et al. (2015) inDrops data, into `klein2015/`: NCBI GEO accession
  GSE65525. The DECENT comparison uses the K562 cells matrix
  (`GSM1599500_K562_cells.csv`) and the ERCC spike-in concentrations.

`dataset_investigation.md` records exactly which files each analysis reads.

## Reproducing

```r
# core simulation studies (no external data needed)
Rscript run.R
Rscript run_v3.R
# ... etc; each driver is self-contained and sources sim.R

# real-data studies (after fetching the data above)
Rscript run_v8_tabula_muris.R
Rscript run_v10_tabula_muris_cross_tissue.R
Rscript run_v11_klein_decent.R

# figures
Rscript make_figures.R
```

Base R is sufficient for the simulation studies. The DECENT comparison requires
the `DECENT` package; `run_v11_klein_decent.R` documents the install and the
units convention (DECENT expects molecules-per-cell, not attomoles).

## Author

Alexander Towell, lex@metafunctor.com, ORCID 0000-0001-6443-9897,
Southern Illinois University Edwardsville.

## License

Code: MIT. The referenced datasets retain their original licenses and should be
obtained from the canonical sources above.
