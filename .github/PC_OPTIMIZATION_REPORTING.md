# PC Optimization Reporting & Visualization

## Summary of Changes

This update adds comprehensive visualization and detailed reporting of the PC optimization process to the sceQTL-flow pipeline.

### New Files Created

#### 1. **R/quarto_reports/pc_optimization_report.qmd**
Comprehensive Quarto report template with:
- **Coarse Stage Visualization**:
  - Table showing all coarse PC values, associations, and gains
  - Gains plot (% improvement vs PC count)
  - Early-stop threshold and ceiling marked
  
- **Fine Stage Visualization**:
  - Table showing fine PC values, associations, and elbow metrics
  - Fine grid plot with elbow tolerance band
  - Selected PC highlighted
  
- **Summary Table**: Selected PC count for each cell type

- **Parameters Table**: All optimization parameters explained

- **Interpretation Guide**: Explains early stopping, ceiling, elbow tolerance, and final selection

### Files Modified

#### 1. **modules/optimize/select_pcs_coarse.nf**
- Added output block with `coarse_summary.csv` emission
- Added R code to generate detailed CSV with columns:
  - `celltype`, `n_pcs`, `n_assoc`
  - `gain`: absolute improvement from previous PC
  - `gain_pct`: percentage improvement
  - `below_threshold`: TRUE if gain < 1% (early-stop candidate)
  - `is_ceiling`: TRUE at the plateau ceiling
  - `is_best`: TRUE for selected coarse PC

#### 2. **modules/optimize/select_pcs.nf**
- Added output block with `fine_summary.csv` emission
- Added R code to generate detailed CSV with columns:
  - `celltype`, `n_pcs`, `n_assoc`
  - `threshold`: elbow tolerance cutoff (98% of max)
  - `within_elbow`: TRUE if within tolerance region
  - `is_selected`: TRUE for final selected PC

#### 3. **modules/reports/organize_pc_optimization.nf** (NEW)
- New process to collect and organize coarse/fine summary CSVs
- Publishes to `{outdir}/optimization/`
- Runs only when `params.optimize_pcs = true`

#### 4. **workflows/matrixeqtl.nf**
- Added import for `organize_pc_optimization` process
- Added collection of `coarse_summary` outputs from `select_pcs_coarse`
- Added collection of `fine_summary` outputs from `select_pcs`
- Added channels for `collected_coarse_summaries` and `collected_fine_summaries`
- Added else-branch placeholder channels for fixed-PC mode
- Added conditional call to `organize_pc_optimization` when `params.optimize_pcs = true`
- Updated final_report invocation to use unified report template with all PC optimization parameters

### Data Flow

```
[Coarse Stage]
count_individuals → coarse PC values
optimize_pcs_coarse → {celltype}_egenes_vs_{n}_coarse.txt (11 files per cell type)
select_pcs_coarse → {celltype}_coarse_summary.csv + pc_values_fine.txt

[Fine Stage]
optimize_pcs_fine → {celltype}_egenes_vs_{n}_fine.txt (20 files per cell type)
select_pcs → {celltype}_fine_summary.csv + {celltype}_{n}_pcs.txt

[Organization]
organize_pc_optimization → copies CSVs to {outdir}/optimization/

[Reporting]
Quarto report reads all CSVs and generates visualizations
```

### Output Files

#### New Files in `{outdir}/optimization/`:
- `{celltype}_coarse_summary.csv` — Coarse grid results with gains
- `{celltype}_fine_summary.csv` — Fine grid results with elbow metrics

#### Unified HTML Report:
- `{outdir}/eQTL_outputs/eqtl_report.html` — Comprehensive report with eQTL results + PC optimization visualizations
  - **Part 1**: eQTL Analysis Report (settings, input overview, eGene counts)
  - **Part 2**: PC Optimization Results (coarse/fine grids, visualizations, parameter documentation)

### CSV Structure

**Coarse Summary Example:**
```
celltype,n_pcs,n_assoc,gain,gain_pct,below_threshold,is_ceiling,is_best
Microglia,2,1245,NA,NA,FALSE,FALSE,FALSE
Microglia,12,2103,858,68.9,FALSE,FALSE,TRUE
Microglia,22,2456,353,16.8,FALSE,FALSE,FALSE
Microglia,32,2589,133,5.4,FALSE,FALSE,FALSE
Microglia,42,2614,25,1.0,FALSE,FALSE,FALSE
Microglia,52,2625,11,0.4,TRUE,TRUE,FALSE
```

**Fine Summary Example:**
```
celltype,n_pcs,n_assoc,threshold,within_elbow,is_selected
Microglia,12,2103,2562.6,FALSE,FALSE
Microglia,14,2287,2562.6,FALSE,FALSE
Microglia,16,2401,2562.6,FALSE,FALSE
Microglia,18,2478,2562.6,FALSE,FALSE
Microglia,20,2542,2562.6,TRUE,FALSE
Microglia,22,2625,2562.6,TRUE,TRUE
```

### Usage

The unified report is automatically generated when:
1. `--optimize_pcs true` (PC optimization enabled) — includes PC optimization section
2. `--optimize_pcs false` (fixed PCs) — report includes only eQTL results
3. `--report true` (reporting enabled) — triggers HTML report generation

**Note**: The report uses the new `R/rmarkdown_reports/unified_final_report.Rmd` template which combines eQTL analysis and PC optimization sections.

To re-render the report after analysis:
```bash
cd /path/to/results
R -e "rmarkdown::render('path/to/unified_final_report.Rmd', output_file='eqtl_report.html', params=list(outdir='.'))"
```

### Visual Elements

The Quarto report generates:
1. **Coarse gains plot** — Shows how associations increase with PC count
   - Red horizontal line = 1% threshold
   - Red vertical line = early-stop ceiling
   - Orange points = below threshold (consecutive stop candidates)
   - Blue points = above threshold

2. **Fine grid plot** — Shows fine-tuned PC selection
   - Light blue band = elbow tolerance region (±2% of max)
   - Blue dots = tested PC values
   - Green triangle = selected PC (smallest within tolerance)

3. **Summary table** — PC selection outcome for all cell types

4. **Parameters table** — Documents all tuning parameters

### Integration Notes

- ✅ CSV generation does not interfere with existing PC selection logic
- ✅ Workflow compatible with both optimized and fixed-PC modes
- ✅ Gracefully handles absence of optimization (no CSVs in fixed-PC mode)
- ✅ Unified report automatically includes/excludes PC optimization section based on `--optimize_pcs` flag
- ✅ All pipeline parameters passed to final report for comprehensive documentation
- ✅ No additional parameters required (uses existing `--optimize_pcs` and `--report` flags)
- ✅ Report uses `R/rmarkdown_reports/unified_final_report.Rmd` as unified template

### Testing

The pipeline has been tested with:
- ✅ PC optimization enabled (`--optimize_pcs true`)
- ✅ PC optimization disabled (`--optimize_pcs false`)
- ✅ Multiple cell types with varying sample sizes
- ✅ Coarse early-stop triggering at different PC counts

