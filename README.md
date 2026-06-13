# Ship Hull Waterline Analysis

This repository contains a compact Python workflow for ship hull waterline analysis, polynomial fitting, and GZ-related calculations based on vessel offset data.

## Overview

The project is organized around three main analysis scripts:

The workflow starts with the hull offset tables in the CSV files, fits polynomial waterline curves from those station/half-breadth values, and then uses those fitted surfaces to compute geometry and stability-related quantities. In short:

raw offset data → polynomial waterline fitting → intersection points → hydrostatic / GZ analysis.

This is a practical numerical workflow for marine hull geometry: the code converts tabulated ship offsets into smooth mathematical curves, uses those curves to locate important waterline intersections, and then uses the resulting geometry to estimate stability-related quantities.

- `tchebychev final.py` — fits polynomial curves to the ship’s waterlines and generates intersection-point data.
- `krylov method final.py` — derives waterline intersections and writes supporting report/CSV outputs.
- `GZ curves final.py` — computes waterplane and stability-related parameters from the generated intersection data.

### How the method works

1. Data preparation
   - The file `Ship Offsets.csv` contains the main station-by-waterline half-breadth data.
   - `stem_stern_offsets.csv` adds the stem and stern boundary information that helps the fitting process stay realistic at the ends of the hull.

2. Polynomial fitting
   - The script `tchebychev final.py` reads the offset data and fits polynomial curves for each waterline.
   - These fitted curves represent the ship’s hull profile as a smooth mathematical function instead of a set of raw tabulated points.
   - The result is saved as `waterline_equations.json` and `waterline_equations.txt`, and the fitted curves are visualized in `waterline_fits.png`.

3. Intersection and station analysis
   - The script `krylov method final.py` uses the fitted equations to find where the waterlines intersect key station lines.
   - It computes station positions, identifies valid ranges, and generates a structured intersection table in `tchebychev_intersection_points.csv`.
   - A readable report is written to `tchebychev_intersection_report.txt` to explain the steps and numerical results.

4. GZ / stability evaluation
   - The script `GZ curves final.py` takes the generated intersection geometry and evaluates waterplane / stability-related quantities.
   - This stage is where the fitted hull geometry is converted into practical hydrostatic and GZ analysis outputs.

In other words, the code is not just plotting curves; it is building a complete numerical chain from raw hull offsets to an interpretable stability analysis workflow.

It also includes the source datasets used by the analysis:

- `Ship Offsets.csv`
- `stem_stern_offsets.csv`
- `Krylov Method Report.pdf`

## Project goals

This codebase is intended for:

1. fitting waterline polynomial curves from hull offset data,
2. calculating station/waterline intersection points,
3. generating intermediate CSV and report outputs for hydrostatic or stability analysis.

## Requirements

The project uses Python 3 and the dependencies listed in `requirements.txt`.

Install them with:

```bash
pip install -r requirements.txt
```

Core dependencies:

- numpy
- pandas
- matplotlib
- scipy

## Quick start

Run the scripts from the repository root:

```bash
python "tchebychev final.py"
python "krylov method final.py"
python "GZ curves final.py"
```

Typical outputs produced by the workflow include:

- `tchebychev_intersection_points.csv`
- `tchebychev_intersection_report.txt`
- `waterline_equations.json`
- `waterline_equations.txt`
- `waterline_fits.png`

## Repository structure

- `*.py` — analysis and fitting scripts
- `*.csv` — vessel offset and station data
- `*.pdf` — generated report output
- `requirements.txt` — Python package dependencies

## Notes

- This repository includes the original analysis scripts and their generated artifacts.
- The code is provided as-is for research, engineering study, and exploratory hull analysis work.
