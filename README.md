# Ship Hull Waterline Analysis

This repository contains a compact Python workflow for ship hull waterline analysis, polynomial fitting, and GZ-related calculations based on vessel offset data.

## Overview

The project is organized around three main analysis scripts:

- `tchebychev final.py` — fits polynomial curves to the ship’s waterlines and generates intersection-point data.
- `krylov method final.py` — derives waterline intersections and writes supporting report/CSV outputs.
- `GZ curves final.py` — computes waterplane and stability-related parameters from the generated intersection data.

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
