# Ship Hull Waterline Analysis

This repository contains a compact Python workflow for ship hull waterline analysis, polynomial fitting, and GZ-related calculations based on vessel offset data.

## Overview

The project is organized around three main analysis scripts:

The workflow starts with the hull offset tables in the CSV files, fits polynomial waterline curves from those station/half-breadth values, and then uses those fitted surfaces to compute geometry and stability-related quantities. In short: raw offset data → polynomial waterline fitting → intersection points → hydrostatic / GZ analysis.

- `tchebychev final.py` — fits polynomial curves to the ship’s waterlines and generates intersection-point data.
- `krylov method final.py` — derives waterline intersections and writes supporting report/CSV outputs.
- `GZ curves final.py` — computes waterplane and stability-related parameters from the generated intersection data.

A brief summary of the method:

- `tchebychev final.py` reads the ship offset tables, builds polynomial fits for each waterline, and writes the intersection-point dataset used downstream.
- `krylov method final.py` uses those fitted equations to locate the ship’s waterline intersections, identify station positions, and generate a structured report.
- `GZ curves final.py` uses the resulting intersection geometry to evaluate the waterplane characteristics and related GZ/stability quantities.

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
