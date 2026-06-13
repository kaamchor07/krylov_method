# Ship Hull Waterline Analysis

This repository contains a small set of Python scripts for hull-waterline analysis, polynomial fitting, and GZ / stability-related calculations based on ship offset data.

## What is included

- `GZ curves final.py` — computes and analyzes waterplane / GZ-related parameters from the generated intersection data.
- `krylov method final.py` — derives waterline intersections and writes the supporting report / CSV outputs.
- `tchebychev final.py` — fits waterline polynomials and generates the intersection-point dataset used by the other scripts.
- `Ship Offsets.csv` and `stem_stern_offsets.csv` — source hull offset data.
- `Krylov Method Report.pdf` — PDF report generated from the workflow.

## Project purpose

The scripts are used to:

1. fit waterline polynomial curves from offset data,
2. calculate station / waterline intersection points,
3. generate intermediate CSV and report outputs for further hydrostatic or stability analysis.

## Requirements

Install the Python dependencies with:

```bash
pip install -r requirements.txt
```

Expected packages:

- numpy
- pandas
- matplotlib
- scipy

## How to run

From the project root:

```bash
python "tchebychev final.py"
python "krylov method final.py"
python "GZ curves final.py"
```

The scripts write generated outputs such as:

- `tchebychev_intersection_points.csv`
- `tchebychev_intersection_report.txt`
- `waterline_equations.json`
- `waterline_equations.txt`
- `waterline_fits.png`

## Notes

- The source files include the original analysis scripts and their generated artifacts.
- The code is provided as-is for research / engineering analysis workflows.
