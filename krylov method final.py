import numpy as np
import pandas as pd
import json
from scipy.optimize import fsolve
import os
script_dir = os.path.dirname(__file__)

with open(os.path.join(script_dir, 'waterline_equations.json'), 'r') as f:
    waterline_equations = json.load(f)

def eval_poly(coeffs, x):
    result = 0
    n = len(coeffs) - 1
    for i, coef in enumerate(coeffs):
        power = n - i
        result += coef * (x ** power)
    return result

wl_1507_key = '15.07'

coeffs_1507 = waterline_equations[wl_1507_key]['coefficients']

def poly_func(x):
    return eval_poly(coeffs_1507, x)

x_search = np.linspace(-50, 300, 10000)
y_search = [poly_func(x) for x in x_search]

sign_changes = []
for i in range(len(y_search) - 1):
    if y_search[i] * y_search[i+1] < 0:
        root = fsolve(poly_func, x_search[i])[0]
        sign_changes.append(root)

zeros = []
for root in sign_changes:
    is_duplicate = False
    for existing_zero in zeros:
        if abs(root - existing_zero) < 0.1:
            is_duplicate = True
            break
    if not is_duplicate:
        zeros.append(root)

zeros.sort()

if len(zeros) < 2:
    print("Error: Could not find 2 zeros for waterline 15.07")
    exit(1)

zero_1 = zeros[0]
zero_2 = zeros[1]

LWL = abs(zero_2 - zero_1)
LWL_half = (zero_1 + zero_2) / 2

print(f"\nFirst zero (forward): {zero_1:.4f}")
print(f"Second zero (aft): {zero_2:.4f}")
print(f"LWL: {LWL:.4f}")
print(f"LWL/2: {LWL_half:.4f}")

lwl_half_fractions = [0.0669, 0.2887, 0.3667, 0.6330, 0.7113, 0.9331]

offsets_from_lwl_half = []
offsets_from_lwl_half.extend([-frac for frac in reversed(lwl_half_fractions)])
offsets_from_lwl_half.extend(lwl_half_fractions)

vertical_lines = [LWL_half + offset * (LWL/2) for offset in offsets_from_lwl_half]

for i, (offset, x_val) in enumerate(zip(offsets_from_lwl_half, vertical_lines)):
    side = "Port" if offset < 0 else "Starboard"
    distance_from_half = abs(offset * (LWL/2))
    print(f"Line {i+1:2d}: {offset:+.4f} × LWL/2 ({distance_from_half:.4f} from LWL/2) → x = {x_val:.4f} ({side})")

intersections = {}
waterline_bounds = {}

for wl_name, wl_data in waterline_equations.items():
    coeffs = wl_data['coefficients']
    intersections[wl_name] = {}
    
    def wl_func(x):
        return eval_poly(coeffs, x)
    
    x_search_range = np.linspace(-50, 300, 5000)
    y_search_range = [wl_func(x) for x in x_search_range]
    
    wl_zeros = []
    for i in range(len(y_search_range) - 1):
        if y_search_range[i] * y_search_range[i+1] < 0:
            root = fsolve(wl_func, x_search_range[i])[0]
            is_new = True
            for existing_root in wl_zeros:
                if abs(root - existing_root) < 0.5:
                    is_new = False
                    break
            if is_new:
                wl_zeros.append(root)
    
    wl_zeros.sort()
    
    if len(wl_zeros) >= 2:
        x_min_valid = wl_zeros[0]
        x_max_valid = wl_zeros[-1]
        waterline_bounds[wl_name] = (x_min_valid, x_max_valid)
        print(f"WL {wl_name}: Valid range x ∈ [{x_min_valid:.4f}, {x_max_valid:.4f}]")
    else:
        waterline_bounds[wl_name] = (-1000, 1000)
        print(f"WL {wl_name}: No clear x-axis intersections, using full range")
    
    x_min_valid, x_max_valid = waterline_bounds[wl_name]
    
    for x_val in vertical_lines:
        if x_val < x_min_valid or x_val > x_max_valid:
            intersections[wl_name][x_val] = 0
        else:
            y_val = eval_poly(coeffs, x_val)
            
            if y_val < 0:
                y_val = 0
            elif 18.10 < y_val <= 20:
                y_val = 18.15
            elif y_val > 20:
                y_val = 0
            
            intersections[wl_name][x_val] = y_val

intersection_data = []
for wl_name in sorted(waterline_equations.keys(), key=float):
    row = {'Waterline': wl_name}
    for i, x_val in enumerate(vertical_lines):
        row[f'Station_{i+1}'] = intersections[wl_name][x_val]
    intersection_data.append(row)

df_intersections = pd.DataFrame(intersection_data)
df_intersections.to_csv(os.path.join(script_dir, 'tchebychev_intersection_points.csv'), index=False)

with open('tchebychev_intersection_report.txt', 'w') as f:
    f.write("=" * 100 + "\n")
    f.write("TCHEBYCHEV INTERSECTION POINTS\n")
    f.write("=" * 100 + "\n\n")
    f.write(f"LWL (from WL 15.07): {LWL:.4f}\n")
    f.write(f"LWL/2 (Midpoint): {LWL_half:.4f}\n")
    f.write(f"Forward zero: {zero_1:.4f}\n")
    f.write(f"Aft zero: {zero_2:.4f}\n\n")
    
    f.write("Vertical Lines:\n")
    f.write("-" * 80 + "\n")
    for i, (offset, x_val) in enumerate(zip(offsets_from_lwl_half, vertical_lines)):
        f.write(f"  Station {i+1:2d}: {offset:+.4f} LWL → x = {x_val:.4f}\n")
    f.write("\n")
    
    for i, x_val in enumerate(vertical_lines):
        offset = offsets_from_lwl_half[i]
        f.write(f"\nStation {i+1} (x = {x_val:.4f}, offset = {offset:+.4f} LWL):\n")
        f.write("-" * 80 + "\n")
        for wl_name in sorted(waterline_equations.keys(), key=float):
            y_val = intersections[wl_name][x_val]
            f.write(f"  WL {wl_name:8s} : y = {y_val:12.6f}\n")
        f.write("\n")

print(f"\nNumber of waterlines: {len(waterline_equations)}")
print(f"Number of stations: {len(vertical_lines)}")
print(f"Total intersection points: {len(waterline_equations) * len(vertical_lines)}")

print("Results saved to:")
print("tchebychev_intersection_points.csv")
print("tchebychev_intersection_report.txt")
 