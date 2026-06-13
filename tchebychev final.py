import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import json
import os

print("SCRIPT STARTED")

script_dir = os.path.dirname(__file__)
df = pd.read_csv(os.path.join(script_dir, 'Ship Offsets.csv'), index_col=0)

stemstern_path = os.path.join(script_dir, 'stem_stern_offsets.csv')
df_stemstern = pd.read_csv(stemstern_path)

stem_waterlines = df_stemstern['WL'].tolist()
stem_offsets = df_stemstern['stem_offset'].tolist()
stern_waterlines = df_stemstern['WL'].tolist()
stern_offsets = df_stemstern['stern_offset'].tolist()

stations = df.index.values
fitted_equations = {}

fig, ax = plt.subplots(figsize=(18, 10))
colors = plt.cm.viridis(np.linspace(0, 1, len(df.columns)))

for idx, col_name in enumerate(df.columns):
    wl = float(col_name)
    print(f"\nProcessing Waterline {wl}...")
    
    offsets_main = df[col_name].values
    stations_main = stations.copy()
    # Trim only consecutive leading and trailing zeros (keep one)
    if wl != 15.07:
        zero_indices = np.where(offsets_main == 0)[0]
        consec_start = 0
        while consec_start + 1 < len(offsets_main) and offsets_main[consec_start] == 0 and offsets_main[consec_start + 1] == 0:
            consec_start += 1
        if consec_start > 0:
            offsets_main = offsets_main[consec_start:]
            stations_main = stations_main[consec_start:]

        zero_indices = np.where(offsets_main == 0)[0]
        consec_end = len(offsets_main) - 1
        while consec_end - 1 >= 0 and offsets_main[consec_end] == 0 and offsets_main[consec_end - 1] == 0:
            consec_end -= 1
        if consec_end < len(offsets_main) - 1:
            offsets_main = offsets_main[:consec_end + 1]
            stations_main = stations_main[:consec_end + 1]
    
    stem_offset = None
    if wl in stem_waterlines:
        stem_idx = stem_waterlines.index(wl)
        stem_offset = stem_offsets[stem_idx]
    
    stern_offset = None
    if wl in stern_waterlines:
        stern_idx = stern_waterlines.index(wl)
        stern_offset = stern_offsets[stern_idx]
    
    all_stations = list(stations_main)
    all_offsets = list(offsets_main)
    
    if stem_offset is not None:
        stem_x = 244.2 + stem_offset
        stem_y = 0.0
        all_stations.append(stem_x)
        all_offsets.append(stem_y)

    
    if stern_offset is not None:
        stern_x = 0.0 + stern_offset
        stern_y = 0.0
        all_stations.insert(0, stern_x)
        all_offsets.insert(0, stern_y)
    
    all_stations = np.array(all_stations)
    all_offsets = np.array(all_offsets)
    
    weights = np.ones(len(all_stations))
    
    if stern_offset is not None:
        weights[0] = 100
    if stem_offset is not None:
        weights[-1] = 500
    
    for i, offset in enumerate(all_offsets):
        if abs(offset - 18.15) < 0.01:
            weights[i] = 100
    
    # Use fixed degree 7 polynomial fit, fallback to degree 6 if needed
    max_degree = 7
    try:
        coeffs = np.polyfit(all_stations, all_offsets, max_degree, w=weights)
    except Exception:
        coeffs = np.polyfit(all_stations, all_offsets, 6)

    poly = np.poly1d(coeffs)
    
    y_pred = poly(all_stations)
    ss_res = np.sum((all_offsets - y_pred) ** 2)
    ss_tot = np.sum((all_offsets - np.mean(all_offsets)) ** 2)
    
    equation_str = f"y = "
    for i, coef in enumerate(coeffs):
        power = len(coeffs) - 1 - i
        if i == 0:
            equation_str += f"{coef:.10e}"
        else:
            sign = "+" if coef >= 0 else ""
            equation_str += f" {sign} {coef:.10e}"
        
        if power > 0:
            equation_str += f"*x^{power}"
    
    fitted_equations[str(wl)] = {
        'coefficients': coeffs.tolist(),
        'equation': equation_str,
        'waterline': float(wl),
        'includes_stem': stem_offset is not None,
        'includes_stern': stern_offset is not None,
        'num_points': len(all_stations)
    }
    
    x_smooth = np.linspace(all_stations.min(), all_stations.max(), 300)
    y_smooth = poly(x_smooth)
    
    ax.plot(x_smooth, y_smooth, color=colors[idx], linewidth=1.8, label=f'WL {wl}', alpha=0.8)
    
    ax.scatter(stations_main, offsets_main, color=colors[idx], s=30, zorder=3, alpha=0.6)
    
    if stern_offset is not None:
        ax.scatter([all_stations[0]], [all_offsets[0]], color=colors[idx], s=80, 
                   marker='^', zorder=4, edgecolors='black', linewidths=1)
    
    if stem_offset is not None:
        ax.scatter([all_stations[-1]], [all_offsets[-1]], color=colors[idx], s=80, 
                   marker='s', zorder=4, edgecolors='black', linewidths=1)

ax.set_xlabel('Station (m)', fontsize=12, fontweight='bold')
ax.set_ylabel('Half-Breadth (m)', fontsize=12, fontweight='bold')
ax.set_title('Waterline Polynomial Fits with Stem and Stern Profiles', fontsize=14, fontweight='bold')
ax.grid(True, alpha=0.3)
ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=9, ncol=1)

plt.tight_layout()
plt.savefig('waterline_fits.png', dpi=300, bbox_inches='tight')
plt.show()

with open(os.path.join(script_dir, 'waterline_equations.json'), 'w') as f:
    json.dump(fitted_equations, f, indent=2)

with open(os.path.join(script_dir, 'waterline_equations.txt'), 'w') as f:
    f.write("=" * 100 + "\n")
    f.write("DEGREE 7 POLYNOMIAL FITS FOR ALL WATERLINES\n")
    f.write("=" * 100 + "\n\n")
    
    for wl_name, data in fitted_equations.items():
        f.write(f"Waterline {wl_name}:\n")
        f.write(f"  Data points: {data['num_points']}\n")
        f.write(f"  Includes stem: {data['includes_stem']}\n")
        f.write(f"  Includes stern: {data['includes_stern']}\n")
        f.write(f"  Equation: {data['equation']}\n")
        f.write(f"  Coefficients:\n")
        for i, coef in enumerate(data['coefficients']):
            power = len(data['coefficients']) - 1 - i
            f.write(f"    a{power} = {coef:.15e}\n")
        f.write("\n" + "-" * 100 + "\n\n")

print("Files saved:")
print("waterline_equations.json")
print("waterline_equations.txt")
print("waterline_fits.png")
