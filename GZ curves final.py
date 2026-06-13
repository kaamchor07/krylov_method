import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import UnivariateSpline
import os
script_dir = os.path.dirname(__file__)

df = pd.read_csv(os.path.join(script_dir, 'tchebychev_intersection_points.csv'))

waterlines = df['Waterline'].values.astype(float)
station_columns = [col for col in df.columns if col.startswith('Station_')]
num_stations = len(station_columns)

BREADTH = 36.3
DRAFT = 15.07
DEPTH = 20.757
CB = 0.81
LENGTH = 251.19
LWL_HALF = 125.595
N_STATIONS = 12
DL = 20.9325

KB0 = 7.535
KG = 13.838
DISPLACEMENT_VOLUME = 111303.09532

print("=" * 120)
print("SHIP PARTICULARS")
print("=" * 120)
print(f"Length (L):              {LENGTH:.2f} m")
print(f"Breadth (B):             {BREADTH:.2f} m")
print(f"Draft (T):               {DRAFT:.2f} m")
print(f"Depth (D):               {DEPTH:.2f} m")
print(f"Block Coefficient (Cb):  {CB:.4f}")
print(f"Number of Stations (N):  {N_STATIONS}")
print(f"Station Spacing (DL):    {DL:.4f} m")
print(f"Displacement Volume (∇): {DISPLACEMENT_VOLUME:.5f} m³")
print(f"KB₀:                     {KB0:.3f} m")
print(f"KG:                      {KG:.3f} m")
print("=" * 120)

def fit_curve(x_data, y_data, degree=6, num_points=100):
    below_cap = np.abs(np.abs(x_data)) < 18.15
    at_cap = np.abs(np.abs(x_data) - 18.15) < 0.01
    
    x_below = x_data[below_cap]
    y_below = y_data[below_cap]
    x_at_cap = x_data[at_cap]
    y_at_cap = y_data[at_cap]
    
    if len(x_below) >= 2:
        if len(x_below) < degree + 1:
            degree = max(1, len(x_below) - 1)
        
        coeffs = np.polyfit(y_below, x_below, degree)
        poly = np.poly1d(coeffs)
        y_smooth_below = np.linspace(y_below.min(), y_below.max(), num_points)
        x_smooth_below = poly(y_smooth_below)
    else:
        y_smooth_below = y_below
        x_smooth_below = x_below
    
    if len(y_at_cap) > 0:
        sign = np.sign(x_at_cap[0])
        if len(y_smooth_below) > 0:
            y_start = y_smooth_below.max()
        else:
            y_start = y_at_cap.min()
        
        y_smooth_cap = np.linspace(y_start, y_at_cap.max(), num_points//4)
        x_smooth_cap = np.full_like(y_smooth_cap, sign * 18.15)
        
        y_smooth = np.concatenate([y_smooth_below, y_smooth_cap])
        x_smooth = np.concatenate([x_smooth_below, x_smooth_cap])
    else:
        y_smooth = y_smooth_below
        x_smooth = x_smooth_below
    
    return x_smooth, y_smooth

def calculate_intersections(rotation_angle, rotation_center_x, rotation_center_y):
    angle_rad = np.radians(rotation_angle)
    slope = np.tan(angle_rad)
    
    horizontal_line_y = 20.757
    horizontal_line_y_bottom = 0
    vertical_line_x_right = 18.15
    vertical_line_x_left = -18.15
    
    intersection_data = []
    
    for i, station_col in enumerate(station_columns):
        station_num = i + 1
        
        half_breadths = df[station_col].values.copy()
        half_breadths[half_breadths > 18.0] = 18.15 

        valid_indices = ~np.isnan(half_breadths) & (half_breadths >= 0)
        valid_waterlines = waterlines[valid_indices]
        valid_half_breadths = half_breadths[valid_indices]
        
        if len(valid_waterlines) > 0:
            for side_type in ['original', 'mirror']:
                if side_type == 'original':
                    if station_num <= 6:
                        x_coords = -valid_half_breadths
                        side_label = "Forward"
                    else:
                        x_coords = valid_half_breadths
                        side_label = "Aft"
                else:
                    if station_num <= 6:
                        x_coords = valid_half_breadths
                        side_label = "Aft (Mirror)"
                    else:
                        x_coords = -valid_half_breadths
                        side_label = "Forward (Mirror)"
                
                x_smooth, y_smooth = fit_curve(x_coords, valid_waterlines, degree=6, num_points=1000)
                
                a = slope
                b = -1
                c = rotation_center_y - slope * rotation_center_x
                
                distances = np.abs(a * x_smooth + b * y_smooth + c) / np.sqrt(a**2 + b**2)
                min_dist_idx = np.argmin(distances)
                x_intersection = x_smooth[min_dist_idx]
                y_intersection = y_smooth[min_dist_idx]
                
                if distances[min_dist_idx] < 0.5:
                    if abs(x_intersection) > vertical_line_x_right:
                        if x_intersection > 0:
                            x_intersection = vertical_line_x_right
                        else:
                            x_intersection = vertical_line_x_left
                        y_intersection = rotation_center_y + slope * (x_intersection - rotation_center_x)
                        intersection_type = 'Vertical Line (±18.15)'
                    else:
                        intersection_type = 'Curve Intersection'
                    
                    distance_from_center = np.sqrt((x_intersection - rotation_center_x)**2 + 
                                                   (y_intersection - rotation_center_y)**2)
                    
                    intersection_data.append({
                        'Station': station_num,
                        'Side': side_label,
                        'X-Coordinate': x_intersection,
                        'Y-Coordinate': y_intersection,
                        'Distance from Center': distance_from_center,
                        'Type': intersection_type
                    })
                else:
                    avg_x = np.mean(x_smooth)
                    y_on_rotated_at_avg_x = rotation_center_y + slope * (avg_x - rotation_center_x)
                    
                    if y_smooth.max() < y_on_rotated_at_avg_x:
                        x_fallback = (horizontal_line_y - rotation_center_y) / slope + rotation_center_x
                        y_fallback = horizontal_line_y
                        fallback_type = 'Horizontal Line (20.757)'
                    else:
                        x_fallback = (horizontal_line_y_bottom - rotation_center_y) / slope + rotation_center_x
                        y_fallback = horizontal_line_y_bottom
                        fallback_type = 'Horizontal Line (0.0 - Baseline)'
                    
                    distance_from_center = np.sqrt((x_fallback - rotation_center_x)**2 + 
                                                   (y_fallback - rotation_center_y)**2)
                    
                    intersection_data.append({
                        'Station': station_num,
                        'Side': side_label,
                        'X-Coordinate': x_fallback,
                        'Y-Coordinate': y_fallback,
                        'Distance from Center': distance_from_center,
                        'Type': fallback_type
                    })
    
    return intersection_data

def calculate_waterplane_parameters(intersection_data, angle):
    a_values = []
    b_values = []
    
    for intersection in intersection_data:
        distance = abs(intersection['Distance from Center'])
        
        if "Aft" in intersection['Side']:
            a_values.append(distance)
        elif "Forward" in intersection['Side']:
            b_values.append(distance)
    
    sum_a_plus_b_raw = sum(a_values) + sum(b_values)
    sum_a2_minus_b2_raw = sum([a**2 for a in a_values]) - sum([b**2 for b in b_values])
    sum_a3_plus_b3_raw = sum([a**3 for a in a_values]) + sum([b**3 for b in b_values])
    
    if sum_a_plus_b_raw != 0:
        Yf = (sum_a2_minus_b2_raw / sum_a_plus_b_raw) / 2
    else:
        Yf = 0
    
    Ix = DL * sum_a3_plus_b3_raw / 3
    DI = DL * sum_a_plus_b_raw * (Yf ** 2)
    Ix_theta = Ix - DI
    BM_theta = Ix_theta / DISPLACEMENT_VOLUME
    Yf_half = Yf / 2
    
    sum_a_plus_b_integrated = DL * sum_a_plus_b_raw
    sum_a2_minus_b2_integrated = DL * sum_a2_minus_b2_raw
    sum_a3_plus_b3_integrated = DL * sum_a3_plus_b3_raw
    
    print(f"\n  Waterplane Parameter Calculations for {angle}°:")
    print(f"  ├────────────────────────────────────────────────────────────")
    print(f"  I.   Σ(a + b)           = {sum_a_plus_b_raw:.6f}")
    print(f"  II.  Σ(a² - b²)         = {sum_a2_minus_b2_raw:.6f}")
    print(f"  III. Σ(a³ + b³)         = {sum_a3_plus_b3_raw:.6f}")
    print(f"  IV.  Yf = (II/I)/2      = ({sum_a2_minus_b2_raw:.6f}/{sum_a_plus_b_raw:.6f})/2 = {Yf:.6f} m")
    print(f"  V.   Ix = DL*III/3      = {DL:.4f}*{sum_a3_plus_b3_raw:.6f}/3 = {Ix:.6f} m⁴")
    print(f"  VI.  DI = DL*I*IV²      = {DL:.4f}*{sum_a_plus_b_raw:.6f}*{Yf**2:.6f} = {DI:.6f} m⁴")
    print(f"  VII. Ixθ = V - VI       = {Ix:.6f} - {DI:.6f} = {Ix_theta:.6f} m⁴")
    print(f"  VIII. BMθ = VII/∇       = {Ix_theta:.6f}/{DISPLACEMENT_VOLUME:.6f} = {BM_theta:.6f} m")
    print(f"  IX.  Yf/2               = {Yf:.6f}/2 = {Yf_half:.6f} m")
    
    return {
        'a_values': a_values,
        'b_values': b_values,
        'sum_a_plus_b_raw': sum_a_plus_b_raw,
        'sum_a2_minus_b2_raw': sum_a2_minus_b2_raw,
        'sum_a3_plus_b3_raw': sum_a3_plus_b3_raw,
        'Yf': Yf,
        'Ix': Ix,
        'DI': DI,
        'Ix_theta': Ix_theta,
        'BM': BM_theta,
        'Yf_half': Yf_half,
        'sum_a_plus_b_integrated': sum_a_plus_b_integrated,
        'sum_a2_minus_b2_integrated': sum_a2_minus_b2_integrated,
        'sum_a3_plus_b3_integrated': sum_a3_plus_b3_integrated
    }

def calculate_new_cog_pdf_method(Yf, delta_theta_deg, current_angle_deg, current_F_x, current_F_y):
    current_angle_rad = np.radians(current_angle_deg)
    delta_theta_rad = np.radians(delta_theta_deg)
    next_angle_rad = current_angle_rad + delta_theta_rad
    
    d = (Yf / 2) * delta_theta_rad
    
    along_old_wl_x = np.cos(current_angle_rad)
    along_old_wl_y = np.sin(current_angle_rad)
    
    B_x = current_F_x + (Yf / 2) * along_old_wl_x
    B_y = current_F_y + (Yf / 2) * along_old_wl_y
    
    along_new_wl_x = np.cos(next_angle_rad)
    along_new_wl_y = np.sin(next_angle_rad)
    
    point_on_aux_x = current_F_x + Yf * along_new_wl_x
    point_on_aux_y = current_F_y + Yf * along_new_wl_y
    
    vec_B_to_point_x = point_on_aux_x - B_x
    vec_B_to_point_y = point_on_aux_y - B_y
    
    dot_product = (vec_B_to_point_x * along_new_wl_x + 
                   vec_B_to_point_y * along_new_wl_y)
    
    O_x = B_x + dot_product * along_new_wl_x
    O_y = B_y + dot_product * along_new_wl_y
    
    return O_x, O_y, d

angles = [0, 5, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90]
F_x = 0
F_y = DRAFT

all_results = []

print("\n" + "=" * 120)
print("CALCULATING BM FOR EACH ANGLE")
print("=" * 120)

for i, angle in enumerate(angles):
    print(f"\n{'='*120}")
    print(f"ANGLE: {angle}° | Current F: ({F_x:.6f}, {F_y:.6f})")
    print(f"{'='*120}")
    
    intersections = calculate_intersections(angle, F_x, F_y)
    waterplane_params = calculate_waterplane_parameters(intersections, angle)
    
    result = {
        'Angle': angle,
        'F_X': F_x,
        'F_Y': F_y,
        'intersections': intersections,
        'waterplane_params': waterplane_params
    }
    all_results.append(result)
    
    if i < len(angles) - 1:
        next_angle = angles[i + 1]
        delta_theta = next_angle - angle
        F_x_new, F_y_new, d = calculate_new_cog_pdf_method(
            waterplane_params['Yf'],
            delta_theta,
            angle,
            F_x,
            F_y
        )
        
        print(f"\n  Next centroid calculation:")
        print(f"    Yf = {waterplane_params['Yf']:.6f} m")
        print(f"    δθ = {delta_theta}°")
        print(f"    New F for {next_angle}°: ({F_x_new:.6f}, {F_y_new:.6f})")
        
        F_x = F_x_new
        F_y = F_y_new

print("\n" + "=" * 120)
print("KRYLOV METHOD INTEGRATION")
print("=" * 120)

n = len(angles)
theta_deg = np.array([r['Angle'] for r in all_results])
theta_rad = np.deg2rad(theta_deg)

BM = np.array([r['waterplane_params']['BM'] for r in all_results])

cos_theta = np.cos(theta_rad)
sin_theta = np.sin(theta_rad)

BM_cos = BM * cos_theta
BM_sin = BM * sin_theta

delta_theta_rad = np.diff(theta_rad)
delta_theta_half = np.zeros(n)
delta_theta_half[1:] = delta_theta_rad / 2

integral_BM_cos = np.zeros(n)
for i in range(1, n):
    integral_BM_cos[i] = integral_BM_cos[i-1] + BM_cos[i-1] + BM_cos[i]

y_theta = delta_theta_half * integral_BM_cos

integral_BM_sin = np.zeros(n)
for i in range(1, n):
    integral_BM_sin[i] = integral_BM_sin[i-1] + BM_sin[i-1] + BM_sin[i]

Z_theta_minus_Zc = delta_theta_half * integral_BM_sin

y_cos = y_theta * cos_theta
Z_sin = Z_theta_minus_Zc * sin_theta
Ic = y_cos + Z_sin

KB_sin = KB0 * sin_theta
KN = Ic + KB_sin
KG_sin = KG * sin_theta
GZ = KN - KG_sin

for i, result in enumerate(all_results):
    result['theta_rad'] = theta_rad[i]
    result['cos_theta'] = cos_theta[i]
    result['sin_theta'] = sin_theta[i]
    result['BM_cos'] = BM_cos[i]
    result['BM_sin'] = BM_sin[i]
    result['integral_BM_cos'] = integral_BM_cos[i]
    result['integral_BM_sin'] = integral_BM_sin[i]
    result['y_theta'] = y_theta[i]
    result['Z_theta_minus_Zc'] = Z_theta_minus_Zc[i]
    result['y_cos'] = y_cos[i]
    result['Z_sin'] = Z_sin[i]
    result['Ic'] = Ic[i]
    result['KB_sin'] = KB_sin[i]
    result['KN'] = KN[i]
    result['GZ'] = GZ[i]

print("\n" + "=" * 200)
print("KRYLOV METHOD CALCULATION TABLE")
print("=" * 200)
print(f"{'θ°':<6} {'θ(rad)':<10} {'BMθ':<10} {'cosθ':<10} {'BM·cosθ':<12} {'∫BMcosθ':<12} {'sinθ':<10} {'BM·sinθ':<12} {'∫BMsinθ':<12} {'yθcosθ':<12} {'Zθsinθ':<12} {'Ic':<12} {'KN':<10} {'GZ':<10}")
print("-" * 200)

for r in all_results:
    print(f"{r['Angle']:<6.0f} "
          f"{r['theta_rad']:<10.6f} "
          f"{r['waterplane_params']['BM']:<10.4f} "
          f"{r['cos_theta']:<10.6f} "
          f"{r['BM_cos']:<12.6f} "
          f"{r['integral_BM_cos']:<12.6f} "
          f"{r['sin_theta']:<10.6f} "
          f"{r['BM_sin']:<12.6f} "
          f"{r['integral_BM_sin']:<12.6f} "
          f"{r['y_cos']:<12.6f} "
          f"{r['Z_sin']:<12.6f} "
          f"{r['Ic']:<12.6f} "
          f"{r['KN']:<10.4f} "
          f"{r['GZ']:<10.4f}")

print("=" * 200)

BM0 = all_results[0]['waterplane_params']['BM']
KM0 = KB0 + BM0
GM0_actual = KM0 - KG
print(f"\n{'='*120}")
print(f"INITIAL STABILITY PARAMETERS")
print(f"{'='*120}")
print(f"BM₀  = {BM0:.6f} m")
print(f"KB₀  = {KB0:.3f} m")
print(f"KM₀  = KB₀ + BM₀ = {KB0:.3f} + {BM0:.6f} = {KM0:.6f} m")
print(f"KG   = {KG:.3f} m")
print(f"GM₀  = KM₀ - KG = {KM0:.6f} - {KG:.3f} = {GM0_actual:.6f} m")
print(f"Status: {'STABLE (GM > 0)' if GM0_actual > 0 else 'UNSTABLE (GM < 0)'}")
print(f"{'='*120}")

results_df = pd.DataFrame({
    'Heel_Angle_deg': theta_deg,
    'Theta_rad': theta_rad,
    'Sum_a_plus_b': [r['waterplane_params']['sum_a_plus_b_raw'] for r in all_results],
    'Sum_a2_minus_b2': [r['waterplane_params']['sum_a2_minus_b2_raw'] for r in all_results],
    'Sum_a3_plus_b3': [r['waterplane_params']['sum_a3_plus_b3_raw'] for r in all_results],
    'Yf_m': [r['waterplane_params']['Yf'] for r in all_results],
    'Ix_m4': [r['waterplane_params']['Ix'] for r in all_results],
    'DI_m4': [r['waterplane_params']['DI'] for r in all_results],
    'Ix_theta_m4': [r['waterplane_params']['Ix_theta'] for r in all_results],
    'BM_m': BM,
    'cos_theta': cos_theta,
    'sin_theta': sin_theta,
    'BM_cos_theta': BM_cos,
    'BM_sin_theta': BM_sin,
    'Integral_BM_cos': integral_BM_cos,
    'Integral_BM_sin': integral_BM_sin,
    'y_theta_m': y_theta,
    'Z_theta_minus_Zc_m': Z_theta_minus_Zc,
    'y_theta_cos': y_cos,
    'Z_theta_sin': Z_sin,
    'Ic_m': Ic,
    'KB_sin_theta': KB_sin,
    'KN_m': KN,
    'GZ_m': GZ
})

output_filename = os.path.join(script_dir, 'krylov_detailed_results.csv')
results_df.to_csv(output_filename, index=False)
print(f"\nResults exported to: {output_filename}")

# Fit spline to GZ curve
# Use UnivariateSpline for smooth interpolation
spline = UnivariateSpline(theta_deg, GZ, s=0.1, k=3)  # s is smoothing factor, k is degree
theta_smooth = np.linspace(0, 90, 200)
gz_smooth = spline(theta_smooth)

# Also keep polynomial fit for comparison
gz_poly_degree = 4
gz_coeffs = np.polyfit(theta_deg, GZ, gz_poly_degree)
gz_poly = np.poly1d(gz_coeffs)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 8))

# Body Plan (Chebyshev) - Modified to plot all 12 waterlines on both sides
for i, station_col in enumerate(station_columns):
    station_num = i + 1
    half_breadths = df[station_col].values.copy()
    half_breadths[half_breadths > 18.0] = 18.15
    
    valid_indices = ~np.isnan(half_breadths) & (half_breadths >= 0)
    valid_waterlines = waterlines[valid_indices]
    valid_half_breadths = half_breadths[valid_indices]
    
    # Filter out zero values - keep last zero and first non-zero onwards
    if len(valid_half_breadths) > 0:
        # Find first non-zero index
        non_zero_mask = valid_half_breadths > 0
        if np.any(non_zero_mask):
            first_non_zero_idx = np.argmax(non_zero_mask)
            # Include the last zero value (if it exists)
            start_idx = max(0, first_non_zero_idx - 1)
            valid_waterlines = valid_waterlines[start_idx:]
            valid_half_breadths = valid_half_breadths[start_idx:]
    
    if len(valid_waterlines) > 0:
        x_smooth, y_smooth = fit_curve(valid_half_breadths, valid_waterlines, degree=5, num_points=300)

        # --- Removed DRAFT clipping block ---
        
        # Color coding based on station number
        # Stations 1-6: blue on right, red on left
        # Stations 7-12: red on right, blue on left
        if station_num <= 6:
            right_color = 'b'
            left_color = 'r'
        else:
            right_color = 'r'
            left_color = 'b'
        
        # Clip curve so it never crosses the midline (HB = 0)
        right_mask = x_smooth >= 0
        x_smooth_right = x_smooth[right_mask]
        y_smooth_right = y_smooth[right_mask]
        # Plot on positive side (right)
        ax1.plot(x_smooth_right, y_smooth_right, right_color + '-', linewidth=1.5, alpha=0.7)

        # Clip mirrored curve so it never crosses midline (HB = 0)
        left_mask = x_smooth >= 0
        x_smooth_left = -x_smooth[left_mask]
        y_smooth_left = y_smooth[left_mask]
        # Plot mirrored on negative side (left)
        ax1.plot(x_smooth_left, y_smooth_left, left_color + '-', linewidth=1.5, alpha=0.7)

ax1.axhline(y=DRAFT, color='green', linestyle='--', linewidth=2, alpha=0.6, label=f'Design WL ({DRAFT} m)')
ax1.axhline(y=0, color='black', linestyle='-', linewidth=1.5, alpha=0.8, label='Baseline')
ax1.axvline(x=0, color='gray', linestyle=':', linewidth=1.5, alpha=0.6)
ax1.set_xlabel('Half-Breadth (m)', fontsize=12, fontweight='bold')
ax1.set_ylabel('Height (m)', fontsize=12, fontweight='bold')
ax1.set_title('Chebyshev Body Plan (All 12 Waterlines Mirrored)', fontsize=13, fontweight='bold')
ax1.grid(True, alpha=0.3)
ax1.legend(fontsize=10)
ax1.set_xlim([-20, 20])
ax1.set_ylim([0, 22])
ax1.set_aspect('equal')

# GZ Curve with spline fit
ax2.plot(theta_deg, GZ, 'ro', linewidth=2.5, markersize=8, label='Calculated GZ', zorder=3)
ax2.plot(theta_smooth, gz_smooth, 'b-', linewidth=2.5, label='Spline Fit', zorder=2)
ax2.axhline(y=0, color='k', linestyle='--', linewidth=1, alpha=0.5)
ax2.axhline(y=0.20, color='g', linestyle='--', linewidth=1.5, alpha=0.7, label='IMO Min (0.20 m)')
ax2.fill_between(theta_deg, 0, GZ, alpha=0.15, color='red', zorder=1)
ax2.fill_between(theta_smooth, 0, gz_smooth, where=(gz_smooth > 0), alpha=0.1, color='blue', zorder=1)
ax2.grid(True, alpha=0.3)
ax2.set_xlabel('Heel Angle (degrees)', fontsize=12, fontweight='bold')
ax2.set_ylabel('GZ (m)', fontsize=12, fontweight='bold')
ax2.set_title('Righting Arm Curve (GZ)', fontsize=13, fontweight='bold')
ax2.legend(fontsize=10, loc='best')
ax2.set_xlim([0, 90])

plt.tight_layout()
plot_filename = os.path.join(script_dir, 'stability_results.png')
plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
print(f"Plot saved to: {plot_filename}")
plt.show()

print("\n" + "=" * 120)
print("KRYLOV METHOD CALCULATIONS COMPLETED")
print("=" * 120)