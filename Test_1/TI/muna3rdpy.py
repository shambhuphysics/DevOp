import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, root_scalar

def birch_murnaghan(V, V0, B0, B0_prime): return (3/2) * B0 * ((V0/V)**(7/3) - (V0/V)**(5/3)) * (1 + (3/4) * (B0_prime - 4) * ((V0/V)**(2/3) - 1))

def bulk_modulus(V, V0, B0, B0_prime): return -V * (-(V0**(5/3) * B0 * (5 * V**(2/3) - 7 * V0**(2/3))) / (2 * V**(10/3)))

V, P_s, P_l = np.loadtxt('VP.dat', comments='#', unpack=True)
TARGET_PRESSURE = 0.5

# Initial estimates
V0_init = V.max() + (V.max() - V.min()) * 0.2
print(f"Birch-Murnaghan EOS Fitting - Target: {TARGET_PRESSURE} GPa")
print(f"Initial estimates for solid: V0={int(V0_init)}, B0=50, B0'=4.5")
print(f"Initial estimates for liquid: V0={int(V0_init)}, B0=50, B0'=4.5")

# Fit both phases
solid_params, _ = curve_fit(birch_murnaghan, V, P_s, p0=[V0_init, 50, 4.5])
liquid_params, _ = curve_fit(birch_murnaghan, V, P_l, p0=[V0_init, 50, 4.5])

# Find volumes at target pressure
vol_s = root_scalar(lambda v: birch_murnaghan(v, *solid_params) - TARGET_PRESSURE, bracket=[V.min(), V.max()*1.5]).root
vol_l = root_scalar(lambda v: birch_murnaghan(v, *liquid_params) - TARGET_PRESSURE, bracket=[V.min(), V.max()*1.5]).root

# Calculate bulk moduli and print results
bulk_s = bulk_modulus(vol_s, *solid_params)
bulk_l = bulk_modulus(vol_l, *liquid_params)
print(f"Solid Phase:\n  V0={solid_params[0]:.2f}, B0={solid_params[1]:.2f}, B0'={solid_params[2]:.2f}\n  Volume at {TARGET_PRESSURE} GPa: {vol_s:.2f} Å³\n  Bulk Modulus: {bulk_s:.2f} GPa\n")
print(f"Liquid Phase:\n  V0={liquid_params[0]:.2f}, B0={liquid_params[1]:.2f}, B0'={liquid_params[2]:.2f}\n  Volume at {TARGET_PRESSURE} GPa: {vol_l:.2f} Å³\n  Bulk Modulus: {bulk_l:.2f} GPa")

# Modern seaborn plotting
sns.set_style('whitegrid'); plt.figure(figsize=(10,6))
V_fit = np.linspace(V.min(), V.max(), 200)
plt.scatter(V, P_s, color='red', s=60, label='Solid Data', zorder=3); plt.scatter(V, P_l, color='blue', s=60, label='Liquid Data', zorder=3)
plt.plot(V_fit, birch_murnaghan(V_fit, *solid_params), 'r:', linewidth=2, label=f'Solid Fit'); plt.plot(V_fit, birch_murnaghan(V_fit, *liquid_params), 'b:', linewidth=2, label=f'Liquid Fit')
plt.xlabel('Volume (Å³)'); plt.ylabel('Pressure (GPa)'); plt.title('Birch-Murnaghan EOS Fitting'); plt.legend(); plt.tight_layout(); plt.show()
