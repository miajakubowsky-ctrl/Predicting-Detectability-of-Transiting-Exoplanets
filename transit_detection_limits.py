import numpy as np
import matplotlib.pyplot as plt
from math import erf

# ----------------- Inputs (USERS MUST EDIT THESE) -----------------
D = 0.127  # telescope aperture in meters
QE = 0.8
exposure_time = 45  # seconds (change to test)
read_noise = 1.5  # electrons per pixel
sky_brightness = 20.29  # mag/arcsec^2
pixel_scale = 0.81  # arcsec/pixel
aperture_radius_pix = 8
dark_current_rate = 0.0005  # e-/pix/sec
full_well = 51.4e3  # electrons
num_caps = 1
vega_zp = 2.42E-09  # erg/m^2/s/Å
eff = 6201.71  # Å
weff = 1253.71  # Å
fwhm_optimistic = 8   # BROAD PSF (optimistic: less concentrated -> fewer saturations)
fwhm_pessimistic = 5  # TIGHT PSF (pessimistic: more concentrated -> more saturation)
magnitudes = np.linspace(6, 15, 300)

# ----------------- Derived -----------------
A = np.pi * (D / 2) ** 2  # collecting area (m^2)
n_pix_aperture = np.pi * aperture_radius_pix**2
dark_current = dark_current_rate * exposure_time  # e-/pix for exposure

# ----------------- PSF fraction function -----------------
def central_pixel_fraction(fwhm_px: float) -> float:
    """Fraction of a 2D Gaussian's total flux that lands in the central pixel
    when the star is centered on that pixel."""
    sigma = fwhm_px / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    a = 0.5 / (np.sqrt(2.0) * sigma)
    return erf(a) ** 2

# ----------------- Optimistic vs Pessimistic (make sure optimistic is BROADER) ----

alpha_opt = central_pixel_fraction(fwhm_optimistic)
alpha_pess = central_pixel_fraction(fwhm_pessimistic)

# ----------------- Photon model -----------------
E_phot = (6.63e-27 * 2.998e10) / (eff * 1e-8)  # erg per photon
PFD = (vega_zp / E_phot) * weff * 10000  # photons/m^2/s for 0-mag star


flux0 = PFD * QE * A * exposure_time            # electrons for a 0-mag star
N_star = flux0 * 10**(-0.4 * magnitudes)         # total star electrons in aperture

# sky
sky_flux_density = PFD * 10**(-0.4 * sky_brightness)
aperture_area_arcsec2 = n_pix_aperture * pixel_scale**2
N_sky = sky_flux_density * QE * A * exposure_time * aperture_area_arcsec2
sky_per_pixel = sky_flux_density * QE * A * exposure_time * (pixel_scale**2)

# noise
read_noise_total_var = (read_noise ** 2) * n_pix_aperture
dc_total = dark_current * n_pix_aperture
noise_total = np.sqrt(N_star + N_sky + dc_total + read_noise_total_var)

# SNR and min detectable depth
snr_single = N_star / noise_total
snr = snr_single * np.sqrt(num_caps)
min_detectable_depth_ppt = (1 / snr) * 1000

# ----------------- Saturation boundaries -----------------
threshold = 0.8 * full_well

def saturation_boundary(alpha):
    """Return the magnitude where peak pixel reaches threshold.
       If sky+dark already exceed threshold, returns -inf (everything saturates)."""
    F_lim = (threshold - (sky_per_pixel + dark_current)) / max(alpha, 1e-12)
    if F_lim <= 0:
        return -np.inf
    else:
        return -2.5 * np.log10(F_lim / max(flux0, 1e-30))

m_sat_opt = saturation_boundary(alpha_opt)   # optimistic boundary (fewer saturations)
m_sat_pess = saturation_boundary(alpha_pess) # pessimistic boundary (more saturations)

# Masks (star is saturated if its magnitude <= boundary)
sat_mask_opt = magnitudes <= m_sat_opt
sat_mask_pess = magnitudes <= m_sat_pess

# ----------------- Plot -----------------
plt.figure(figsize=(10, 6))
plt.plot(magnitudes, min_detectable_depth_ppt, label='1σ Detection Limit')

# Shade optimistic region (definitely saturated)
xmin, xmax = magnitudes.min(), magnitudes.max()
if np.isfinite(m_sat_opt):
    right = np.clip(m_sat_opt, xmin, xmax)
    if right > xmin:
        plt.axvspan(xmin, right, alpha=0.25, color='red', label='Saturated when FWHM = ' + str(fwhm_optimistic))

# Put black Xs for pessimistic saturation (subsampled to ~0.5 mag spacing)
if sat_mask_pess.any():
    # step size in indices for ~0.5 mag
    mag_step = 0.5
    idx_step = int(np.round(mag_step / (magnitudes[1] - magnitudes[0])))
    pess_indices = np.where(sat_mask_pess)[0][::idx_step]

    plt.scatter(magnitudes[pess_indices], min_detectable_depth_ppt[pess_indices],
                marker='x', s=50, color='black', label='Saturated when FWHM = ' + str(fwhm_pessimistic))

plt.xlabel('Star Magnitude')
plt.ylabel('Min Detectable Transit Depth (ppt)')
plt.title(f'Predicted Exoplanet Transit Detection Limits for a 5-Inch Telescope\n'
          f'({exposure_time}s exposures)')
plt.gca().invert_xaxis()
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

