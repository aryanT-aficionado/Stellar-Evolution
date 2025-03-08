# Stellar Evolution Modeling Project

## Overview
This project aims to model the life cycle of stellar bodies by analyzing data from the **Gaia Archive**. The goal is to extract and analyze key stellar parameters (e.g., mass, luminosity, temperature) to simulate stellar evolution.

---

## 🚀 Progress Steps

### ✅ Step 1: Initial Data Exploration
- Explored the **Gaia Archive** to identify useful catalogs and data sets.
- Successfully pulled data using the `astroquery` library.
- Extracted specific columns from the large Gaia dataset:
    - **ra**, 
    - **dec**, 
    - **parallax**, 
    - **parallax_error**, 
    - **phot_g_mean_mag**,  -- Apparent magnitude in G-band
    - **phot_bp_mean_mag**, -- BP-band magnitude
    - **phot_rp_mean_mag**, -- RP-band magnitude
    - **bp_rp**, -- BP - RP color (temperature proxy)
    - **teff_gspphot**, -- Effective temperature (from Gaia Photometry)
    - **luminosity_gspphot**, -- Luminosity (from Gaia Photometry)
    - **radius_gspphot**,-- Stellar radius (from Gaia Photometry)
    - **radial_velocity**, -- Line-of-sight velocity (if available)
    - **ruwe**, -- Renormalized Unit Weight Error (data quality flag)
    - **phot_variable_flag**, -- Star’s variability flag
    - **feh_gspphot** -- Metallicity [Fe/H] (from Gaia Photometry)

| Column               | Description                                                                 | Project Relevance                                                                 |
|----------------------|-----------------------------------------------------------------------------|-----------------------------------------------------------------------------------|
| `ra`, `dec`          | Right Ascension and Declination (sky coordinates in degrees).              | Spatial mapping of stars for cluster identification.                              |
| `parallax`           | Parallax (milliarcseconds). Convert to distance: `distance_pc = 1000/parallax`. | Calculate distances to stars. Critical for luminosity calculations.               |
| `parallax_error`     | Uncertainty in parallax measurement.                                       | Filter stars with reliable distances (`parallax_error/parallax < 0.2`).           |
| `phot_g_mean_mag`    | Apparent magnitude in Gaia’s broad "G" band.                               | Calculate absolute magnitude (`M_G = G - 5 log10(distance) + 5`).                 |
| `phot_bp_mean_mag`   | Blue Photometer (BP) magnitude.                                            | Stellar temperature/color analysis.                                               |
| `phot_rp_mean_mag`   | Red Photometer (RP) magnitude.                                             | Stellar temperature/color analysis.                                               |
| `bp_rp`              | BP - RP color index.                                                       | Proxy for temperature (redder = cooler, bluer = hotter).                          |
| `radial_velocity`    | Line-of-sight velocity (km/s).                                             | Kinematic studies (e.g., Galactic structure, cluster membership).                 |
| `ruwe`               | Renormalized Unit Weight Error (astrometric quality flag).                | Filter high-quality data (`ruwe < 1.4` = reliable positions).                     |
| `teff_gspphot`       | Effective temperature (Kelvin) from Gaia Photometry.                       | Direct input for Hertzsprung-Russell (HR) diagrams and stellar evolution modeling.|
| `phot_variable_flag` | Flag indicating stellar variability.                                       | Identify pulsating stars (e.g., Cepheids, RR Lyrae) for phase classification.     |

---

## 📂 Future Steps
- Clean and preprocess the data.  
- Develop a mathematical model for stellar evolution.  