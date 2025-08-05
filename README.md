# Statistical Modeling of Stellar Evolution with Gaia Data

> *"The stars are not for the models alone, but for the data that shapes them."*

This project combines astrophysical insight with modern statistical learning to build data-driven models of stellar evolution using high-precision observations from the **Gaia Data Release 3 (DR3)**. By integrating techniques from survival analysis, clustering, and Bayesian inference, we aim to move beyond purely theoretical frameworks and ground our understanding of how stars live and die in empirical reality.

## 1. Introduction and Motivation 🦕

### 1.1 Project Rationale 🧠

Stellar evolution theory has long been guided by physical models derived from stellar structure equations and nuclear physics. While powerful, these models often rely on simplifying assumptions—about convection, mass loss, metallicity, and binary interactions—that are difficult to test observationally at scale.

Enter **Gaia**: a revolutionary space observatory providing astrometric, photometric, and spectroscopic data for over **1.8 billion stars**. With parallaxes, proper motions, and multi-band photometry, Gaia enables us to construct detailed Hertzsprung-Russell (HR) diagrams with unprecedented accuracy and sample size.

This project leverages this wealth of data to:

- **Build empirical models** of stellar evolution stages,
- **Validate theoretical predictions**—like the mass-lifetime relation—using rigorous statistical methods,
- And create a **fully reproducible pipeline** that connects raw Gaia queries to publication-ready figures and inference.

We treat stars not just as points in an HR diagram, but as outcomes of stochastic processes shaped by mass, age, and environment—opening the door to probabilistic modeling in stellar astrophysics.

### 1.2 Core Objectives 💡

The project is structured around four key goals:

1. **Construct a clean, reliable stellar catalog** from Gaia DR3, applying astrometric and photometric quality cuts.
2. **Classify stars into evolutionary phases** (e.g., pre-main-sequence, main sequence, red giant, white dwarf) using unsupervised clustering techniques (e.g., Gaussian Mixture Models, DBSCAN).
3. **Model main-sequence lifetimes** using survival analysis, treating the transition off the main sequence as an "event time" problem.
4. **Validate the theoretical mass-lifetime relation** via Bayesian hierarchical modeling, quantifying uncertainties and deviations from canonical power laws.

Ultimately, this work aims to establish a statistical framework that future studies can extend—making stellar evolution not just a theoretical narrative, but a quantitatively testable science.


---

## 2. Data Acquisition and Preprocessing 📊

### 2.1 Data Collection

This project is built on high-precision astrometric and photometric data from **Gaia Data Release 3 (DR3)**, accessed via the [Gaia Archive](https://gea.esac.esa.int/archive/) using the `astroquery` Python package. 
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

### 2.2 Preprocessing Pipeline

Raw Gaia data requires careful cleaning and transformation before scientific analysis. Our preprocessing pipeline consists of the following steps:

#### 📏 Distance Estimation

We compute distances using the inverse of parallax, with appropriate quality control:

$$
d\ \text{(pc)} = \frac{1000}{\varpi}, \quad \text{where } \varpi > 0
$$

⚠️ *Note:* While more sophisticated Bayesian distance estimators (e.g., [Bailer-Jones 2021](https://ui.adsabs.harvard.edu/abs/2021AJ....161..147B/abstract)) exist, we use the simple inversion here as a baseline, with strict filtering to minimize bias.

#### 🔆 Absolute Magnitude Calculation

We correct apparent G-band magnitude for distance to obtain absolute magnitude:

$$
M_G = G - 5 \log_{10}(d) + 5
$$

where:
- $ G = \text{phot\_g\_mean\_mag} $
- $ d = \text{distance in parsecs} $

This places all stars on a common luminosity scale, enabling direct comparison in the Hertzsprung-Russell diagram.

#### 🧹 Quality Filtering

To ensure data reliability, we applied the following cuts:

1. **Parallax significance**: $ \frac{\varpi}{\sigma_\varpi} > 5 $ → removes high-uncertainty parallaxes
2. **Positive parallax only**: $ \varpi > 0 $
3. **Effective temperature range**: $ 2500\ \text{K} < T_{\text{eff}} < 15000\ \text{K} $ → excludes outliers and poorly fitted stars
4. **Astrometric quality**: $ \text{ruwe} < 1.4 $ → removes sources with poor astrometric fits
5. **Photometric variability**: Excluded stars flagged as `phot_variable_flag != 'NOT_AVAILABLE' AND 'VARIABLE'`

Additionally, we excluded stars with missing `feh_gspphot` or `luminosity_gspphot`, as these are critical for physical modeling.

#### ✅ Final Dataset Summary

| Stage | Star Count |
|------|------------|
| Initial query result | 100,000 |
| After quality cuts | **70,118** |
| Reduction | ~29.9% |

This 30% reduction reflects the trade-off between sample size and data integrity — a deliberate choice to prioritize **measurement reliability** over completeness.


---

## 📂 Future Steps
- Clean and preprocess the data.  
- Develop a mathematical model for stellar evolution.  