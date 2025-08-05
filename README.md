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

## 3. HR Diagram and Evolutionary Classification

The Hertzsprung-Russell (HR) diagram is the cornerstone of stellar astrophysics — a plot of luminosity versus temperature (or color) that reveals the evolutionary state of stars. In this section, we construct a high-fidelity HR diagram from our cleaned Gaia DR3 sample and use unsupervised clustering to objectively identify major evolutionary phases.

### 3.1 HR Diagram Construction

We begin by visualizing the classical HR diagram using **absolute G-band magnitude** ($M_G$) on the y-axis (inverted, as standard) and **BP − RP color** as a proxy for effective temperature on the x-axis.

To enhance physical interpretability, the first plot colors each star by its **effective temperature** ($T_{\text{eff}}$) from `teff_gspphot`, using a reversed `viridis` colormap so that **hotter (bluer) stars appear in yellow/blue**, and **cooler (redder) stars in deep red**. A logarithmic normalization ensures balanced color distribution across the wide range of stellar temperatures.

![HR Diagram Colored by Temperature](photos/output_hr_init.png)

This temperature-colored view clearly reveals major sequences:
- The **main sequence** diagonal stretching from hot, bright blue stars to cool, faint red dwarfs
- A dense **red giant branch** rising vertically from the lower right
- A prominent **horizontal branch** and **asymptotic giant branch**
- A tight sequence of **white dwarfs** in the lower-left corner

In a second visualization, we highlight **stellar density** in the HR diagram using kernel density estimation (KDE). Each point is colored by local point density, revealing regions of high stellar concentration — particularly along the main sequence and red clump — while suppressing noise from sparse or scattered outliers. This "density map" view helps identify natural groupings in the data, guiding our choice of clustering strategy.

![HR Diagram Density Map](photos/k_means.png)

> 🔭 **Why This Matters**  
> Unlike theoretical HR diagrams, this is a *data-driven* representation of stellar populations within 1 kpc of the Sun. It reflects real observational biases, completeness limits, and Galactic structure — grounding our analysis in empirical reality.

### 3.2 Clustering Methodology

To move from visual inspection to quantitative classification, we applied **K-means clustering** in the $(BP - RP, M_G)$ plane to group stars into distinct evolutionary populations.

#### Why K-means?

While several clustering algorithms were considered — including **Gaussian Mixture Models (GMM)** and **DBSCAN** — we selected **K-means** for the following reasons:

- **Interpretability**: K-means produces compact, spherical clusters ideal for identifying well-separated sequences like the main sequence, giants, and white dwarfs.
- **Scalability**: With ~70k stars, K-means is computationally efficient and deterministic with fixed initialization.
- **Geometric alignment**: The HR diagram features elongated but relatively convex structures, which K-means can approximate well with sufficient clusters.
- **Reproducibility**: Unlike DBSCAN (sensitive to density variations) or GMM (assumes elliptical distributions), K-means offers stable results across runs when $k$ is well-chosen.

We standardized the features and use the **silhouette analysis** to determine the optimal number of clusters. A value of $k = 5$ provided the best balance between intra-cluster cohesion and astrophysical meaning,aligning with major evolutionary stages.

! [K Means](photos/k_means.png)

The resulting clusters were then analyzed for their photometric and physical characteristics to assign evolutionary labels.

### 3.3 Phase Classification Results

The five clusters identified by K-means correspond closely to canonical stellar evolutionary phases. Below is the classification summary:

| Evolutionary Phase | Count | Percentage | Characteristic Features |
|--------------------|-------|------------|--------------------------|
| **Supergiants**     | 18,040 | 25.73% | Highly luminous ($M_G < -1$), cool ($BP - RP > 1.2$), likely massive stars in late stages (e.g., red supergiants or AGB stars) |
| **Main Sequence**   | 16,955 | 24.18% | Forms the diagonal band from upper-left to lower-right; spans $0.2 < BP - RP < 2.5$, representing core hydrogen-burning stars |
| **White Dwarfs**    | 13,799 | 19.68% | Faint ($M_G > 10$) and hot ($BP - RP < 0.5$); concentrated in the lower-left, marking the final stage of low- to intermediate-mass stars |
| **Giants**          | 13,556 | 19.33% | Bright ($-1 < M_G < 3$) and red ($BP - RP > 1.0$); includes red giant branch (RGB) and red clump stars |
| **Subgiants**       | 7,768  | 11.08% | Transition population between main sequence and giants; located in the "hook" region of the HR diagram |

---

## 📂 Future Steps
- Clean and preprocess the data.  
- Develop a mathematical model for stellar evolution.  