
# Sample Size Considerations for Differential Item Functioning

**Abhishek Dey & Tilaye Yeshanew** — Renaissance Learning *(under review)*

Simulation study formalising the data-generating process underlying Winsteps-style DIF detection and characterising how each detection rule behaves as sample size grows. Key results: p-value-based flagging exceeds 10% FPR at n = 750 and continues inflating; Winsteps contrast estimates converge to within 0.13 logits of the true effect by n = 3,750 with FPR at 0%.

---

## Scripts

Run in the order listed.

**`generate_dif_data.py`**: simulates item response data under a Rasch model with known group-specific difficulty shifts. Output serves as ground truth throughout.

**`sample_size_effects_sim_anchor.R`**: Monte Carlo simulation across sample sizes n ∈ [50, 10,000]. Outputs a labelled CSV (~300 MB) with per-item contrast estimates, p-values, and classification flags across replications. Group labels are M/F reflecting the study's original application; these are interchangeable with any binary group variable.

**`dif_script_anchor.R`**:single-sample DIF analysis pipeline implementing the Winsteps-style contrast statistic and Welch t-test.

**`plot_style.py` / `dif_analysis_repo.ipynb`:** reference plots for the simulation results. Plotting code was AI-assisted; statistical analysis and markdown commentary are our own.

---

## Requirements

R: `TAM`, `dplyr`, `tidyr`

Python: `numpy`, `pandas`, `matplotlib`
