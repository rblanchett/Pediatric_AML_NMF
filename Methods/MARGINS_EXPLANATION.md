# NMF Group Assignment and Margin Interpretation

## Overview

The assignment of patients to the *k* = 23 NMF model was performed by utlizing the highest coefficient in the NMF coefficient matrix (*H*).

Although this approach provides a discrete NMF group assignment for every patient, it does not directly indicate how strongly a patient is associated with their assigned group relative to other components.

To address this, two complementary measures were considered:

1. **Absolute Margin (`Margin`)** — the difference between the highest and second-highest NMF component coefficients for each patient.
2. **Relative Margin (`Margin_Relative`)** — the same difference normalized by the total NMF component coefficient magnitude for that patient.

These measures help characterize the separation between competing component assignments and the extent to which a patient's NMF representation is dominated by a single component.

Importantly, neither measure represents a statistical confidence interval, assignment probability, or formal measure of classification accuracy.

---
![NMF margin distribution](Figures/Original_23_Group_Assignment_Margin.png)

### Primary group assignment

Each patient was assigned to the NMF component with the largest coefficient in their corresponding column of *H*.

In R, the assignment was calculated using:

```r
max_index <- apply(H, 2, which.max)

NMF_Group <- paste0("nmf", max_index)

NMF_Score <- H[cbind(max_index, seq_len(ncol(H)))]
```

The `NMF_Score` therefore represents the coefficient associated with the patient's assigned NMF component.

Although the model contained 23 components, only 21 components received at least one patient assignment. The remaining two components were retained in the model but were not the highest-scoring component for any patient.

---

## Absolute Margin

The absolute margin quantifies the difference between the highest and second-highest NMF component coefficients for each patient.

**Margin(j) = H₍₁₎ⱼ − H₍₂₎ⱼ**

where:

- H₍₁₎ⱼ is the largest coefficient for patient *j*.
- H₍₂₎ⱼ is the second-largest coefficient for patient *j*.

The calculation was implemented as:

```r
NMF_Score <- apply(H, 2, max)

Second_Score <- apply(
  H,
  2,
  function(x) sort(x, decreasing = TRUE)[2]
)

Margin <- NMF_Score - Second_Score
```

### Example

Consider two patients with the following three largest component coefficients:

| Patient | Highest | Second-highest | Third-highest | Absolute Margin |
|---|---:|---:|---:|---:|
| Patient A | 0.80 | 0.20 | 0.10 | 0.60 |
| Patient B | 0.80 | 0.75 | 0.10 | 0.05 |

Both patients would be assigned to their highest-scoring component.

However, Patient A exhibits a larger difference between the first- and second-ranked components, whereas Patient B has two components with nearly equivalent coefficients.

### Interpretation

**Large absolute margin**

A large positive margin indicates that the coefficient of the assigned component substantially exceeds that of the second-ranked component.

**Small absolute margin**

A small margin indicates that the two highest-scoring components have similar coefficients.

Such patients may have methylation profiles represented by meaningful contributions from multiple NMF components.

### Important limitation

The absolute margin is expressed in the units of the fitted *H* coefficients.

Consequently, its magnitude depends on component scaling and the NMF model parameterization. It should not be interpreted as an intrinsic biological distance or compared indiscriminately across independently fitted NMF models.

---

## Relative Margin

Because absolute margins depend on the magnitude of the component coefficients, an additional normalized metric can be used to describe the difference between the two strongest components relative to the patient's total coefficient magnitude.

The relative margin is defined as:

**Margin_Relative(j) = (H₍₁₎ⱼ − H₍₂₎ⱼ) / Σᵢ Hᵢⱼ**

Equivalently:

**Margin_Relative(j) = Margin(j) / Total H Score(j)**

where the denominator is the sum of the coefficients across all 23 NMF components for patient *j*.

This calculation corresponds to first normalizing each patient's *H* coefficients by their column sum and then subtracting the two largest normalized coefficients.

Example implementation:

```r
Total_H <- colSums(H)

Margin_Relative <- ifelse(
  Total_H > 0,
  Margin / Total_H,
  NA_real_
)
```

### Example

Suppose a patient has the following NMF coefficient characteristics:

- Highest coefficient: 0.80
- Second-highest coefficient: 0.30
- Sum of all 23 coefficients: 2.00

The absolute margin is:

**0.80 − 0.30 = 0.50**

The relative margin is:

**0.50 / 2.00 = 0.25**

Thus, the difference between the two largest component coefficients accounts for 25% of the total coefficient magnitude for that patient.

### Interpretation

**Higher relative margin**

A higher relative margin indicates that the assigned component exceeds the next-strongest component by a greater fraction of the patient's total component-coefficient magnitude.

**Lower relative margin**

A lower relative margin indicates less separation between the two strongest components relative to the total component-coefficient magnitude.

---

## Relationship Between Margins and NMF Group Identity

NMF group assignments were determined by the maximum coefficient in each patient's *H* vector.

Margin calculations did not alter these assignments.

Instead, the margins provide additional information about the numerical distinction between the first- and second-ranked components.

This distinction is particularly important because a discrete group label may give the impression that each patient belongs exclusively to one molecular subgroup.

NMF, however, represents each patient's methylation profile as a non-negative combination of multiple components.

Therefore:

- The assigned NMF group identifies the patient's highest-coefficient component.
- The absolute margin describes its coefficient separation from the nearest competing component.
- The relative margin describes that separation in relation to the patient's overall coefficient magnitude.
- Neither measure indicates that contributions from other NMF components are biologically absent.

Patients assigned to the same NMF group may consequently exhibit substantial variation in margin values and in their broader patterns of component coefficients.

This variation could be informative when evaluating within-group heterogeneity.

- `01_Original_23_Group_Analysis` — Original NMF group assignments and associated analyses.
- `16_NMF_H_UMAP_Interactive_Explorer` — H-based UMAP visualizations and interactive exploration.
- `Figures/` — Curated visualizations associated with the repository.

