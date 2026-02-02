# SUBXII-and-GSUBXII
The Sine-Unit Burr XII (SUBXII) distribution and a GSUBXII GAMLSS regression model with statistical properties, inference, and applications.
# GSUBXII GAMLSS Regression Model

This repository contains the theoretical development and implementation of the
**Sine-Unit Burr XII (SUBXII)** distribution and the **GSUBXII regression model**
within the **GAMLSS (Generalized Additive Models for Location, Scale and Shape)**
framework.

The proposed model is suitable for modeling response variables defined on the
unit interval (0, 1), allowing flexible distributional assumptions and covariate
effects on multiple parameters.

---

## 📌 Main Contributions
- Definition of the **SUBXII distribution**
- Derivation of statistical properties
- Parameter estimation via maximum likelihood
- Proposal of the **GSUBXII regression model under GAMLSS**
- Simulation studies
- Real data applications

---

## 📊 Model Description
The GSUBXII model extends the SUBXII distribution by incorporating covariates
through structured predictors for location, scale, and shape parameters,
following the GAMLSS framework.

This allows greater flexibility compared to classical regression models,
particularly for bounded data.

---

## 💻 Software and Requirements
The implementation is carried out in **R**, using the following main packages:
- `gamlss`
- `gamlss.dist`
- `stats`
- `ggplot2`

Additional packages may be required for simulations and visualization.

---

## ▶️ Usage
The repository includes scripts for:
- Distribution functions
- Parameter estimation
- Simulation studies
- Regression model fitting
- Diagnostic analysis

Detailed examples are provided within the R scripts.

---

## 📁 Repository Structure
