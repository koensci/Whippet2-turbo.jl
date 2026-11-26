# Statistical Methods: Covariate-Aware Beta-Binomial GLM in `whippet-delta.jl` (`Whippet2-turbo`)

This document describes the statistical model implemented in `src/diff_glm.jl` and `bin/whippet-delta.jl` of **Whippet2-turbo** (`v2.0.0-fork.3`), a performance-optimized and statistically upgraded fork of `timbitz/Whippet.jl` (Sterne-Weiler et al., *Molecular Cell* 2018).

---

## 1. Motivation

Upstream `whippet-delta.jl` compares two groups by pooling reads into group-level Beta posteriors and estimating $P(|\Delta\Psi| > \delta)$ via 10,000 Monte Carlo draws. While fast for simple two-group comparisons, this formulation has two limitations in real biological cohorts:

1. **No biological overdispersion modeling:** Pooling reads across replicates treats biological replicates as technical replicates of a single Binomial/Beta process, inflating false-positive rates under sample-to-sample biological variance.
2. **No covariate adjustment:** Confounders such as library strandedness, sequencing batch, donor/patient pairing, or continuous clinical covariates cannot be controlled in a two-group comparison.

`Whippet2-turbo` replaces the default `whippet-delta.jl` engine (preserving the original method under `--legacy`) with a multithreaded **Beta-Binomial generalized linear model (GLM)** featuring **Cox–Reid Empirical Bayes dispersion shrinkage**, **finite-sample $F$-test likelihood-ratio p-values**, and **Benjamini–Hochberg FDR**.

---

## 2. Continuous Beta-Binomial Model

For each splice-graph event $i$ (whitelisted canonical Whippet motifs: `CE`, `AA`, `AD`, `RI`, `AF`, `AL`, `TS`, `TE`) and each valid sample $j \in \{1, \dots, m_i\}$ passing `--min-reads` ($N_{ij} \ge r_{\min}$), let $\Psi_{ij} \in [0, 1]$ be the quantified Percent Spliced In and $N_{ij} > 0$ be `Total_Reads`. Because `whippet-quant.jl` emits fractional read counts (from bias correction and EM multi-mapper assignment), inclusion and exclusion counts are represented continuously:

$$k_{ij} = \Psi_{ij} N_{ij}, \qquad r_{ij} = N_{ij} - k_{ij}.$$

We model the mean inclusion proportion via a logit link,

$$\mu_{ij} = \text{logistic}(x_j^\top \beta_i) = \frac{1}{1 + \exp(-x_j^\top \beta_i)},$$

where $x_j \in \mathbb{R}^p$ is sample $j$'s row of the full-rank design matrix $X$ (supporting categorical factors, continuous numeric predictors, and `*`/`:` interaction terms). Parameterizing Beta-Binomial overdispersion $\phi_i \in (0, 1)$ by log-precision $\theta_i = \ln M_i$ where $M_i = (1 - \phi_i)/\phi_i$ (so $\text{Var}(k_{ij}) = N_{ij}\mu_{ij}(1-\mu_{ij})[1 + (N_{ij} - 1)\phi_i]$), the continuous Beta-Binomial log-likelihood (omitting the constant log-binomial term that cancels across all derivatives and likelihood-ratio tests) is:

$$\ell(\beta_i, M_i) = \sum_{j=1}^{m_i} \Big[\ln\Gamma(k_{ij} + M_i\mu_{ij}) + \ln\Gamma(r_{ij} + M_i(1-\mu_{ij})) - \ln\Gamma(N_{ij} + M_i) + \ln\Gamma(M_i) - \ln\Gamma(M_i\mu_{ij}) - \ln\Gamma(M_i(1-\mu_{ij}))\Big].$$

### 2.1 Coefficient Estimation (Damped Fisher Scoring)
For fixed precision $M_i$, $\hat{\beta}_i(M_i)$ is computed via damped Fisher scoring (IRLS) with score and observed/Fisher weight for observation $j$:

$$u_{ij} = M_i \mu_{ij}(1-\mu_{ij})\Big[\psi(k_{ij} + a_{ij}) - \psi(a_{ij}) - \psi(r_{ij} + b_{ij}) + \psi(b_{ij})\Big],$$

$$w_{ij} = \max\!\left(10^{-12},\; M_i^2 \mu_{ij}^2 (1-\mu_{ij})^2 \Big[\psi_1(a_{ij}) - \psi_1(k_{ij} + a_{ij}) + \psi_1(b_{ij}) - \psi_1(r_{ij} + b_{ij})\Big]\right),$$

where $a_{ij} = M_i\mu_{ij}$, $b_{ij} = M_i(1-\mu_{ij})$, and $\psi, \psi_1$ are the digamma and trigamma functions. A light $L_2$ ridge penalty ($\lambda = 10^{-4}$ on non-intercept terms; $10^{-6}$ on the intercept), maximum step clamping ($\|\Delta\beta\|_\infty \le 5$), Cholesky solves with adaptive diagonal damping, and monotone backtracking line search guarantee finite convergence even under complete separation ($\Psi = 1.0$ vs. $0.0$).

---

## 3. Three-Pass Cox–Reid Empirical Bayes Dispersion Shrinkage

To stabilize dispersion estimation across events with small sample sizes while avoiding overfitting:

1. **Pass 1 — Lockstep Cohort Loading & Per-Event Rank Reduction:**
   Valid samples are filtered per event (`--min-reads`, `--min-samples`). For each event's valid subset $X_{\text{sub}}$, column-pivoted QR decomposition (`qr(A, ColumnNorm())`) retains the intercept first, drops any nuisance covariate columns that become collinear due to per-event missingness (e.g., an absent batch level), and skips events only if the tested contrast itself loses identifiability.

2. **Pass 2 — Raw Cox–Reid Adjusted Profile Likelihood (APL) on Polymorphic Events:**
   Constitutive/boundary events ($\sum_j k_{ij} < 1$, $\sum_j r_{ij} < 1$, pooled $\bar{\Psi}_i \le 0.01$ or $\ge 0.99$, or $\max_j \Psi_{ij} - \min_j \Psi_{ij} < 10^{-4}$) have a flat Beta-Binomial log-likelihood in $M_i$, where the Cox–Reid $-\tfrac{1}{2}\ln\det(X^\top W X)$ penalty alone produces an artificial interior peak near $\phi \approx 0.64$. Such events are excluded from Pass 2 raw dispersion estimation and Pass 3 MAP optimization (receiving their stratum prior $\mu_\theta$ directly). For every polymorphic event with residual degrees of freedom $d_{\text{res}, i} = m_i - p_i > 0$, raw log-precision $\hat{\theta}_i$ maximizes the Cox–Reid APL (McCarthy, Chen & Smyth, *NAR* 2012) over $\theta \in [-4.0, 11.5]$ ($\phi \in [10^{-5}, 0.982]$) via coarse grid + golden-section search:

   $$\text{APL}(\theta) = \ell(\hat{\beta}_i(e^\theta), e^\theta) - \tfrac{1}{2}\ln\det\!\left(X_i^\top W_i(\theta) X_i + 10^{-8} I\right).$$

   For interior optima, local sampling variance $s_i^2 = 1 / (-\text{APL}''(\hat{\theta}_i))$ is estimated by 3-point finite differences.

3. **Event-Class Stratified Empirical Bayes Priors with Finite-df Bias Correction:**
   Because exon-body tandem UTR events (`TS`, `TE`) exhibit substantially higher biological and positional coverage dispersion than splice-junction events (`CE`, `AA`, `AD`, `RI`, `AF`, `AL`), `whippet-delta.jl` fits separate robust normal priors $\theta_i \sim \mathcal{N}(\mu_\theta, \sigma_\theta^2)$ for **Spliced** vs. **Tandem UTR (`TS`/`TE`)** polymorphic events by default (falling back to the combined polymorphic prior if a stratum has $< 20$ interior events; disable via `--no-stratify-eb`). Within each stratum, $\hat{\sigma}_{\text{obs}} = 1.4826 \cdot \text{MAD}(\hat{\theta})$, and because finite residual degrees of freedom bias sample dispersion medians downward ($\hat{\phi} \propto \chi^2_d / d$, shifting $\hat{\theta}$ upward), the prior mean applies Smyth's exact $\log\chi^2$ median bias correction:

   $$\mu_\theta = \text{median}(\hat{\theta}) + \ln\!\left(\frac{\text{median}(\chi^2_{\bar{d}_{\text{res}}})}{\bar{d}_{\text{res}}}\right), \qquad \sigma_\theta^2 = \text{clamp}\!\left(\hat{\sigma}_{\text{obs}}^2 - \text{median}(s^2),\; 0.08,\; 9.0\right),$$

   with effective prior degrees of freedom $d_0 = \text{clamp}(2 / \sigma_\theta^2,\; 1,\; 30)$.

4. **Pass 3 — MAP Dispersion Shrinkage & Finite-Sample $F$-LRT:**
   Each event's shrunk log-precision $\tilde{\theta}_i$ maximizes the penalized APL posterior $\text{APL}(\theta) - (\theta - \mu_\theta)^2 / (2\sigma_\theta^2)$ (yielding $\tilde{\phi}_i = 1 / (1 + e^{\tilde{\theta}_i})$). With precision fixed at $\tilde{M}_i = e^{\tilde{\theta}_i}$, full ($X_{\text{full}}$) and null ($X_{\text{null}}$, dropping the $\text{df}$ tested columns) models are fitted to obtain the likelihood-ratio statistic $\text{LRT}_i = 2[\ell(\hat{\beta}_{\text{full}}, \tilde{M}_i) - \ell(\hat{\beta}_{\text{null}}, \tilde{M}_i)]$. To account for uncertainty in the shrunk dispersion estimate (analogous to `edgeR` quasi-likelihood $F$-tests), p-values are computed from the finite-sample $F$-distribution:

   $$p_i = P\!\left(F_{\text{df},\; d_{\text{res}, i} + d_0} \ge \frac{\text{LRT}_i}{\text{df}}\right),$$

   reducing to $\chi^2_{\text{df}}$ when `--fixed-dispersion` is supplied. Multiple-testing adjustment uses Benjamini–Hochberg `FDR` across all tested events (preserving `NaN` for untested events).

---

## 4. Effect-Size Reporting

- **Two-group (`-a`/`-b`) and explicit 2-level contrasts (`--test var:A_vs_B`):** Covariate-adjusted marginal inclusion means are computed over the observed covariate distribution of samples in groups A and B:
  $$\hat{\Psi}_B = \frac{1}{n_{AB}} \sum_{j \in A \cup B} \text{logistic}(x_{0,j}^\top \hat{\beta}_0), \qquad \hat{\Psi}_A = \frac{1}{n_{AB}} \sum_{j \in A \cup B} \text{logistic}(x_{0,j}^\top \hat{\beta}_0 + \hat{\beta}_{\text{test}}), \qquad \Delta\Psi = \hat{\Psi}_A - \hat{\Psi}_B,$$
  ensuring `DeltaPsi` and `LogOddsRatio` ($\hat{\beta}_{\text{test}}$, with Fisher standard error `SE`) always share identical sign.
- **Continuous predictors (`--test age`):** Reports the average marginal slope $\overline{\partial\Psi / \partial x} = \frac{1}{m_i}\sum_j \hat{\mu}_{ij}(1-\hat{\mu}_{ij})\hat{\beta}_{\text{test}}$ alongside `LogOddsRatio` and `SE`.
- **Multi-level factors ($K > 2$ levels without `_vs_`):** Performs a $(K-1)$-df omnibus LRT.

---

## 5. Calibration Summary

In 10,000-event Beta-Binomial cohort simulations under heterogeneous biological overdispersion ($\phi \sim \text{LogNormal}$, 5 vs. 5 replicates), the Cox–Reid Empirical Bayes Beta-Binomial GLM controlled null Type-I error (**FPR = `0.045`** at $\alpha = 0.05$) and realized Benjamini–Hochberg FDR (**empirical FDR = `0.0484`** at nominal `0.05`), compared to **34–49% empirical FDR** under unshrunk pooled Monte Carlo testing (`--legacy`). Full unit and integration verification is maintained in `test/runtests.jl` (601,418 assertions).

---

## 6. Example Methods Text for Manuscripts

> Alternative splicing quantification and differential splicing analysis were performed using **Whippet2-turbo** (`v2.0.0-fork.3`, https://github.com/koensci/Whippet2-turbo.jl), a multithreaded fork of Whippet (Sterne-Weiler et al., *Molecular Cell* 2018). Reads were quantified with `whippet-quant.jl` (`--biascorrect`), and differential splicing was tested using the modified `whippet-delta.jl` module introduced in Whippet2-turbo—which replaces Whippet's original pooled Beta Monte Carlo comparison with a covariate-aware Beta-Binomial generalized linear model (`--formula "..." --test ...`). Briefly, the GLM models fractional inclusion counts ($k = \Psi \cdot N$) via damped Fisher scoring, estimates event-wise overdispersion via Cox–Reid Adjusted Profile Likelihood with genome-wide Empirical Bayes shrinkage and Smyth finite-df $\log\chi^2$ median bias correction, and evaluates significance using finite-sample $F(\text{df}, d_{\text{res}} + d_0)$ likelihood-ratio tests with Benjamini–Hochberg FDR control (full statistical formulation in `docs/GLM_METHODS.md`).

---

## References

1. Sterne-Weiler T., Weatheritt R.J., Best A.J., Ha K.C.H., Blencowe B.J. (2018). Efficient and accurate quantitative profiling of alternative splicing patterns of any complexity on a laptop. *Molecular Cell*, 72(1):187–200.e6.
2. McCarthy D.J., Chen Y., Smyth G.K. (2012). Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation. *Nucleic Acids Research*, 40(10):4288–4297.
3. Lund S.P., Nettleton D., McCarthy D.J., Smyth G.K. (2012). Detecting differential expression in RNA-sequence data using quasi-likelihood with shrunken dispersion estimates. *Statistical Applications in Genetics and Molecular Biology*, 11(5).
