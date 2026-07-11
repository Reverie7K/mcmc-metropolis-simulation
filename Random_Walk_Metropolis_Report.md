# Random Walk Metropolis Algorithm for Sampling a Laplace Distribution

**Course:** ST2195 – Programming for Data Science
**Language:** Python 3 · `NumPy` · `SciPy` · `Matplotlib`
**Reproducibility:** Random seed = 42 throughout

---

## Table of Contents

1. [Overview](#1-overview)
2. [The Algorithm](#2-the-algorithm)
3. [Part (a) — Histogram, KDE & Monte Carlo Estimates](#3-part-a--histogram-kde--monte-carlo-estimates-n--10000-s--1)
4. [Part (b) — R-hat Convergence Diagnostic](#4-part-b--r-hat-convergence-diagnostic-n--2000-j--4)
5. [Summary of Findings](#5-summary-of-findings)
6. [References](#6-references)

---

## 1. Overview

This report implements the **Random Walk Metropolis** algorithm — a variant of Markov Chain Monte Carlo (MCMC) — to draw samples from the probability density function

$$f(x) = \frac{1}{2} e^{-|x|}, \qquad x \in \mathbb{R}$$

This is the well-known **Laplace (double exponential) distribution**, with mean 0 and variance 2 (≈ 1.4142 standard deviation). Direct sampling from $f(x)$ is not straightforward, so the Metropolis algorithm instead constructs a Markov chain whose stationary distribution is $f(x)$, producing a sequence $x_0, x_1, \dots, x_N$. The values $x_1, \dots, x_N$ are then used for every estimate that follows.

---

## 2. The Algorithm

**Step 1.** Choose an initial value $x_0$, a chain length $N$, and a proposal standard deviation $s$.

**Step 2.** For $i = 1, \dots, N$:

1. Draw a candidate $x^{\ast} \sim \text{Normal}(x_{i-1}, s)$.
2. Compute the acceptance ratio $r(x^{\ast}, x_{i-1}) = \dfrac{f(x^{\ast})}{f(x_{i-1})}$.
3. Draw $u \sim \text{Uniform}(0, 1)$.
4. If $u < r(x^{\ast}, x_{i-1})$, accept: set $x_i = x^{\ast}$. Otherwise reject: set $x_i = x_{i-1}$.

**Numerical stability tip:** rather than comparing $u < r$ directly — which risks numerical underflow for very small exponentials — the implementation uses the equivalent **log-ratio criterion**:

$$\log u < \log r(x^{\ast}, x_{i-1}) = \log f(x^{\ast}) - \log f(x_{i-1})$$

---

## 3. Part (a) — Histogram, KDE & Monte Carlo Estimates (N = 10,000, s = 1)

Starting from $x_0 = 0$, at each step a candidate $x^{\ast}$ is drawn from a Normal distribution centred at the current value with standard deviation $s = 1$. The generated samples were used to build a histogram and a kernel density estimate (KDE), with the true $f(x)$ overlaid for comparison.

![Histogram, KDE, and true density overlay](assets/fig1_histogram_kde.png)

The histogram and KDE both track the true Laplace density closely, confirming the chain is sampling from the correct target distribution.

### Monte Carlo estimates

| Quantity | Monte Carlo estimate (seed = 42) | Theoretical value |
|---|---|---|
| Mean | 0.01946 | 0.00000 |
| Standard deviation | 1.49280 | 1.41421 |

> **Conclusion (1a):** The sample mean and standard deviation closely approximate their theoretical counterparts. The small discrepancies are expected random variation for a finite sample of $N = 10{,}000$.

---

## 4. Part (b) — R-hat Convergence Diagnostic (N = 2000, J = 4)

The estimates in Part (a) assume the chain has already **converged** to its stationary distribution. To check this formally, the **Gelman–Rubin $\hat{R}$ diagnostic** was used: $J$ independent chains are run from different starting points, and the variation *within* each chain is compared to the variation *between* chains. Values of $\hat{R}$ close to 1 indicate convergence; the conventional threshold is $\hat{R} < 1.05$.

For $J$ chains, each of length $N$:

$$M_j = \frac{1}{N}\sum_{i=1}^{N} x_i^{(j)} \qquad V_j = \frac{1}{N}\sum_{i=1}^{N} \left(x_i^{(j)} - M_j\right)^2$$

$$W = \frac{1}{J}\sum_{j=1}^{J} V_j \qquad M = \frac{1}{J}\sum_{j=1}^{J} M_j \qquad B = \frac{1}{J}\sum_{j=1}^{J} \left(M_j - M\right)^2$$

$$\hat{R} = \sqrt{\frac{B + W}{W}}$$

Four chains ($J = 4$) were run for $N = 2000$ iterations each, starting from $x_0 \in \{0, 1, 2, 3\}$.

| Setting | Configuration | $\hat{R}$ | Convergence check |
|---|---|---|---|
| $s = 0.001$ | $N = 2000,\ J = 4$ | 56.95180 | ❌ No — well above the 1.05 threshold |
| $s = 0.738$ | $N = 2000,\ J = 4$ | 1.0001 | ✅ Yes — comfortably below the threshold |

![R-hat versus proposal step size](assets/fig2_rhat_vs_stepsize.png)

At the very small step size $s = 0.001$, the chain moves in tiny increments and cannot escape the neighbourhood of its starting value within $N = 2000$ iterations — so the four chains remain distinguishable from one another, and $\hat{R}$ is far above the threshold. As $s$ increases, the chains mix far more freely and $\hat{R}$ drops rapidly toward 1.

> **Conclusion (1b):** Convergence requires a sufficiently large proposal step size. Very small $s$ produces slow mixing because the chain barely moves between iterations. For this target distribution, the optimal region for the proposal standard deviation falls approximately within $s \in [0.5, 1.0]$.

---

## 5. Summary of Findings

- The Random Walk Metropolis sampler recovers the Laplace($0$, 1) distribution accurately: Monte Carlo mean and standard deviation both closely match the theoretical values at $N = 10{,}000$.
- Convergence is **highly sensitive to the proposal step size $s$**: a step size that is too small ($s = 0.001$) causes the chain to get "stuck" near its starting value, giving a badly inflated $\hat{R}$ (56.95); a well-chosen step size ($s \approx 0.7$–$1.0$) yields $\hat{R} \approx 1$, indicating good mixing and convergence across independently-initialised chains.
- This highlights a core practical lesson of MCMC methods: the target distribution can be recovered correctly *only if* the proposal mechanism is tuned well enough to let the chain explore the full support of $f(x)$.

---

## 6. References

1. Gelman, A. and Rubin, D.B. (1992). *Inference from iterative simulation using multiple sequences.* Statistical Science, 7(4), 457–472.
2. Gilks, W.R., Richardson, S. and Spiegelhalter, D.J. (1996). *Markov Chain Monte Carlo in Practice.* Chapman and Hall, London.
3. Harris, C.R. et al. (2020). *Array programming with NumPy.* Nature, 585, 357–362.
4. Virtanen, P. et al. (2020). *SciPy 1.0: Fundamental Algorithms for Scientific Computing in Python.* Nature Methods, 17, 261–272.
5. Hunter, J.D. (2007). *Matplotlib: A 2D Graphics Environment.* Computing in Science and Engineering, 9(3), 90–95.
