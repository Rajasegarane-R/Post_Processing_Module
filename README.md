# Post_Processing_Module

This module extracts modal parameters — damping ratio, damped frequency, and natural frequency — from a given vibration signal using the **analytic signal (Hilbert transform)** of the input signal.
---
## What It Does

Given a raw, uniformly sampled vibration signal, the module:

1. **Pads** both ends of the signal with polynomially interpolated points to suppress edge effects from FFT-based filtering and the Hilbert transform (`stitching`).
2. **Filters** the padded signal around the frequency of interest using a zero-phase Butterworth band-pass or low-pass filter, and computes its analytic signal (`filtering`, `analytical_signal`).
3. **Extracts** the instantaneous modal parameters by using total-variation regularised differentiation (`Post_Processing`).

The method is based on writing the analytic signal in a log-polar form:
```math
P(t) = log(A(t)) + j\phi(t)
```
where $A(t)$ is the instantaneous amplitude and $\phi(t)$ is the unwrapped instantaneous phase. 

Incase of a ring-down signal $x(t) = A_0 e^{-\zeta\omega_{n}t} e^{j\omega_{d}t}$, the time derivative satisfies:
```math
\dot{P} = -\zeta\omega_{n} + j\omega_d
```
Using `tvdiff`, we get :
- $\omega_d = imag(\dot{P}$) — instantaneous damped natural frequency
- $\omega_n = abs(\dot{P})$  — instantaneous undamped natural frequency
- $\zeta = -cos(\angle(\dot{P}))$ — instantaneous damping ratio

---

## Dependencies

| Package | Purpose |
|---|---|
| `FFTA.jl` | Fast Fourier Transform to compute Hilbert transform|
| `DSP.jl` | Butterworth filter design and zero-phase filtering |
| `Polynomials.jl` | Polynomial fitting for edge padding |
| `NoiseRobustDifferentiation.jl` | Total-variation regularised differentiation (`tvdiff`) |

---

## Installation

**1. Install Julia** if you have not already.

**2. Install the required packages.** Open Julia REPL and run:

```julia
import Pkg
Pkg.add(["FFTA", "DSP", "Polynomials", "NoiseRobustDifferentiation"])
```
Alternatively, if you have the `Project.toml` from this repository:
```julia
import Pkg
Pkg.activate(".")      # activate the project environment
Pkg.instantiate()      # install all dependencies in Project.toml
```
---
## Usage

### Loading the Module
```julia
include("Post_Processing_Module.jl")
using.Post_Processing_Module
```

### Example

```julia
result = Post_Processing(x, t, N_Butter_worth, filt_bw, f_filt, n_add, n_interpol, diff_iter, diff_alpha)
```

### Parameters

| Parameter | Type | Description |
|---|---|---|
| `x` | `Vector{Float64}` | Uniformly sampled raw input signal |
| `t` | `Vector{Float64}` | Uniformly sampled time vector |
| `N_Butter_worth` | `Int` | Butterworth filter order |
| `filt_bw` | `Real` | Filter passband bandwidth (in [Hz]) |
| `f_filt` | `Real` | Filter frequency of interest (in [Hz]) |
| `n_add` | `Int` | Number of additional padding points to add at each end of `x` |
| `n_interpol` | `Int` | Number of additional points to add at each end of `x` |
| `diff_iter` | `Int` |  Number of iterations for total-variation regularised differentiation |
| `diff_alpha` | `Real` | Regularisation parameter for total-variation regularised differentiation |

### Return Value — `Final_Result`

`Post_Processing` returns a struct named `Final_Result` with the following fields:

| Field | Type | Description |
|---|---|---|
|`x_filter` |`Vector{Float64}`|Filtered signal, same length as `x`|
| `a_filter` |`Vector{ComplexF64}`|Analytic signal of `x_filter`|
| `zeta` |`Vector{Float64}`|Instantaneous damping ratio|
| `omega_d`|`Vector{Float64}`|Instantaneous damped frequency (in [rad/s])|
| `omega_n`|`Vector{Float64}`|Instantaneous natural frequency (in [rad/s])|
| `x_stitch`|`Vector{Float64}`|Fully stitched signal|
| `t_stitch`|`Vector{Float64}`|Fully stitched time vector|
| `x_s_add_1`|`Vector{Float64}`|Additional points added at the start of the signal|
| `t_s_add_1`|`Vector{Float64}`|Time vector corresponding to `x_s_add_1`|
|`x_s_add_2`|`Vector{Float64}`|Additional points added at the end of the signal|
| `t_s_add_2`|`Vector{Float64}`|Time vector corresponding to `x_s_add_2`|
| `x_s_filt`|`Vector{Float64}`|Filtered stitched signal|
| `a_s_filt`|`Vector{ComplexF64}`|Analytic signal of filtered stitched signal|

---

## In-REPL Help

Once the module is loaded, full documentation for any function in this module can be accessed via `?`:

```julia
?Post_Processing_Module.Post_Processing
?Post_Processing_Module.Final_Result
?Post_Processing_Module.stitching
?Post_Processing_Module.filtering
?Post_Processing_Module.analytical_signal
```
---

## Test Case Code

A worked example script `Post_Processing_Test.jl` is provided. It demonstrates the full workflow on a synthetic two-mode free-decay signal with added Gaussian noise:

```math
x(t) = a_{1}e^{-\zeta_{1}\omega_{1}t}cos(\omega_{1}t) + a_{2}e^{-\zeta_{2}\omega_{2}t}cos(\omega_{2}t) + \text{noise}
```

Running the script produces four figures:

| Figure | Description |
|---|---|
| Frequency Response | FFT magnitude and phase of raw vs filtered signal |
| Time Series | Raw signal, filtered signal, and instantaneous amplitude envelope |
| Stitched Signal | How the edge padding looks before and after stitching |
| Modal Results | Instantaneous `zeta`, `omega_d`, `omega_n` vs displacement amplitude |

### Additional Dependencies for the Example

The example script requires these additional packages on top of the core dependencies if you are installing manually :

```julia
import Pkg
Pkg.add(["GLMakie","LaTeXStrings","Distributions"])
```
Alternatively, if you have the `Project.toml` file from this repository, it is already in use.

### Accessing Individual Field Results

```julia
x_f  = results.x_filter   # filtered signal
a_f  = results.a_filter   # analytic signal of x_f
z    = results.zeta        # damping ratio
o_d   = results.omega_d     # damped frequency
o_n   = results.omega_n     # natural frequency
```
