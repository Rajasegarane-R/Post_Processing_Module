# Post_Processing_Module

Filters a given uniformly sampled raw input signal around a frequency of interest (in [Hz]) and computes its instantaneous modal parameters (damping ratio, damped frequency, and natural frequency) from the resulting filtered analytic signal.

# Pipeline
1. [`stitching`](@ref): pads both ends of the input signal with polynomially interpolated additional points to suppress the edge effects that arise from FFT-based filtering and the Hilbert transform at the signal boundaries.
2. [`filtering`](@ref): filters the frequency of interest (in [Hz]) out of the padded signal using a Butterworth low-pass/band-pass filter, and returns the filtered signal along with its analytic signal.
3. [`analytical_signal`](@ref): computes the analytic signal of a given uniformly sampled real signal; used internally by [`filtering`](@ref).
4. [`Post_Processing`](@ref): ties the above together to obtain the filtered signal and its analytic signal, then computes the instantaneous modal parameters (damping ratio, damped frequency, and natural frequency) using total-variation regularised differentiation (`tvdiff`).
5. [`Final_Result`](@ref) : holds the final results of [`Post_Processing`](@ref) and other necessary results.

# Exports
- [`Post_Processing`](@ref)
- [`Final_Result`](@ref)
