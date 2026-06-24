"""
    Post_Processing_Module

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

"""
module Post_Processing_Module

export Post_Processing , Final_Result

using LinearAlgebra
using FFTA, DSP
using Polynomials , NoiseRobustDifferentiation

#  Analytical Signal Function

"""
    analytical_signal(x) -> a

Computes the analytic signal `a(t) = x(t) + j*x_H(t)` for a given discrete-time signal `x(t)`, where `x_H(t)` is the Hilbert transform of `x(t)`.

`x_H(t)` is obtained by phase-shifting the positive-frequency components of `x` by -90 degrees and the negative-frequency components by +90 degrees, carried out in the frequency domain via Fast Fourier Transform (FFT). 

# Arguments
- `x:: Vector{Float64}`: uniformly sampled, real-valued discrete input signal.

# Returns
- `a :: Vector{ComplexF64}`: complex-valued analytic signal, same length as `x`.`abs.(a)` gives the instantaneous amplitude envelope of `x`, and `angle.(a)` gives its instantaneous phase.

"""
function analytical_signal(x)
    N = length(x) ; 
    X = fft(x) ;
    H = zeros(ComplexF64,N) ;
    if iseven(N)
        H[2:Int(N/2)] .= -im ;
        H[(Int(N/2)+2):end] .= im ;  
    else
        H[2:Int((N+1)/2)] .= -im ;
        H[(Int((N+1)/2)+1):end] .= im ;
    end
    X_h = X .* H ;
    x_H = ifft(X_h) ;
    a = complex.(x) .+ (im.*x_H) ; 
    return a 
end



# Filtering Function

"""
    filtering(omega_filt, delta, N_butter_worth, f_s, x) -> x_filter, a_estimate

Passes the input signal `x` through a zero-phase Butterworth filter targeting angular frequency `omega_filt` (in [rad/s]), then returns the filtered signal together with its analytic signal, computed via [`analytical_signal`](@ref).

If the lower edge of the passband, `omega_filt/(2*pi) - delta` (in [Hz]), is less than or equal to zero, a low-pass Butterworth filter of order `N_butter_worth` with cutoff frequency `omega_filt/(2*pi) + delta` (in [Hz]) is used instead. Otherwise, a band-pass Butterworth filter of order `N_butter_worth`, centred on the frequency of interest `omega_filt/(2*pi)` (in [Hz]) with half-bandwidth `delta`, is used.

# Arguments
- `omega_filt :: Real`: angular frequency of interest (in [rad/s]).
- `delta :: Real`: passband half-bandwidth (in [Hz]).
- `N_butter_worth :: Int`: Butterworth filter order.
- `f_s :: Real`: sampling frequency (in [Hz]).
- `x :: Vector{Float64}`: uniformly sampled, real-valued discrete input signal.

# Returns
- `x_filter :: Vector{Float64}`: filtered signal.
- `a_estimate :: Vector{ComplexF64}`: analytic signal of `x_filter`(see [`analytical_signal`](@ref)).

"""
function filtering(omega_filt , delta , N_butter_worth , f_s , x)
    if ( ((omega_filt)/(2*pi)) - delta ) <= 0
        f_lp = ((omega_filt)/(2*pi)) + delta ;
        flt = digitalfilter(Lowpass(f_lp),Butterworth(N_butter_worth),fs = f_s);
    else
        f_bp_1 = ((omega_filt)/(2*pi)) - delta ;
        f_bp_2 = ((omega_filt)/(2*pi)) + delta ;
        flt = digitalfilter(Bandpass(f_bp_1,f_bp_2),Butterworth(N_butter_worth),fs = f_s);
    end

    x_filter = filtfilt(flt,x)   ; #Filtered Signal

    a_estimate = analytical_signal(x_filter)  # Getting Analytical signal

    return x_filter , a_estimate
end


# Stitching Function

"""
    stitching(n_additional, n_points, t_original, x_original) -> x, t, x_add_1, t_add_1, x_add_2, t_add_2

Pad both ends of `x_original` with `n_additional` synthetic points, then mirrors the padded signal about each end.

At the start of the signal, `x_original` is mirror-reflected to form `reverse(x_original)`, and `n_additional` synthetic points (`x_add_1`) are inserted in the gap between `reverse(x_original)` and `x_original`. The original time vector `t_original` is shifted on each side to make room for this insertion; `t_add_1` is the corresponding time vector for `x_add_1`. The same is done at the end of the signal, producing `x_add_2` and `t_add_2`. The final stitched signal is `x = [reverse(x_original); x_add_1; x_original; x_add_2]`, with `t` the matching concatenated time vector.

`x_add_1` and `x_add_2` are computed locally: at each end, the nearest `n_points` samples of `x_original` are mirror-reflected about that boundary, a polynomial of order `2*((2*n_points)-1)` is fit to the resulting `2*n_points` points, and the polynomial is evaluated at `t_add_1` and `t_add_2` to produce `x_add_1` and `x_add_2` respectively.

# Arguments
- `n_additional :: Int`: number of additional points to add at each end of the signal.
- `n_points :: Int`: number of points taken from each end of the mirrored signal for the polynomial fit.
- `t_original :: Vector{Float64}`: original uniformly sampled time vector.
- `x_original :: Vector{Float64}`: uniformly sampled, real-valued discrete input signal.

# Returns
- `x :: Vector{Float64}`: fully stitched signal.
- `t :: Vector{Float64}`: fully stitched time vector.
- `x_add_1 :: Vector{Float64}`: additional points added at the start of the signal.
- `t_add_1 :: Vector{Float64}`: time vector corresponding to `x_add_1`.
- `x_add_2:: Vector{Float64}`: additional points added at the end of the signal.
- `t_add_2 :: Vector{Float64}`: time vector corresponding to `x_add_2`.

!!! note
    The output time vector `t` is longer than `t_original` and is shifted relative to it, so don't assume `t[1] == t_original[1]`.

"""
function stitching(n_additional,n_points,t_original,x_original)
    dt = t_original[2] - t_original[1] ; 
    order = 2*((2*n_points)-1) ;
    t_interpol = t_original[1:n_points] .+ (((n_additional+1).*dt)./2) ;
    t_interpol = vcat( -reverse(t_interpol) , t_interpol );
    x_interpol = x_original[1:n_points] ;
    x_interpol = vcat( reverse(x_interpol) , x_interpol ) ;
    fit = Polynomials.fit(t_interpol,x_interpol,order) ;

    t_add_1 = range( -(((n_additional+1).*dt)./2) .+ dt , stop = (((n_additional+1).*dt)./2) .- dt , length = n_additional ) ;
    x_add_1 = fit.(t_add_1);

    x_interpol = x_original[end-n_points+1:end] ;
    x_interpol = vcat(x_interpol,reverse(x_interpol));
    fit = Polynomials.fit(t_interpol,x_interpol,order) ;

    t_add_2 = t_add_1 .+ t_original[end] .+ ((n_additional+1).*dt) ;
    x_add_2 = fit.(t_add_1);


    x = vcat(reverse(x_original),x_add_1,x_original,x_add_2) ;
    t = vcat( -reverse(t_original .+ (((n_additional+1).*dt)./2)) , t_add_1 , (t_original .+ (((n_additional+1).*dt)./2)) , t_add_2  ) ;

    return x , t , x_add_1 , t_add_1 , x_add_2 , t_add_2
end



# Final Results

"""
     Final_Result

This structure holds the final results of [`Post_Processing`](@ref) function.

# Fields
- `x_filter :: Vector{Float64}` : filtered signal, same length as `x`.
- `a_filter :: Vector{ComplexF64}` : analytic signal of `x_filter`.
- `zeta :: Vector{Float64}` : instantaneous damping ratio.
- `omega_d :: Vector{Float64}` : instantaneous damped frequency (in [rad/s]).
- `omega_n :: Vector{Float64}` : instantaneous natural frequency (in [rad/s]).
- `x_stitch :: Vector{Float64}` : fully stitched signal
- `t_stitch :: Vector{Float64}` : fully stitched time vector
- `x_s_add_1 :: Vector{Float64}`: additional points added at the start of the signal.
- `t_s_add_1 :: Vector{Float64}`: time vector corresponding to `x_s_add_1`.
- `x_s_add_2 :: Vector{Float64}`: additional points added at the end of the signal.
- `t_s_add_2 :: Vector{Float64}`: time vector corresponding to `x_s_add_2`.
- `x_s_filt :: Vector{Float64}`: filtered stitched signal.
- `a_s_filt :: Vector{ComplexF64}`: analytic signal of filtered stitched signal.
"""
struct Final_Result
    x_filter :: Vector{Float64}
    a_filter :: Vector{ComplexF64}
    zeta :: Vector{Float64}
    omega_d :: Vector{Float64}
    omega_n :: Vector{Float64}
    x_stitch :: Vector{Float64}
    t_stitch :: Vector{Float64}
    x_s_add_1 :: Vector{Float64}
    t_s_add_1 :: Vector{Float64}
    x_s_add_2:: Vector{Float64}
    t_s_add_2 :: Vector{Float64}
    x_s_filt :: Vector{Float64}
    a_s_filt :: Vector{ComplexF64}
end


# Main Function

"""
    Post_Processing(x, t, N_Butter_worth, filt_bw, f_filt, n_add, n_interpol, diff_iter, diff_alpha) -> Final_Result

Filters the uniformly sampled raw input signal `x` around the frequency of interest `f_filt` to obtain a filtered signal `x_filter` and its analytic signal `a_filter`, then estimates the instantaneous damping ratio `zeta`, damped natural frequency `omega_d`, and undamped natural frequency `omega_n` from `a_filter`.

# Method
1. From the time vector `t`, compute the sampling frequency `f_s` and signal length `N`. Convert the filter frequency of interest `f_filt` (in [Hz]) to angular frequency `Omega_filt` (in [rad/s]).
2. Pad `x` with `n_add` synthetic points at each end, built from the nearest `n_interpol` samples, using [`stitching`](@ref), to get the stitched signal `x_stitch`.
3. Band-pass/low-pass filter `x_stitch` around `Omega_filt` using [`filtering`](@ref), giving the filtered stitched signal `x_s_filt` and its analytic signal `a_s_filt`. Crop the padding back off to get `x_filter` and `a_filter`.
4. From `a_filter`, compute the instantaneous amplitude `A` and unwrapped instantaneous phase `phi`.
5. Form `P = log(A) + i*phi` and differentiate it with total-variation regularised differentiation (`tvdiff`) to get `P_dot`. 
     For a signal `x(t) = A0*exp(-zeta*omega_n*t)*exp(j*omega_d*t)`,
     `P_dot = -zeta*omega_n + i*omega_d`, so:
     - `omega_d = imag(P_dot)`
     - `omega_n = abs(P_dot)`
     - `zeta = -cos(angle(P_dot))`

# Arguments
- `x :: Vector{Float64}`: uniformly sampled raw input signal.
- `t :: Vector{Float64}`: uniformly sampled time vector.
- `N_Butter_worth:: Int`: Butterworth filter order.
- `filt_bw:: Real`: filter passband bandwidth (in [Hz]).
- `f_filt:: Real`: filter frequency of interest (in [Hz]).
- `n_add:: Int`: number of additional points to add at each end of `x`.
- `n_interpol:: Int`: number of points taken from each end of the mirrored signal for the polynomial fit.
- `diff_iter:: Int`: number of iterations for total-variation regularised differentiation.
- `diff_alpha:: Real`: regularisation parameter for total-variation regularised differentiation.

# Returns
[`Final_Result`](@ref) : an immutable structure that holds the following results.
- `x_filter :: Vector{Float64}` : filtered signal, same length as `x`.
- `a_filter :: Vector{ComplexF64}` : analytic signal of `x_filter`.
- `zeta :: Vector{Float64}` : instantaneous damping ratio.
- `omega_d :: Vector{Float64}` : instantaneous damped frequency (in [rad/s]).
- `omega_n :: Vector{Float64}` : instantaneous natural frequency (in [rad/s]).
- `x_stitch :: Vector{Float64}` : fully stitched signal
- `t_stitch :: Vector{Float64}` : fully stitched time vector
- `x_s_add_1 ::Vector{Float64}`: additional points added at the start of the signal.
- `t_s_add_1 :: Vector{Float64}`: time vector corresponding to `x_s_add_1`.
- `x_s_add_2 :: Vector{Float64}`: additional points added at the end of the signal.
- `t_s_add_2 :: Vector{Float64}`: time vector corresponding to `x_s_add_2`.
- `x_s_filt :: Vector{Float64}`: filtered stitched signal.
- `a_s_filt :: Vector{ComplexF64}`: analytic signal of filtered stitched signal.

# Example
```julia
     result = Post_Processing(x, t, 1, 0.2, 1.5, 1, 1, 100, 1e-05);
     result.x_filter ;
     result.a_filter ;
     result.zeta ;
     result.omega_d ;
     result.omega_n ;
```

"""
function Post_Processing(x,t,N_Butter_worth,filt_bw,f_filt,n_add,n_interpol,diff_iter,diff_alpha)

    dt = t[2] - t[1] ; # 
    f_s = 1/dt ; #Samping Frequency
    N = length(x) ;
    Omega_filt = f_filt*(2*pi) ;

    x_stitch , t_stitch , x_s_add_1 , t_s_add_1 , x_s_add_2 , t_s_add_2 = stitching(n_add , n_interpol , t , x) ;
    x_s_filt , a_s_filt = filtering(Omega_filt , (filt_bw/2) , N_Butter_worth , f_s , x_stitch );

    x_filter = x_s_filt[N+n_add+1 : end-n_add] ;
    a_filter = a_s_filt[N+n_add+1 : end-n_add] ;

    A = abs.(a_filter) ;
    phi = unwrap(angle.(a_filter)) ;

    P = log.(A) .+ (im.*phi) ;
    P_dot = tvdiff(P , diff_iter , diff_alpha ; dx = dt) ;

    zeta = -cos.(angle.(P_dot)) ;
    omega_d = imag.(P_dot) ;
    omega_n = abs.(P_dot) ;

    return Final_Result(x_filter,a_filter,zeta,omega_d,omega_n,x_stitch,t_stitch,x_s_add_1,t_s_add_1,x_s_add_2,t_s_add_2,x_s_filt,a_s_filt)
    
end

end
