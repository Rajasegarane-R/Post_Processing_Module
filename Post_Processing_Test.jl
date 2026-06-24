include("Post_Processing_Module.jl")
using LinearAlgebra , FFTA
using GLMakie , LaTeXStrings , CairoMakie
using Random , Distributions
using .Post_Processing_Module

GLMakie.closeall() ;
set_theme!(
    theme_latexfonts();
    fontsize = 30,
    linewidth = 2,
    Legend = (
        fontsize = 25,
    )
)

#%%

zeta_1 = 0.5e-2 ; zeta_2 = 0.8e-3 ;
a_1 = 1.5 ; a_2 = -1e-02 ;
omega_1 = 3*pi  ; omega_2 = 8*pi
f_Nyquist = 2*(max(omega_1,omega_2)/(2*pi)) ;

# Filter Control Parameters
dt = (1/f_Nyquist)/1.5;
N_butter_worth = 1 ;
filter_bw = 0.2 ;
f_filt = omega_1/(2*pi) ;
n_additional = 1 ;
n_points = 1 ;
n_iter = 100 ;
alpha = 1e-05 ; 

# Given Input Signal
t_original = 0:dt:100 ;
f_s = 1/dt ;  N = length(t_original);
omega = (0:fld(N,2))*((2*pi*f_s)/N) ;

Random.seed!(4)
N_sig = 2.5e-03  ;
x_original = (a_1.*exp.(-zeta_1.*omega_1.*t_original).*cos.(omega_1.*t_original)) .+ (a_2.*exp.(-zeta_2.*omega_2.*t_original).*cos.(omega_2.*t_original)) .+ rand(Normal(0,N_sig),length(t_original))  ;

#%%
results = Post_Processing_Module.Post_Processing(x_original, t_original,
    N_butter_worth, filter_bw, f_filt,
    n_additional, n_points, n_iter, alpha);

X = fft(x_original) ; 
X_filter = fft(results.x_filter);
#%%
fig_FRF = Figure(size=(1280,720));
fig_t = Figure(size=(1280,720));
fig_s = Figure(size=(1280,720));
fig_results = Figure(size = (720,1280),fontsize=18);

#Plotting FRF
ax_FRF_mag = Axis(fig_FRF[1,1],
    xlabel = L"\text{Frquency} \; \omega\; [\text{rad}\text{s}^{-1}]",
    ylabel = L"\text{Magnitude}\;|X(j\omega)|" ,
    yscale = log10
);
lines!(ax_FRF_mag,omega,abs.(X[1:Int64(fld(N,2))+1]),color= :black,label=L"|X(j\omega)|");
lines!(ax_FRF_mag,omega,abs.(X_filter[1:Int64(fld(N,2))+1]),color=:red,label=L"|X_\text{filter}(j\omega)|",linewidth=3);
axislegend(ax_FRF_mag,position=:rt);

ax_FRF_ph = Axis(fig_FRF[1,2],
    xlabel = L"\text{Frquency} \; \omega\; [\text{rad}\text{s}^{-1}]",
    ylabel = L"\text{Phase Difference}\;\angle(X(j\omega)) \; [^\circ]" ,
);
lines!(ax_FRF_ph,omega,rad2deg.(angle.(X[1:Int64(fld(N,2))+1])),color= :black,label=L"\angle(X(j\omega))");
lines!(ax_FRF_ph,omega,rad2deg.(angle.(X_filter[1:Int64(fld(N,2))+1])),color= :red,label=L"\angle(X_{filter}(j\omega))",linewidth=3);
axislegend(ax_FRF_ph,position=:rt);

Label(fig_FRF[1,1:2,Top()], "Frequency Response",padding = (0, 0, 40, 0),font = :bold)

        
# Plotting time domain data
ax_time = Axis(fig_t[1,1],
    xlabel = L"\text{Time} \; t \; [\text{s}]",
    ylabel = L"\text{Signal}\; x(t)",
    title = "Time Series Data"
)
scatterlines!(ax_time,t_original,x_original,color=:black,label=L"x(t)");
scatterlines!(ax_time,t_original,results.x_filter,color=:blue,label=L"x_\text{filter}(t)");
scatterlines!(ax_time,t_original,abs.(results.a_filter),linewidth=3.5,color=:red,label=L"|a(t)|");
axislegend(ax_time,position=:rt);

        
# Plotting stitched time  domain data
ax_s = Axis(fig_s[1,1],
    xlabel = L"\text{Time} \; t \; [\text{s}]",
    ylabel = L"\text{Signal}\; x(t)",
    title =  "Stitched Time Series Data "
);
scatterlines!(ax_s,vcat(results.t_stitch,(t_original .+ 1.5* (((n_additional+1).*dt)) .+ t_original[end])),vcat(results.x_stitch,reverse(x_original)) ,color=:black,label=nothing );
scatterlines!(ax_s,results.t_stitch,results.x_stitch ,color=:red,label=nothing );
scatterlines!(ax_s,results.t_s_add_1,results.x_s_add_1,color=:blue,label="Additional Points",markersize = 15);
scatterlines!(ax_s,results.t_s_add_2,results.x_s_add_2,color=:blue,label=nothing,markersize = 15);
axislegend(ax_s,position=:rt); 


ax_z = Axis(fig_results[1,1],
    ylabel = L"\textbf{Damping Ratio} \; \zeta",
    # title = "Damping Ratio vs Amplitude"
    xscale = log10
)
scatter!(ax_z,abs.(results.a_filter),results.zeta,color = :black);
vlines!(ax_z,3N_sig)

ax_o_d = Axis(fig_results[2,1],
    ylabel = L"\textbf{Damping Frequency} \; f_d \; [\text{Hz}]",
    #title = "Damping Frequency vs Amplitude"
     xscale = log10
) 
scatter!(ax_o_d,abs.(results.a_filter),(results.omega_d)/(2*pi),color = :black);
vlines!(ax_o_d,3N_sig)

ax_o_n = Axis(fig_results[3,1],
    xlabel = L"\textbf{Displacement Amplitude} \; a \; [\text{m}]",
    ylabel = L"\textbf{Natural Frequency} \; f_n \; [\text{Hz}]",
    #title = "Natural Frequency vs Amplitude"
    xscale = log10
)
scatter!(ax_o_n,abs.(results.a_filter),(results.omega_n)/(2*pi),color = :black);
vlines!(ax_o_n,3N_sig)
linkxaxes!(ax_z,ax_o_d,ax_o_n)
Label(fig_results[1,1,Top()], "Modal Results",padding = (0, 0, 20, 0),font = :bold,fontsize = 25)

#%%
display(GLMakie.Screen(),fig_FRF);
display(GLMakie.Screen(),fig_t);
display(GLMakie.Screen(),fig_s);
display(GLMakie.Screen(),fig_results);
