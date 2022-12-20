using Plots
using LaTeXStrings

bias_plot = plot(q_vals, bias_potential, title=L"Bias-Potential \\ $\beta$ = 7.2 \\ 24 x 24", label=L"$V(Q_{cont})$", linewidth=1)
xlabel!(L"$Q_{cont}$")
ylabel!(L"$V(Q_{cont})$")
#savefig(bias_plot,"bias_pot.pdf")

cont_plot = plot(1:N_sweeps,Q_cont_meta, title=L"$Q_{cont}$ \\ $\beta$ = 7.2 \\ 24 x 24", label=L"$Q_{cont}$ + Meta", linewidth=1)
#plot!(cont_plot,1:N_sweeps,Q_cont,label=L"$Q_{cont}$", linewidth=1)
xlabel!("iteration")
ylabel!(L"$Q_{cont}$")
#savefig(cont_plot,"q_cont_compare.pdf")

top_plot = plot(1:N_sweeps,Q_top_meta, title=L"$Q_{top}$ \\ $\beta$ = 7.2 \\ 24 x 24", label=L"$Q_{top}$ + Meta", linewidth=1)
#plot!(top_plot,1:N_sweeps,Q_top,label=L"$Q_{top}$", linewidth=1)
xlabel!("iteration")
ylabel!(L"$Q_{top}$")
#savefig(top_plot,"q_top_compare.pdf")