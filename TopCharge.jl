### Function for determination of topological charge ###

function topological_charge(configs::Array{Float64,4})
    Q = complex(zeros(N_sweeps))
    for i=1:N_sweeps    
        for s=1:N_s
            for t=1:N_t
                @inbounds Q[i]+=log(plaquette(configs[:,:,:,i],s,t))
            end
        end
    end
    return imag.(Q)./(2*pi)
end;

function top_charge(links::Array{Float64,3})
    Q = complex(0)
    for s=1:N_s
        for t=1:N_t
            @inbounds Q+=log(plaquette(links,s,t))
        end
    end
    return round(Int64,imag(Q)/(2*pi))
end;