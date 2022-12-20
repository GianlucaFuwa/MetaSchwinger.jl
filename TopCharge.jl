### Functions for determination of topological charge ###

function topological_charge(configs)
    N_s,N_t,dir,N_prod = size(configs)
    Q = complex(zeros(N_prod))
    for i=1:N_prod    
        for s=1:N_s
            for t=1:N_t
                Q[i]+=log(plaquette(configs[:,:,:,i],s,t))
            end
        end
    end
    return round.(Int,imag.(Q)./(2*pi))
end;
