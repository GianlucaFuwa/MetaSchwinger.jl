### Functions for determination of continuous charge and its difference for changes in 1 link-value ###

function cont_charge(links::Array{Float64,3})
    charge = 0 
    for s=1:N_s
        for t=1:N_t
            @inbounds charge+=plaquette(links,s,t)
        end
    end
    return imag(charge)/(2*pi)
end;

function cont_charge_full(configs::Array{Float64,4})
    Q = complex(zeros(N_sweeps))
    for i=1:N_sweeps    
        for s=1:N_s
            for t=1:N_t
                @inbounds Q[i]+=plaquette(configs[:,:,:,i],s,t)
            end
        end
    end
    return imag.(Q)./(2*pi)
end;

function dcont_charge(links::Array{Float64,3},s::Int64,t::Int64,dir::Int64,dx::Float64)
    @inline →(n) = mod1(n+1,N_s) 
    @inline ←(n) = mod1(n-1,N_s)  

    if dir == 1
        a = links[s,t,1]+dx + links[→(s),t,2]    - links[s,→(t),1] - links[s,t,2]
        b = links[s,←(t),1] + links[→(s),←(t),2] - links[s,t,1]-dx - links[s,←(t),2]
        c = links[s,t,1]    + links[→(s),t,2]    - links[s,→(t),1] - links[s,t,2]
        d = links[s,←(t),1] + links[→(s),←(t),2] - links[s,t,1]    - links[s,←(t),2]
    elseif dir == 2
        a = links[←(s),t,1] + links[s,t,2]+dx - links[←(s),→(t),1] - links[←(s),t,2]
        b = links[s,t,1]    + links[→(s),t,2] - links[s,→(t),1]    - links[s,t,2]-dx
        c = links[←(s),t,1] + links[s,t,2]    - links[←(s),→(t),1] - links[←(s),t,2]
        d = links[s,t,1]    + links[→(s),t,2] - links[s,→(t),1]    - links[s,t,2]
    end

    (sin(a) + sin(b) - sin(c) - sin(d)) / (2*pi)
end
