### Metropolis-Algorithm involving Metadynamics ###

function sweep_meta!(links::Array{Float64,3},bias::Array{Float64,1})   #Metropolis-Sweep Funktion
    Q = cont_charge(links)
    for s = 1:N_s
        for t = 1:N_t
            for dir=1:d
                dU = (rand()-0.5)*2*ε
                Q_prop = Q + dcont_charge(links,s,t,dir,dU)
                old_ind = grid_ind(Q)
                prop_ind = grid_ind(Q_prop)
                ΔS = daction(links,s,t,dir,dU)  #Berechne Wirkungsunterschied
                ΔV = bias[prop_ind]-bias[old_ind]
                ΔV_pen = penalty_potential(Q_prop)-penalty_potential(Q)
                if rand()<exp(-ΔS-ΔV-ΔV_pen)   #Annahme nach gegebener Wahrscheinlichkeit
                    links[s,t,dir] += dU
                    Q = Q_prop
                    global acc += 1
                    update_bias!(bias,Q)
                end
            end
        end
    end
    return links
end;

function metropolis_meta(start::Array{Float64,3},bias::Array{Float64,1})
    links = deepcopy(start)
    Q_top = zeros(Int64,Int(N_sweeps/n_skip),1)
    Q_cont = zeros(Float64,Int(N_sweeps/n_skip),1)
    Q_squared_weighted = zeros(Float64,Int(N_sweeps/n_skip),1)
    denom = 0
    for l = 1:N_therm
        links = sweep!(links)
    end
    i = 1
    for l = 1:N_sweeps  #Führe Simulation durch
        links = sweep_meta!(links,bias)
        if l%n_skip == 0
            Qc = cont_charge(links) 
            @inbounds Q_top[i] = top_charge(links)
            @inbounds Q_cont[i] = Qc
            @inbounds Q_squared_weighted[i] = Qc^2*exp(-bias[grid_ind(Qc)]-penalty_potential(Qc)) 
            denom += exp(-bias[grid_ind(Qc)]-penalty_potential(Qc))
            i += 1 
        end
    end
    display(denom)
    Q_squared_weighted = Q_squared_weighted./denom
    println("acc. rate = $(acc/(N_sweeps*N_s*N_t*d))")
    return Q_top,Q_cont,Q_squared_weighted
end;

### Metropolis-Algorithm without Metadynamics ###

function sweep!(links::Array{Float64,3})   #Metropolis-Sweep Funktion
    for x = 1:N_s
        for t = 1:N_t
            for dir=1:d
                dx = (rand()-0.5)*2*ε
                @inbounds ΔS = daction(links,x,t,dir,dx)   #Berechne Wirkungsunterschied
                if rand()<exp(-ΔS)   #Annahme nach gegebener Wahrscheinlichkeit
                    @inbounds links[x,t,dir] += dx
                end
            end
        end
    end
    return links
end;

function metropolis(start::Array{Float64,3})
    links = deepcopy(start)  #Startkonfiguration
    Q_top = zeros(Int64,Int(N_sweeps/n_skip),1)
    Q_cont = zeros(Float64,Int(N_sweeps/n_skip),1)
    Q_squared_weighted = zeros(Float64,Int(N_sweeps/n_skip),1)
    denom = 0
    for l = 1:N_therm
        links = sweep!(links)
    end
    for l = 1:N_sweeps  #Führe Simulation durch
        links = sweep!(links)
        if l%n_skip == 0
            Qc = cont_charge(links) 
            @inbounds Q_top[i] = top_charge(links)
            @inbounds Q_cont[i] = Qc
            @inbounds Q_squared_weighted[i] = Qc^2*exp(-bias[grid_ind(Qc)]-penalty_potential(Qc)) 
            denom += exp(-bias[grid_ind(Qc)]-penalty_potential(Qc))
            i += 1 
        end
    end
    Q_squared_weighted = Q_squared_weighted./denom
    return Q_top,Q_cont,Q_squared_weighted
end;

function taueffint(obs::Array{Float64,1})
    Gammat = autocor(obs, 1:length(obs)-1)
    W = findfirst(x->x<=0, Gammat[:])
    taueff = Array{Float64}(undef,W)
    for l = 1:W
        @inbounds taueff[l] = 1/log(abs(Gammat[l]/Gammat[l+1]))
    end
    tauint = 1/2
    for l=1:W
        @inbounds tauint += Gammat[l]
    end
    tauint = round(Int,tauint)
    return tauint, taueff, W
end;
