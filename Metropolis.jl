### Metropolis-Algorithm involving Metadynamics ###

function sweep_meta!(links,beta,epsilon,bias_potential,Q_max,dq,w)   #Metropolis-Sweep Funktion
    N_s,N_t,d = size(links)
    Q = cont_charge(links)
    for s = 1:N_s
        for t = 1:N_t
            for dir=1:d
                dU = (rand()-0.5)*2*epsilon
                Q_prop = Q + dcont_charge(links,[s,t],dir,dU)
                old_ind = grid_ind(Q,dq,Q_max)
                prop_ind = grid_ind(Q_prop,dq,Q_max)
                DeltaS = daction(links,[s,t],dir,beta,dU)  #Berechne Wirkungsunterschied
                DeltaV = bias_potential[prop_ind]-bias_potential[old_ind]
                if rand()<exp(-DeltaS-DeltaV)   #Annahme nach gegebener Wahrscheinlichkeit
                    links[s,t,dir] += dU
                    Q = Q_prop
                    update_bias!(bias_potential,Q,Q_max,dq,w)
                end
            end
        end
    end
    return links
end;

function metropolis_meta(start,N_s,N_t,beta,epsilon,Q_max,dq,w,N_therm,N_sweeps,bias_potential)
    links = deepcopy(start)
    configs = Array{Float64}(undef,N_s,N_t,2,N_sweeps)   #Leeres Array mit Werten für Position jeder Komponente nach jedem Update
    for l = 1:N_therm
        links = sweep!(links,beta,epsilon)
    end
    for l = 1:N_sweeps  #Führe Simulation durch
        links = sweep_meta!(links,beta,epsilon,bias_potential,Q_max,dq,w)
        configs[:,:,:,l] = links
    end
    return configs
end;

### Metropolis-Algorithm without Metadynamics ###

function sweep!(links,beta,epsilon)   #Metropolis-Sweep Funktion
    N_s,N_t,d = size(links)
    for x = 1:N_s
        for t = 1:N_t
            for dir=1:d
                dx = (rand()-0.5)*2*epsilon
                DeltaS = daction(links,[x,t],dir,beta,dx)   #Berechne Wirkungsunterschied
                if rand()<exp(-DeltaS)   #Annahme nach gegebener Wahrscheinlichkeit
                    links[x,t,dir] += dx
                end
            end
        end
    end
    return links
end;

function metropolis(N_sweeps,N_therm,N_s,N_t,m,epsilon)
    links = zeros(N_s,N_t,2)  #Startkonfiguration
    configs = Array{Float64}(undef,N_s,N_t,2,N_sweeps)   #Leeres Array mit Werten für Position jeder Komponente nach jedem Update
    for l = 1:N_therm
        links = sweep!(links,m,epsilon)
    end
    for l = 1:N_sweeps  #Führe Simulation durch
        links = sweep!(links,beta,epsilon)
        configs[:,:,:,l] = links
    end
    return configs
end;

function taueffint(configs, length::Int, obs=nothing)
    if obs == nothing
        means = mean(configs,dims=[1,2,3])
        Gammat = autocor(means[:], 1:length)
    else
        Gammat = autocor(obs(configs), 1:length)
    end
    W = findfirst(x->x<=0, Gammat[:])
    taueff = Array{Float64}(undef,W)
    for l = 1:W
        global taueff[l] = 1/log(abs(Gammat[l]/Gammat[l+1]))
    end
    tauint = 1/2
    for l=1:W
        tauint += Gammat[l]
    end
    tauint = round(Int,tauint)
    return tauint, taueff, W
end;
