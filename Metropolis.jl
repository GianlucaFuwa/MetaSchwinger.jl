### Metropolis-Algorithm involving Metadynamics ###

function sweep_meta!(links::Array{Float64,3},bias::Array{Float64,1},is_static::Bool=false)   #Metropolis-Sweep Funktion
    Q = cont_charge(links)
    for eo in [0,1]
        for s = (eo+1):2:N_s
            for t = 1:N_t
                dU = (rand()-0.5)*2*ε
                old_ind = grid_ind(Q)
                Q_prop = Q + dcont_charge(links,s,t,2,dU)
                prop_ind = grid_ind(Q_prop)
                ΔS = daction(links,s,t,2,dU)  #Berechne Wirkungsunterschied
                ΔV = bias[prop_ind]-bias[old_ind]
                ΔV_pen = penalty_potential(Q_prop)-penalty_potential(Q)
                if rand()<exp(-ΔS-ΔV-ΔV_pen)   #Annahme nach gegebener Wahrscheinlichkeit
                    links[s,t,2] += dU
                    Q = Q_prop
                    global acc += 1
                    if is_static == false
                        update_bias!(bias,Q)
                    end
                end
            end
        end
    end
    for eo in [0,1]
        for t = (eo+1):2:N_t
            for s = 1:N_s
                dU = (rand()-0.5)*2*ε
                old_ind = grid_ind(Q)
                Q_prop = Q + dcont_charge(links,s,t,1,dU)
                prop_ind = grid_ind(Q_prop)
                ΔS = daction(links,s,t,1,dU)  #Berechne Wirkungsunterschied
                ΔV = bias[prop_ind]-bias[old_ind]
                ΔV_pen = penalty_potential(Q_prop)-penalty_potential(Q)
                if rand()<exp(-ΔS-ΔV-ΔV_pen)   #Annahme nach gegebener Wahrscheinlichkeit
                    links[s,t,1] += dU
                    Q = Q_prop
                    global acc += 1
                    if is_static == false
                        update_bias!(bias,Q)
                    end
                end
            end
        end
    end
    return links,bias
end;

function metropolis_meta(start::Array{Float64,3},bias::Array{Float64,1},is_static::Bool=false)
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
        links,bias = sweep_meta!(links,bias,is_static)
        if l%n_skip == 0
            Qc = cont_charge(links) 
            @inbounds Q_top[i] = top_charge(links)
            @inbounds Q_cont[i] = Qc
            @inbounds Q_squared_weighted[i] = Qc^2*exp(bias[grid_ind(Qc)]+penalty_potential(Qc)) 
            denom += exp(bias[grid_ind(Qc)]+penalty_potential(Qc))
            i += 1 
        end
    end
    display(denom)
    Q_squared_weighted = Q_squared_weighted./denom
    println("acc. rate = $(100*acc/(N_sweeps*N_s*N_t*d))%")
    return Q_top,Q_cont,Q_squared_weighted,bias
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
