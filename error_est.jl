function taueffint(obs)
    Gammat = autocor(obs, 1:length(obs)-1)
    W = findfirst(x->x<=0, Gammat[:])
    taueff = Array{Float64,1}(undef,W)
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

function bootstrap(obs,N_boot::Int64,weighted::Bool=false)
    N = round(Int,N_sweeps/n_skip)
    B, = taueffint(obs)
    NdivB = round(Int,N/B)
    A = zeros(N_boot,1)
    B = 2*round(Int,B)+1
    for i = 1:N_boot
        boot_obs = []
        for l = 1:NdivB
            r = rand(1:N)
            q = range(r,r+B-1)
            boot_obs = append!(boot_obs,obs[mod1.(q,N)])
        end
        if weighted==false
            A[i] = mean(boot_obs)
        elseif weighted==true
            A[i] = sum(boot_obs)
        end
    end
    obs_mean = mean(A)
    obs_std  = std(A)*sqrt(N_boot-1)
    return obs_mean,obs_std
end;