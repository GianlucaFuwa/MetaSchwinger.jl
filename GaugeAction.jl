function plaquette(links,s,t)
    N_s,N_t,dir = size(links)
    p = links[s,t,1]+links[mod1(s+1,N_s),t,2]-links[s,mod1(t+1,N_t),1]-links[s,t,2]
    return exp(p*im)
end

function action_full(configs,beta)
    N_s,N_t,d,N_sweeps = size(configs)
    s = zeros(N_sweeps)
    for i = 1:N_sweeps
        for x = 1:N_s
            for t = 1:N_t
                s[i] += real(beta*(1-plaquette(configs[:,:,:,i],x,t)))
            end
        end
    end
    return s
end;

function action(links,beta)
    N_s,N_t,d = size(links)
    s = 0
    for x = 1:N_s
        for t = 1:N_t
            s += real(beta*(1-plaquette(links,x,t)))
        end
    end
    return s
end;

function daction(links,n,dir,beta,dx)
    N_s,N_t,d = size(links)
    if dir==1
        n_p_mu = mod1.(n.+[1,0],N_s)
        n_p_nu = mod1.(n.+[0,1],N_s)
        n_n_mu = mod1.(n.-[1,0],N_s)
        n_n_nu = mod1.(n.-[0,1],N_s)
        n_nu_mu = mod1.(n.+[1,-1],N_s)
    elseif dir==2
        n_p_nu = mod1.(n.+[1,0],N_s)
        n_p_mu = mod1.(n.+[0,1],N_s)
        n_n_nu = mod1.(n.-[1,0],N_s)
        n_n_mu = mod1.(n.-[0,1],N_s)
        n_nu_mu = mod1.(n.+[-1,1],N_s)
    end
    nu = dir%2+1
    l=links[n[1],n[2],nu]+links[n_p_nu[1],n_p_nu[2],dir]-links[n_p_mu[1],n_p_mu[2],nu]
    r=-links[n_n_nu[1],n_n_nu[2],nu]+links[n_n_nu[1],n_n_nu[2],dir]+links[n_nu_mu[1],n_nu_mu[2],nu]
    A_mu=exp(l*im)+exp(r*im)
    link_prop=links[n[1],n[2],dir]+dx
    ds=-beta*real(  (exp(link_prop*im)-exp(links[n[1],n[2],dir]*im))*conj(A_mu)  )
    return ds
end;

function gauge_trafo(links,trafo)
    links_trafo[:,:,1] = trafo.+links[:,:,1].-circshift(trafo,(-1,0))
    links_trafo[:,:,2] = trafo.+links[:,:,2].-circshift(trafo,(0,-1))
    return links_trafo
end;