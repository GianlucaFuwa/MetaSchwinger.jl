function plaquette(links::Array{Float64,3},s::Int64,t::Int64)
    @inline →(n) = mod1(n+1,N_s) 
    p = links[s,t,1]+links[→(s),t,2]-links[s,→(t),1]-links[s,t,2]
    return exp(p*im)
end

function action_full(configs::Array{Float64,4})
    s = zeros(N_sweeps)
    for i = 1:N_sweeps
        for x = 1:N_s
            for t = 1:N_t
                @inbounds s[i] += real(β*(1-plaquette(configs[:,:,:,i],x,t)))
            end
        end
    end
    return s
end;

function action(links::Array{Float64,3})
    s = 0
    for x = 1:N_s
        for t = 1:N_t
            @inbounds s += real(β*(1-plaquette(links,x,t)))
        end
    end
    return s
end;

function staple(links::Array{Float64,3},s::Int64,t::Int64,dir::Int64)
    @inline →(n) = mod1(n+1,N_s) 
    @inline ←(n) = mod1(n-1,N_s) 
    if dir == 1
        @inbounds l = links[s,t,2]     + links[s,→(t),1] - links[→(s),t,2]
        @inbounds r = -links[s,←(t),2] + links[s,←(t),1] + links[→(s),←(t),2]
    elseif dir == 2
        @inbounds l = links[s,t,1]     + links[→(s),t,2] - links[s,→(t),1]
        @inbounds r = -links[←(s),t,1] + links[←(s),t,2] + links[←(s),→(t),1]
    end
    return exp(l*im)+exp(r*im)
end

function daction(links::Array{Float64,3},s::Int64,t::Int64,dir::Int64,dx::Float64)
    link_old = links[s,t,dir]
    A = staple(links,s,t,dir)
    ds = -β*real(  (exp((link_old+dx)*im)-exp(link_old*im))*conj(A)  )
    return ds
end;

function gauge_trafo(links::Array{Float64,3},trafo::Array{Float64,2})
    links_trafo[:,:,1] = trafo.+links[:,:,1].-circshift(trafo,(-1,0))
    links_trafo[:,:,2] = trafo.+links[:,:,2].-circshift(trafo,(0,-1))
    return links_trafo
end;