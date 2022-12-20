### Functions for determination of continuous charge and its difference for changes in 1 link-value ###

function cont_charge(links)
    N_s,N_t,dir = size(links)
    charge = 0 
    for s=1:N_s
        for t=1:N_t
            charge+=plaquette(links,s,t)
        end
    end
    return imag(charge)/(2*pi)
end;

function cont_charge_full(configs)
    N_s,N_t,dir,N_sweeps = size(configs)
    Q = complex(zeros(N_sweeps))
    for i=1:N_sweeps    
        for s=1:N_s
            for t=1:N_t
                Q[i]+=plaquette(configs[:,:,:,i],s,t)
            end
        end
    end
    return imag.(Q)./(2*pi)
end;

function dcont_charge(links,n,dir,dx)
    
    N_s,N_t,d = size(links)
    ns = n[1]
    nt = n[2]
    NT = nt%N_t +1  
    NS = ns%N_s +1  
    TN = (nt + N_t -2)%N_t +1   
    SN = (ns + N_s -2)%N_s +1  

    if dir == 1
        a = links[ns,nt,1]+dx + links[NS,nt,2] - links[ns,NT,1] - links[ns,nt,2]
        b = links[ns,TN,1] + links[NS,TN,2] - links[ns,nt,1]-dx - links[ns,TN,2]
        c = links[ns,nt,1] + links[NS,nt,2] - links[ns,NT,1] - links[ns,nt,2]
        d = links[ns,TN,1] + links[NS,TN,2] - links[ns,nt,1] - links[ns,TN,2]
    elseif dir == 2
        a = links[SN,nt,1] + links[ns,nt,2]+dx - links[SN,NT,1] - links[SN,nt,2]
        b = links[ns,nt,1] + links[NS,nt,2] - links[ns,NT,1] - links[ns,nt,2]-dx
        c = links[SN,nt,1] + links[ns,nt,2] - links[SN,NT,1] - links[SN,nt,2]
        d = links[ns,nt,1] + links[NS,nt,2] - links[ns,NT,1] - links[ns,nt,2]
    end

    (sin(a) + sin(b) - sin(c) - sin(d)) / (2*pi)
end
