### Function to find point on potential grid for given Q ###
function grid_ind(q::Float64)
    grid_index = (q+Q_max)/dq+0.5
    return round(Int,grid_index)
end;

### Function to update bias potentia for given Q ###
function update_bias!(bias_potential::Array{Float64,1},q::Float64)
    grid_index = grid_ind(q)
    grid_q = -Q_max+(grid_index-1)*dq

    if grid_index >= 1 && grid_index < length(bias_potential)
        bias_potential[grid_index] += w*(1-(q-grid_q)/dq)
        bias_potential[grid_index+1] += w*(q-grid_q)/dq
    else
        boundary_count += 1
    end
end;

function penalty_potential(q::Float64)
    if abs(q) >= Q_thresh
        p_pot = k*min((q+Q_thresh)^2,(q-Q_thresh)^2)
    else 
        p_pot = 0
    end
    return p_pot
end;


