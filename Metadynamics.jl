function grid_ind(q,dq,Q_max)
    grid_index = (q+Q_max)/dq+0.5
    return round(Int,grid_index)
end;

function update_bias!(bias_potential,q,Q_max,dq,w)
    grid_index = grid_ind(q,dq,Q_max)
    grid_q = -Q_max+(grid_index-1)*dq

    bias_potential[grid_index] += w*(1-(q-grid_q)/dq)
    #if grid_index <= grid_ind(Q_max,dq) && grid_index >= grid_ind(-Q_max,dq)
    #    bias_potential[grid_index+1] += w+(q-grid_q)/dq
    #elseif grid_index >= grid_ind(Q_max,dq) || grid_index <= grid_ind(-Q_max,dq)
    #    boundary_count += 1
    #end
end;