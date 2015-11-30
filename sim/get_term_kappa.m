function term_kappa = get_term_kappa(i,j,p,tau,taup)

    id = find(j == p);
    if ~isempty(id)
        neighbor = i(id);
        tau_neighbor = tau(neighbor);
    else
        tau_neighbor = [];
    end
    
    if isempty(tau_neighbor)
        
        term_kappa = 0;
    else
        term_kappa = sum((taup - tau_neighbor).^2); 
    end

end