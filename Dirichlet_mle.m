function alpha = Dirichlet_mle(r, k)
% DIRICHLET_MLE  Maximum likelihood esimation of the Dirichelt parameter

alpha = ones(k, 1);
ave_r = mean(mean((log(r)), 3))';

y = psi(sum(alpha)) + ave_r;
for ii = 1:1000
    alpha1 = alpha;
    
    % Initialization
    alpha = exp(y)+1/2;
    i = find(y <= -2.22);
    alpha(i) = -1./(y(i) - psi(1));
    for jj = 1:100
        alpha2 = alpha;
        alpha = alpha - (psi(alpha) - y)./psi(1, alpha);
        if max(abs(alpha - alpha2)) < 1e-6
            break;
        end
    end
    
    if max(abs(alpha1 - alpha)) < 1e-6
        break;
    end
    
    y = psi(sum(alpha)) + ave_r;
end


