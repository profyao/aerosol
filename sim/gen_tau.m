function tau = gen_tau(kappa,XDim_r,YDim_r)

Q = igmrfprec([XDim_r, YDim_r], 1);
R = chol(Q+eye(XDim_r*YDim_r)*1e-3);

z = 1/sqrt(kappa) * randn(XDim_r*YDim_r,1) ;
tau = R\z + 0.3;

%imagesc(reshape(tau,XDim_r,YDim_r)),colorbar

end

