function ret_val = p(y,X,pos_vec)
cs = 0;
for l = 1:6
    pi_l = pos_vec(:,l);
    nrm = vecnorm(X-pi_l',2,2);     % Norm of each row
    cs = cs + (y(l) - 90 + 10*3*log10(nrm)).^2; 
end
ret_val = 1/(2*pi*1.5^2)^3*exp(-1/(2*1.5^2)*cs);
end
