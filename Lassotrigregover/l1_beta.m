function beta = l1_beta(w,A,f,lambda,L,mu)
% For the coefficients of $\ell_1-$regularized approximation
Af = A'*diag(w)*f;
ALA = A'*diag(w)*A;
for i = 0:L
    if 2*Af(i+1)>lambda*mu(i+1)
        beta(i+1) = (2*Af(i+1)-lambda*mu(i+1))/(2*ALA(i+1,i+1));
    elseif 2*Af(i+1)<-lambda*mu(i+1)
        beta(i+1) = (2*Af(i+1)+lambda*mu(i+1))/(2*ALA(i+1,i+1));
    else
        beta(i+1) = 0;
    end
end
end