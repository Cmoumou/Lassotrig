function beta = l2_beta(w,A,f,lambda,L,mu)
 Af = A'*diag(w)*f;
for i = 1:L+1
       beta(i) = Af(i)/(1+lambda*mu(i)^2);
end
end
