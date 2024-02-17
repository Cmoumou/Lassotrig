% IMPORTANT! Please add the Chebfun Toolbox 
% onto path before running this demo

% Lasso trignometric interpolation on the unit circle

% Please note that our plotting are set for accommodating the plotting 
% requirments in our paper. If you choose to test another function, you may
% need to modify plotting setting.
% you can try different values of parameters N, but
% you should remember the trigonoemtric interpolation condintion which requires L= N-1.
clear all;close all
N=501; 
L=N-1;
x=(2*pi)/N:(2*pi)/N:2*pi;
xx=-3*pi:6*pi/2000:3*pi;

 W1=ones(1,N);
w=((2*pi)/N)*W1;

%% function values

 
   example_idx=1
switch example_idx
  case 1
       OO=chebop(0,2*pi)
OO.op=@(x,u) 0.001*diff(u,2)+0.001*diff(u)-cos(x).*u; 
OO.bc='periodic';
oo=OO\0.1;
G_f1=oo(x);
GG_f1=oo(xx);
case 2    
     G_f1 = exp(cos(x)); GG_f1 = exp(cos(xx));
 case 3
       G_f1=cos(50*x).*(1+cos(5*x)/5);GG_f1=cos(50*xx).*(1+cos(5*xx)/5);
    case 4
     f=cheb.gallerytrig('wavepacket');G_f1=f(x);GG_f1=f(xx);
    case 5
        f=cheb.gallerytrig('tsunami');G_f1=f(x);GG_f1=f(xx);
    case 6 
           f=cheb.gallerytrig('random');G_f1=f(x);GG_f1=f(xx);
     case 7 
           f=cheb.gallerytrig('gibbs');G_f1=f(x);GG_f1=f(xx);
      case 8 
          G_f1 = exp(cos(x))+sin(30*x); GG_f1 = exp(cos(xx))+sin(30*xx);         
end




%% add noise



[Y_f1,NOISE_f1]=noisegen(G_f1,10);
[YY_f1,NOISE1_f1]=noisegen(GG_f1,10);







%% producing matrix of least squares problem
for l = 1:L+1
      for j = 1:N
          
    if mod(l, 2) == 0
        A(j,l) =sin(((l)/2)*(x(j)))/sqrt(pi);
  
    else
      A(j,l) = cos(((l-1)/2)*(x(j)))/sqrt(pi);
   
    end
      end
   end
   A(:,1) =1/sqrt(2*pi);
      for l = 1:L+1
      for j = 1:length(xx)
       
    if mod(l, 2) == 0
       C(j,l) =sin(((l)/2)*(xx(j)))/sqrt(pi);
  
    else
      C(j,l) = cos(((l-1)/2)*(xx(j)))/sqrt(pi);
   
    end
      end
    end
  
C(:,1) =1/sqrt(2*pi);
%% parameter choice
 lambda=0.0001:0.0001:0.3;
   for k=1:length(lambda)
   beta_f1 = l1_beta(w,A,Y_f1',lambda(k),L,ones(L+1));
 p_f1=beta_f1*C'; 
  er(k)=norm(p_f1-GG_f1);
   end
 for i=1:L+1
ccc(i)=sin(pi*(i/L))/(pi*(i/L));
end
   [ma,t]=min(er);
   lambdaopt_f1=lambda(t);
   
   %%  coefficients of Lasso trigonometric interpolation 
 
   
   Beta_f1= l1_beta(w,A,Y_f1',0,L,ones(L+1));  
   beta_f1= l1_beta(w,A,Y_f1',lambdaopt_f1,L,ones(L+1));
   
    %%  Lanczos sigma factor
   sigma=5;
   Lanczos_beta_f1=beta_f1*(diag(ccc))^(sigma);
      
 
     P_f1=Beta_f1*C';  
 p_f1=beta_f1*C'; 
 Lanczos_p_f1= Lanczos_beta_f1*C';

 realated_error_P_f1=abs(P_f1-GG_f1);
 realated_error_p_f1=abs(p_f1-GG_f1);
 realated_error_Lanczos_f1=abs(Lanczos_p_f1-GG_f1);
  

%% Figure 1
Color = [215,25,28;
0 0 128;
254,204,92;
102, 0, 204;
255,255,255]/255;

fontsize_baseline = 20;
fontsize_baselinet = 25;
fontsize_baselinea = 15;
figure(1)
axes('position',[0.03 0.58 0.22 0.37]), 
plot(xx,GG_f1,'linewidth',1.2,'color','k'), box on,...
    set(gca, 'fontsize', fontsize_baselinea),...
    xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),...%ylabel('$f(x)$','interpreter','latex', 'fontsize', fontsize_baseline),...
    title('exact function','interpreter','latex', 'fontsize', fontsize_baselinet),...
     grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),axis([-2*pi,2*pi,-10,8.5]),
axes('position',[0.03 0.07 0.22 0.37]), 
plot(xx,YY_f1,'linewidth',1.2,'color','k'),box on,...
    set(gca, 'fontsize', fontsize_baselinea),...
    xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),...%ylabel('$f(x)$','interpreter','latex', 'fontsize', fontsize_baseline),...
    title('noisy function','interpreter','latex', 'fontsize', fontsize_baselinet),...
     grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),axis([-2*pi,2*pi,-10,8.5]),
axes('position',[0.2755 0.58 0.22 0.37]),
plot(xx,P_f1,'linewidth',1.2,'color','k'), box on,...
    set(gca, 'fontsize', fontsize_baselinea),...
    xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),...%ylabel('$f(x)$','interpreter','latex', 'fontsize', fontsize_baseline),...
    title('$\mathcal{T}_n f $ with $N=501$','interpreter','latex', 'fontsize', fontsize_baselinet),  grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),axis([-2*pi,2*pi,-10,8.5]),
axes('position',[0.2755 0.07 0.22 0.37]), 
plot(xx,realated_error_P_f1,'linewidth',1.2,'color','k'),set(gca, 'fontsize', fontsize_baselinea),...
    set(gca, 'fontsize', fontsize_baselinea),...
    xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),...%ylabel('Absolute error', 'interpreter','latex','fontsize', fontsize_baseline),...
    title('error','interpreter','latex','fontsize', fontsize_baselinet),box on, grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),axis([-2*pi,2*pi,0 2.5]),
axes('position',[0.521 0.58 0.22 0.37]),
plot(xx,p_f1,'linewidth',1.2,'color','k'),box on,set(gca, 'fontsize', fontsize_baselinea), ...
    set(gca, 'fontsize', fontsize_baselinea),...
     xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),...%ylabel('$f(x)$','interpreter','latex', 'fontsize', fontsize_baseline),...
     title('$\mathcal{T}_n^{\lambda} f $ with $N=501$','interpreter','latex', 'fontsize', fontsize_baselinet),...
     grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),axis([-2*pi,2*pi,-10,8.5]),
axes('position',[0.521 0.07 0.22 0.37]), set(gca, 'fontsize', fontsize_baselinea), 
plot(xx,realated_error_p_f1,'linewidth',1.2,'color','k'),set(gca, 'fontsize', fontsize_baselinea),...
    set(gca, 'fontsize', fontsize_baselinea),...
    xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),...%ylabel('Absolute error','interpreter','latex', 'fontsize', fontsize_baseline),...
    title('error','interpreter','latex', 'fontsize', fontsize_baselinet),box on,axis([-2*pi,2*pi,0,2.5]), grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off')
 axes('position',[0.7675 0.58 0.22 0.37]), set(gca, 'fontsize', fontsize_baselinea), 
plot(xx, Lanczos_p_f1,'linewidth',1.2,'color','k'),set(gca, 'fontsize', fontsize_baselinea),...
    set(gca, 'fontsize', fontsize_baselinea),...
    xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),...%ylabel('Absolute error','interpreter','latex', 'fontsize', fontsize_baseline),...
    title('$\mathcal{T}_n^{\lambda\sigma} f $ with $N=501$','interpreter','latex', 'fontsize', fontsize_baselinet),box on,axis([-2*pi,2*pi,-10,8.5]), grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off')
axes('position',[0.7675 0.07 0.22 0.37]), set(gca, 'fontsize', fontsize_baselinea), 
plot(xx, realated_error_Lanczos_f1,'linewidth',1.2,'color','k'),set(gca, 'fontsize', fontsize_baselinea),...
    set(gca, 'fontsize', fontsize_baselinea),...
    xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),...%ylabel('Absolute error','interpreter','latex', 'fontsize', fontsize_baseline),...
    title('error','interpreter','latex', 'fontsize', fontsize_baselinet),box on,axis([-2*pi,2*pi,0,2.5]), grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off')