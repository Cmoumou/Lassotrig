clear all;close all
N=501; 
L=500;
K=randn(1,2);
x=(2*pi)/N:(2*pi)/N:2*pi;
xx=-5*pi:10*pi/2000:5*pi;
 W1=ones(1,N);
w=((2*pi)/N)*W1;
%% 五个函数赋值
 G_f5=tanh(50*sin(x));
%G_f5=x.*cos(x).*sin(50*x);
%G_f5=(2/pi)*atan(tan(x/2));
GG_f5=tanh(50*sin(xx));
%GG_f5=xx.*cos(xx).*sin(50*xx);
%GG_f5=(2/pi)*atan(tan(xx/2));



%% 函数加噪



[Y_f5,NOISE_f5]=noisegen(G_f5,30);
[YY_f5,NOISE1_f5]=noisegen(GG_f5,30);







%% 最小二乘矩阵生成
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
%% 参数选取
 lambda=0.0001:0.0001:0.3;
   for k=1:length(lambda)
   beta_f5 = l1_beta(w,A,Y_f5',lambda(k),L,ones(L+1));
 p_f5=beta_f5*C'; 
  er(k)=norm(p_f5-GG_f5);
   end
     for i=1:L+1
ccc(i)=sin(pi*(i/L))/(pi*(i/L));
end
   [ma,t]=min(er);
   lambdaopt_f5=lambda(t);
   
   %%  三角多项式系数确定
 
   
        Beta_f5= l1_beta(w,A,Y_f5',0,L,ones(L+1));  
   beta_f5= l1_beta(w,A,Y_f5',lambdaopt_f5,L,ones(L+1));
   Lanczos_beta_f5=beta_f5*(diag(ccc))^(30);
      
 
     P_f5=Beta_f5*C';  
 p_f5=beta_f5*C'; 
 Lanczos_p_f5= Lanczos_beta_f5*C';

 realated_error_P_f5=abs(P_f5-GG_f5);
 realated_error_p_f5=abs(p_f5-GG_f5);
 realated_error_Lanczos_f5=abs(Lanczos_p_f5-GG_f5);
  

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
plot(xx,GG_f5,'linewidth',1.2,'color','k'), box on,...
    set(gca, 'fontsize', fontsize_baselinea),...
    xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),...%ylabel('$f(x)$','interpreter','latex', 'fontsize', fontsize_baseline),...
    title('Exact function','interpreter','latex', 'fontsize', fontsize_baselinet),...
     grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),axis([-5*pi,5*pi,-1.2,1.2]),
axes('position',[0.03 0.07 0.22 0.37]), 
plot(xx,YY_f5,'linewidth',1.2,'color','k'),box on,...
    set(gca, 'fontsize', fontsize_baselinea),...
    xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),...%ylabel('$f(x)$','interpreter','latex', 'fontsize', fontsize_baseline),...
    title('Noisy function','interpreter','latex', 'fontsize', fontsize_baselinet),...
     grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),axis([-5*pi,5*pi,-1.2,1.2]),
axes('position',[0.2755 0.58 0.22 0.37]),
plot(xx,P_f5,'linewidth',1.2,'color','k'), box on,...
    set(gca, 'fontsize', fontsize_baselinea),...
    xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),...%ylabel('$f(x)$','interpreter','latex', 'fontsize', fontsize_baseline),...
    title('$\mathcal{T}_n f $ with $N=501$','interpreter','latex', 'fontsize', fontsize_baselinet),  grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),axis([-5*pi,5*pi,-1.2,1.2]),
axes('position',[0.2755 0.07 0.22 0.37]), 
plot(xx,realated_error_P_f5,'linewidth',1.2,'color','k'),set(gca, 'fontsize', fontsize_baselinea),...
    set(gca, 'fontsize', fontsize_baselinea),...
    xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),...%ylabel('Absolute error', 'interpreter','latex','fontsize', fontsize_baseline),...
    title('Error','interpreter','latex','fontsize', fontsize_baselinet),box on, grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),axis([-5*pi,5*pi,0,0.5]),
axes('position',[0.521 0.58 0.22 0.37]),
plot(xx,p_f5,'linewidth',1.2,'color','k'),box on,set(gca, 'fontsize', fontsize_baselinea), ...
    set(gca, 'fontsize', fontsize_baselinea),...
     xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),...%ylabel('$f(x)$','interpreter','latex', 'fontsize', fontsize_baseline),...
     title('$\mathcal{T}_n^{\lambda} f $ with $N=501$','interpreter','latex', 'fontsize', fontsize_baselinet),...
     grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),axis([-5*pi,5*pi,-1.2,1.2]),
axes('position',[0.521 0.07 0.22 0.37]), set(gca, 'fontsize', fontsize_baselinea), 
plot(xx,realated_error_p_f5,'linewidth',1.2,'color','k'),set(gca, 'fontsize', fontsize_baselinea),...
    set(gca, 'fontsize', fontsize_baselinea),...
    xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),...%ylabel('Absolute error','interpreter','latex', 'fontsize', fontsize_baseline),...
    title('Error','interpreter','latex', 'fontsize', fontsize_baselinet),box on,axis([-5*pi,5*pi,0,0.5]), grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off')
 axes('position',[0.7675 0.58 0.22 0.37]), set(gca, 'fontsize', fontsize_baselinea), 
plot(xx, Lanczos_p_f5,'linewidth',1.2,'color','k'),set(gca, 'fontsize', fontsize_baselinea),...
    set(gca, 'fontsize', fontsize_baselinea),...
    xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),...%ylabel('Absolute error','interpreter','latex', 'fontsize', fontsize_baseline),...
    title('$\mathcal{T}_n^{\lambda\sigma} f $ with $N=501$','interpreter','latex', 'fontsize', fontsize_baselinet),box on,axis([-5*pi,5*pi,-1.2,1.2]), grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off')
axes('position',[0.7675 0.07 0.22 0.37]), set(gca, 'fontsize', fontsize_baselinea), 
plot(xx, realated_error_Lanczos_f5,'linewidth',1.2,'color','k'),set(gca, 'fontsize', fontsize_baselinea),...
    set(gca, 'fontsize', fontsize_baselinea),...
    xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),...%ylabel('Absolute error','interpreter','latex', 'fontsize', fontsize_baseline),...
    title('Error','interpreter','latex', 'fontsize', fontsize_baselinet),box on,axis([-5*pi,5*pi,0,0.5]), grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off')