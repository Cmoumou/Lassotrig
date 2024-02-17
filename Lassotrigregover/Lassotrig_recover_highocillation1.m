clear all;close all
N=501; 
L=500;
x=(2*pi)/N:(2*pi)/N:2*pi;
xx=-9:0.0001:9;

 W1=ones(1,N);
w=((2*pi)/N)*W1;
OO=chebop(0,2*pi)
OO.op=@(x,u) 0.001*diff(u,2)+0.001*diff(u)-cos(x).*u; 
OO.bc='periodic';
oo=OO\0.1;
%% 五个函数赋值

G_f4=oo(x);
GG_f4=oo(xx);





%% 函数加噪



[Y_f4,NOISE_f4]=noisegen(G_f4,10);
[YY_f4,NOISE1_f4]=noisegen(GG_f4,10);







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
   beta_f4 = l1_beta(w,A,Y_f4',lambda(k),L,ones(L+1));
 p_f4=beta_f4*C'; 
  er(k)=norm(p_f4-GG_f4);
   end
     for i=1:L+1
ccc(i)=sin(pi*(i/L))/(pi*(i/L));
end
   [ma,t]=min(er);
   lambdaopt_f4=lambda(t);
   
   %%  三角多项式系数确定
 
   
        Beta_f4= l1_beta(w,A,Y_f4',0,L,ones(L+1));  
   beta_f4= l1_beta(w,A,Y_f4',lambdaopt_f4,L,ones(L+1));
   Lanczos_beta_f4=beta_f4*(diag(ccc))^(5);
      
 
     P_f4=Beta_f4*C';  
 p_f4=beta_f4*C'; 
 Lanczos_p_f4= Lanczos_beta_f4*C';

 realated_error_P_f4=abs(P_f4-GG_f4);
 realated_error_p_f4=abs(p_f4-GG_f4);
 realated_error_Lanczos_f4=abs(Lanczos_p_f4-GG_f4);
  

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
plot(xx,GG_f4,'linewidth',1.2,'color','k'), box on,...
    set(gca, 'fontsize', fontsize_baselinea),...
    xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),...%ylabel('$f(x)$','interpreter','latex', 'fontsize', fontsize_baseline),...
    title('Exact function','interpreter','latex', 'fontsize', fontsize_baselinet),...
     grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),axis([-2*pi,2*pi,-10,8.5]),
axes('position',[0.03 0.07 0.22 0.37]), 
plot(xx,YY_f4,'linewidth',1.2,'color','k'),box on,...
    set(gca, 'fontsize', fontsize_baselinea),...
    xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),...%ylabel('$f(x)$','interpreter','latex', 'fontsize', fontsize_baseline),...
    title('Noisy function','interpreter','latex', 'fontsize', fontsize_baselinet),...
     grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),axis([-2*pi,2*pi,-10,8.5]),
axes('position',[0.2755 0.58 0.22 0.37]),
plot(xx,P_f4,'linewidth',1.2,'color','k'), box on,...
    set(gca, 'fontsize', fontsize_baselinea),...
    xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),...%ylabel('$f(x)$','interpreter','latex', 'fontsize', fontsize_baseline),...
    title('$\mathcal{T}_n f $ with $N=501$','interpreter','latex', 'fontsize', fontsize_baselinet),  grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),axis([-2*pi,2*pi,-10,8.5]),
axes('position',[0.2755 0.07 0.22 0.37]), 
plot(xx,realated_error_P_f4,'linewidth',1.2,'color','k'),set(gca, 'fontsize', fontsize_baselinea),...
    set(gca, 'fontsize', fontsize_baselinea),...
    xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),...%ylabel('Absolute error', 'interpreter','latex','fontsize', fontsize_baseline),...
    title('Error','interpreter','latex','fontsize', fontsize_baselinet),box on, grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),axis([-2*pi,2*pi,0,2.5]),
axes('position',[0.521 0.58 0.22 0.37]),
plot(xx,p_f4,'linewidth',1.2,'color','k'),box on,set(gca, 'fontsize', fontsize_baselinea), ...
    set(gca, 'fontsize', fontsize_baselinea),...
     xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),...%ylabel('$f(x)$','interpreter','latex', 'fontsize', fontsize_baseline),...
     title('$\mathcal{T}_n^{\lambda} f $ with $N=501$','interpreter','latex', 'fontsize', fontsize_baselinet),...
     grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),axis([-2*pi,2*pi,-10,8.5]),
axes('position',[0.521 0.07 0.22 0.37]), set(gca, 'fontsize', fontsize_baselinea), 
plot(xx,realated_error_p_f4,'linewidth',1.2,'color','k'),set(gca, 'fontsize', fontsize_baselinea),...
    set(gca, 'fontsize', fontsize_baselinea),...
    xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),...%ylabel('Absolute error','interpreter','latex', 'fontsize', fontsize_baseline),...
    title('Error','interpreter','latex', 'fontsize', fontsize_baselinet),box on,axis([-2*pi,2*pi,0,2.5]), grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off')
 axes('position',[0.7675 0.58 0.22 0.37]), set(gca, 'fontsize', fontsize_baselinea), 
plot(xx, Lanczos_p_f4,'linewidth',1.2,'color','k'),set(gca, 'fontsize', fontsize_baselinea),...
    set(gca, 'fontsize', fontsize_baselinea),...
    xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),...%ylabel('Absolute error','interpreter','latex', 'fontsize', fontsize_baseline),...
    title('$\mathcal{T}_n^{\lambda\sigma} f $ with $N=501$','interpreter','latex', 'fontsize', fontsize_baselinet),box on,axis([-2*pi,2*pi,-10,8.5]), grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off')
axes('position',[0.7675 0.07 0.22 0.37]), set(gca, 'fontsize', fontsize_baselinea), 
plot(xx, realated_error_Lanczos_f4,'linewidth',1.2,'color','k'),set(gca, 'fontsize', fontsize_baselinea),...
    set(gca, 'fontsize', fontsize_baselinea),...
    xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),...%ylabel('Absolute error','interpreter','latex', 'fontsize', fontsize_baseline),...
    title('Error','interpreter','latex', 'fontsize', fontsize_baselinet),box on,axis([-2*pi,2*pi,0,2.5]), grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off')