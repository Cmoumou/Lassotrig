clear all;close all
N=501;
N1=501;
L=500;
xx=-3*pi:6*pi/2000:3*pi;
K=randn(1,2);



x=(2*pi)/N:(2*pi)/N:2*pi;
x1=(2*pi)/N1:(2*pi)/N1:2*pi;
 W1=ones(1,N);
w=((2*pi)/N)*W1;
   %% f_1(x) add noise      
G_f1= exp(cos(x));
G1_f1= exp(cos(x1));
GG_f1=  exp(cos(xx));
[Y_f1,NOISE_f1]=noisegen(G_f1,10);
[YY_f1,NOISE1_f1]=noisegen(G1_f1,10);
[yy_f1,NOISE2_f1]=noisegen(GG_f1,10);
delta_f1=norm(NOISE2_f1);

%% f_2(x) add noise
G_f2= K(1)*cos(50*x).*(1+cos(5*x)/5)+K(2)*cos(50*x+4*sin(5*x));
G1_f2= K(1)*cos(50*x1).*(1+cos(5*x1)/5)+K(2)*cos(50*x1+4*sin(5*x1));
GG_f2= K(1)*cos(50*xx).*(1+cos(5*xx)/5)+K(2)*cos(50*xx+4*sin(5*xx));
[Y_f2,NOISE_f2]=noisegen(G_f2,10);
[YY_f2,NOISE1_f2]=noisegen(G1_f2,10);
[yy_f2,NOISE2_f2]=noisegen(GG_f2,10);
delta_f2=norm(NOISE2_f2);

 %% f_3(x) with h=3 add noise 
G_f3=(2*(3)/pi)*asin(sin(x)); 
GG_f3=(2*(3)/pi)*asin(sin(xx));
G1_f3=(2*(3)/pi)*asin(sin(x1));
[Y_f3,NOISE_f3]=noisegen(G_f3,10);
[YY_f3,NOISE1_f3]=noisegen(G1_f3,10);
[yy_f3,NOISE2_f3]=noisegen(GG_f3,10);
delta_f3=norm(NOISE2_f3);
%% matrix A 
 for l = 1:L+1
      for j = 1:N
          
    if mod(l, 2) == 0
        A1(j,l) =sin(((l)/2)*(x(j)))/sqrt(pi);
  
    else
      A1(j,l) = cos(((l-1)/2)*(x(j)))/sqrt(pi);
   
    end
      end
   end
   A1(:,1) =1/sqrt(2*pi);
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
%% 参数选取
 lambda=0.07:0.0001:0.3;
   for k=1:length(lambda)
   beta_f1 = l1_beta(w,A1,Y_f1',lambda(k),L,ones(L+1));
   beta_f2 = l1_beta(w,A1,Y_f2',lambda(k),L,ones(L+1));
   beta_f3 = l1_beta(w,A1,Y_f3',lambda(k),L,ones(L+1));
   p_f1=beta_f1*C'; p_f2=beta_f2*C'; p_f3=beta_f3*C'; 
   er(k,1)=norm(p_f1-GG_f1);er(k,2)=norm(p_f2-GG_f2);er(k,3)=norm(p_f3-GG_f3);
   end
     for i=1:L+1
ccc(i)=sin(pi*(i/L))/(pi*(i/L));
end
   [ma,t]=min(er);
   lambda_f1=lambda(t(1));lambda_f2=lambda(t(2));lambda_f3=lambda(t(3));



%% determin regularization solution
beta_f1 = l2_beta(w,A1,Y_f1',lambda_f1,L,ones(L+1));
beta_f1=beta_f1*(diag(ccc))^(0);
Beta_f1 = l2_beta(w,A1,Y_f1',0,L,ones(L+1));
beta_f2 = l2_beta(w,A1,Y_f2',lambda_f2,L,ones(L+1));
Beta_f2 = l2_beta(w,A1,Y_f2',0,L,ones(L+1));
beta_f3 = l2_beta(w,A1,Y_f3',lambda_f3,L,ones(L+1));
beta_f3=beta_f3*(diag(ccc))^(0);
Beta_f3 = l2_beta(w,A1,Y_f3',0,L,ones(L+1));

%% approximation trigonometric polynomial
 p_f1=beta_f1*C';
 P_f1=Beta_f1*C';
 p_f2=beta_f2*C';
 P_f2=Beta_f2*C';
 p_f3=beta_f3*C';
 P_f3=Beta_f3*C';
%% plot

fontsize_baselinet=12;
fontsize_baseline=10;
Color = [215,25,28;
0 0 128;
254,204,92;
171,217,233;
44,123,182]/255;
%图第一行
          axes('position',[0.015 0.71 0.153 0.262]),
 plot(xx,GG_f1,'color',Color(1,:),'linewidth',1.2); set(gca, 'fontsize', fontsize_baseline),set(gca, 'XMinorGrid', 'off'),set(gca, 'YMinorGrid', 'off'), grid on; 
 title('Exact $f_1(x)$' ,'interpreter','latex', 'fontsize', fontsize_baselinet)
xlabel('$x$ ','interpreter','latex', 'fontsize', fontsize_baselinet),...
% ylabel('$f(x)$', 'interpreter','latex','fontsize', fontsize_baselinet),...
 axis([-3*pi,3*pi, -1,4])
          axes('position',[0.181 0.71 0.153 0.262]),
plot(xx,yy_f1,'color',Color(1,:),'linewidth',1.2); set(gca, 'fontsize', fontsize_baseline),set(gca, 'XMinorGrid', 'off'),set(gca, 'YMinorGrid', 'off'), grid on; 
title('$f_1(x)$ with noise' ,'interpreter','latex', 'fontsize', fontsize_baselinet)
xlabel('$x$ ','interpreter','latex', 'fontsize', fontsize_baselinet),...
% ylabel('$f(x)$', 'interpreter','latex','fontsize', fontsize_baselinet),...
axis([-3*pi,3*pi, -1,4])
          axes('position',[0.347 0.71 0.153 0.262]),
plot(xx,GG_f2,'color',Color(2,:),'linewidth',1.2); set(gca, 'fontsize', fontsize_baseline),set(gca, 'XMinorGrid', 'off'),set(gca, 'YMinorGrid', 'off'), grid on; 
 title('Exact $f_2(x)$' ,'interpreter','latex', 'fontsize', fontsize_baselinet)
xlabel('$x$ ','interpreter','latex', 'fontsize', fontsize_baselinet),...
% ylabel('$f(x)$', 'interpreter','latex','fontsize', fontsize_baselinet),...
axis([-2,2, -1.5,1.5])
         axes('position',[0.513 0.71 0.153 0.262]),
 plot(xx,yy_f2,'color',Color(2,:),'linewidth',1.2); set(gca, 'fontsize', fontsize_baseline),set(gca, 'XMinorGrid', 'off'),set(gca, 'YMinorGrid', 'off'), grid on; 
title('$f_2(x)$ with noise' ,'interpreter','latex', 'fontsize', fontsize_baselinet)
xlabel('$x$ ','interpreter','latex', 'fontsize', fontsize_baselinet),...
% ylabel('$f(x)$', 'interpreter','latex','fontsize', fontsize_baselinet),...
axis([-2,2, -1.5,1.5])
         axes('position',[0.679 0.71 0.153 0.262]),
   plot(xx,GG_f3,'color',Color(5,:),'linewidth',1.2); set(gca, 'fontsize', fontsize_baseline),set(gca, 'XMinorGrid', 'off'),set(gca, 'YMinorGrid', 'off'), grid on; 
 title('Exact $f_3(x)$' ,'interpreter','latex', 'fontsize', fontsize_baselinet)
xlabel('$x$ ','interpreter','latex', 'fontsize', fontsize_baselinet),...
% ylabel('$f(x)$', 'interpreter','latex','fontsize', fontsize_baselinet),...
 axis([-3*pi,3*pi, -4,4])
        axes('position',[0.845 0.71 0.153 0.262]),
plot(xx,yy_f3,'color',Color(5,:),'linewidth',1.2); set(gca, 'fontsize', fontsize_baseline),set(gca, 'XMinorGrid', 'off'),set(gca, 'YMinorGrid', 'off'), grid on; 
title('$f_3(x)$ with noise' ,'interpreter','latex', 'fontsize', fontsize_baselinet)
xlabel('$x$ ','interpreter','latex', 'fontsize', fontsize_baselinet),...
% ylabel('$f(x)$', 'interpreter','latex','fontsize', fontsize_baselinet),...
 axis([-3*pi,3*pi, -4,4])
 
 
 
 
 
 
 
 
 
 %图第二行
          axes('position',[0.015 0.376 0.153 0.262]),
 plot(xx,P_f1,'color',Color(1,:),'linewidth',1.2); set(gca, 'fontsize', fontsize_baseline),set(gca, 'XMinorGrid', 'off'),set(gca, 'YMinorGrid', 'off'), grid on; 
  title('Recover $f_1(x)$ by trigonometric interpolation' ,'interpreter','latex', 'fontsize', fontsize_baselinet)
xlabel('$x$ ','interpreter','latex', 'fontsize', fontsize_baselinet),...
% ylabel('$f(x)$', 'interpreter','latex','fontsize', fontsize_baselinet),...
 axis([-3*pi,3*pi, -1,4])
         axes('position',[0.181 0.376 0.153 0.262]),
  plot(xx,abs(P_f1-GG_f1),'color',Color(1,:),'linewidth',1.2);set(gca, 'fontsize', fontsize_baseline), set(gca, 'XMinorGrid', 'off'),set(gca, 'YMinorGrid', 'off'), grid on; 
title('Error' ,'interpreter','latex', 'fontsize', fontsize_baselinet)
xlabel('$x$ ','interpreter','latex', 'fontsize', fontsize_baselinet),...
% ylabel('Absolute error', 'interpreter','latex','fontsize', fontsize_baselinet),... 
axis([-3*pi,3*pi, 0,1.8])
        axes('position',[0.347  0.376 0.153 0.262]),
plot(xx,P_f2,'color',Color(2,:),'linewidth',1.2); set(gca, 'fontsize', fontsize_baseline), set(gca, 'fontsize', fontsize_baseline),set(gca, 'XMinorGrid', 'off'),set(gca, 'YMinorGrid', 'off'), grid on; 
title('Recover $f_2(x)$ by trigonometric interpolation' ,'interpreter','latex', 'fontsize', fontsize_baselinet)
xlabel('$x$ ','interpreter','latex', 'fontsize', fontsize_baselinet),...
% ylabel('$f(x)$', 'interpreter','latex','fontsize', fontsize_baselinet),...
axis([-2,2, -1.5,1.5])
       axes('position',[0.513  0.376 0.153 0.262]),
plot(xx,abs(P_f2-GG_f2),'color',Color(2,:),'linewidth',1.2); set(gca, 'fontsize', fontsize_baseline),set(gca, 'XMinorGrid', 'off'),set(gca, 'YMinorGrid', 'off'), grid on; 
title('Error' ,'interpreter','latex', 'fontsize', fontsize_baselinet)
xlabel('$x$ ','interpreter','latex', 'fontsize', fontsize_baselinet),...
%ylabel('Absolute error', 'interpreter','latex','fontsize', fontsize_baselinet),... 
axis([-2,2, 0,1])
      axes('position',[0.679  0.376 0.153 0.262]),
  plot(xx,P_f3,'color',Color(5,:),'linewidth',1.2); set(gca, 'fontsize', fontsize_baseline),set(gca, 'XMinorGrid', 'off'),set(gca, 'YMinorGrid', 'off'), grid on; 
 title('Recover $f_3(x)$ by trigonometric interpolation' ,'interpreter','latex', 'fontsize', fontsize_baselinet)
xlabel('$x$ ','interpreter','latex', 'fontsize', fontsize_baselinet),...
% ylabel('$f(x)$', 'interpreter','latex','fontsize', fontsize_baselinet),...
 axis([-3*pi,3*pi, -4,4])
     axes('position',[0.845  0.376 0.153 0.262]),
 plot(xx,abs(P_f3-GG_f3),'color',Color(5,:),'linewidth',1.2); set(gca, 'fontsize', fontsize_baseline),set(gca, 'XMinorGrid', 'off'),set(gca, 'YMinorGrid', 'off'), grid on; 
title('Error' ,'interpreter','latex', 'fontsize', fontsize_baselinet)
xlabel('$x$ ','interpreter','latex', 'fontsize', fontsize_baselinet),...
% ylabel('Absolute error', 'interpreter','latex','fontsize', fontsize_baselinet),... 
 axis([-3*pi,3*pi, 0,2])


%图第三行
     axes('position',[0.015 0.05 0.153 0.262]),
plot(xx,p_f1,'color',Color(1,:),'linewidth',1.2);  set(gca, 'fontsize', fontsize_baseline),set(gca, 'XMinorGrid', 'off'),set(gca, 'YMinorGrid', 'off'), grid on; 
title('Recover $f_1(x)$ by LTI' ,'interpreter','latex', 'fontsize', fontsize_baselinet)
xlabel('$x$ ','interpreter','latex', 'fontsize', fontsize_baselinet),...
%ylabel('$f(x)$', 'interpreter','latex','fontsize', fontsize_baselinet),...
axis([-3*pi,3*pi, -1,4])
     axes('position',[0.181  0.05 0.153 0.262]),
plot(xx,abs(p_f1-GG_f1),'color',Color(1,:),'linewidth',1.2);set(gca, 'fontsize', fontsize_baseline),  set(gca, 'XMinorGrid', 'off'),set(gca, 'YMinorGrid', 'off'), grid on; 
title('Error' ,'interpreter','latex', 'fontsize', fontsize_baselinet)
xlabel('$x$ ','interpreter','latex', 'fontsize', fontsize_baselinet),...
% ylabel('Absolute error', 'interpreter','latex','fontsize', fontsize_baselinet),... 
axis([-3*pi,3*pi, 0,1.8]),
     axes('position',[0.347 0.05 0.153 0.262]),
plot(xx,p_f2,'color',Color(2,:),'linewidth',1.2);set(gca, 'fontsize', fontsize_baseline), set(gca, 'XMinorGrid', 'off'),set(gca, 'YMinorGrid', 'off'), grid on; 
title('Recover $f_2(x)$ by LTI' ,'interpreter','latex', 'fontsize', fontsize_baselinet)
xlabel('$x$ ','interpreter','latex', 'fontsize', fontsize_baselinet),...
% ylabel('$f(x)$', 'interpreter','latex','fontsize', fontsize_baselinet),...
axis([-2,2, -1.5,1.5])
     axes('position',[0.513 0.05 0.153 0.262]),
plot(xx,abs(p_f2-GG_f2),'color',Color(2,:),'linewidth',1.2); set(gca, 'fontsize', fontsize_baseline), set(gca, 'XMinorGrid', 'off'),set(gca, 'YMinorGrid', 'off'), grid on; 
title('Error' ,'interpreter','latex', 'fontsize', fontsize_baselinet)
xlabel('$x$ ','interpreter','latex', 'fontsize', fontsize_baselinet),...
% ylabel('Absolute error', 'interpreter','latex','fontsize', fontsize_baselinet),... 
axis([-2,2, 0,1])
    axes('position',[0.679 0.05 0.153 0.262]),
plot(xx,p_f3,'color',Color(5,:),'linewidth',1.2); set(gca, 'fontsize', fontsize_baseline),set(gca, 'XMinorGrid', 'off'),set(gca, 'YMinorGrid', 'off'), grid on; 
title('Recover $f_3(x)$ by LTI' ,'interpreter','latex', 'fontsize', fontsize_baselinet)
xlabel('$x$ ','interpreter','latex', 'fontsize', fontsize_baselinet),...
ylabel('$f(x)$', 'interpreter','latex','fontsize', fontsize_baselinet),...
 axis([-3*pi,3*pi, -4,4])
    axes('position',[0.845 0.05 0.153 0.262]),
plot(xx,abs(p_f3-GG_f3),'color',Color(5,:),'linewidth',1.2);  set(gca, 'fontsize', fontsize_baseline), set(gca, 'XMinorGrid', 'off'),set(gca, 'YMinorGrid', 'off'), grid on; 
title('Error' ,'interpreter','latex', 'fontsize', fontsize_baselinet)
xlabel('$x$ ','interpreter','latex', 'fontsize', fontsize_baselinet),...
%ylabel('Absolute error', 'interpreter','latex','fontsize', fontsize_baselinet),... 
 axis([-3*pi,3*pi, 0,2])
 