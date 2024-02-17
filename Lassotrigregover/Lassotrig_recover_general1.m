clear all;close all
N=501;
N1=501;
L=500;
xx=-3*pi:6*pi/2000:3*pi;
K=randn(1,2);



x=-pi+(2*pi)/N:(2*pi)/N:pi;
x1=-pi+(2*pi)/N1:(2*pi)/N1:pi;
 W1=ones(1,N);
w=((2*pi)/N)*W1;
   %% f_1(x) add noise      
G_f1= sin(x);
G1_f1= sin(x1);
GG_f1=  sin(xx);
[Y_f1,NOISE_f1]=noisegen(G_f1,10);
[YY_f1,NOISE1_f1]=noisegen(G1_f1,10);
[yy_f1,NOISE2_f1]=noisegen(GG_f1,10);
delta_f1=norm(NOISE2_f1);

%% f_2(x) add noise
G_f2= sawtooth(x,1/2);
G1_f2=sawtooth(x1,1/2);
GG_f2= sawtooth(xx,1/2);
[Y_f2,NOISE_f2]=noisegen(G_f2,10);
[YY_f2,NOISE1_f2]=noisegen(G1_f2,10);
[yy_f2,NOISE2_f2]=noisegen(GG_f2,10);
delta_f2=norm(NOISE2_f2);

 %% f_3(x) with h=3 add noise 
G_f3=sawtooth(x,0.99); 
GG_f3=sawtooth(xx,0.99); 
G1_f3=sawtooth(x1,0.99); 
[Y_f3,NOISE_f3]=noisegen(G_f3,10);
[YY_f3,NOISE1_f3]=noisegen(G1_f3,10);
[yy_f3,NOISE2_f3]=noisegen(GG_f3,10);
delta_f3=norm(NOISE2_f3);

 %% f_4(x) add noise
G_f4=tanh(50*sin(x));
GG_f4=tanh(50*sin(xx));
G1_f4=tanh(50*sin(x1));
[Y_f4,NOISE_f4]=noisegen(G_f4,10);
[YY_f4,NOISE1_f4]=noisegen(G1_f4,10);
[yy_f4,NOISE2_f4]=noisegen(GG_f4,10);
delta_f4=norm(NOISE2_f4);


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
    beta_f4 = l1_beta(w,A1,Y_f4',lambda(k),L,ones(L+1));
   p_f1=beta_f1*C'; p_f2=beta_f2*C'; p_f3=beta_f3*C'; p_f4=beta_f4*C'
   er(k,1)=norm(p_f1-GG_f1);er(k,2)=norm(p_f2-GG_f2);er(k,3)=norm(p_f3-GG_f3);er(k,4)=norm(p_f4-GG_f4);
   end
     for i=1:L+1
ccc(i)=sin(pi*(i/L))/(pi*(i/L));
end
   [ma,t]=min(er);
   lambda_f1=lambda(t(1));lambda_f2=lambda(t(2));lambda_f3=lambda(t(3));lambda_f4=lambda(t(4));



%% determin regularization solution
beta_f1 = l1_beta(w,A1,Y_f1',lambda_f1,L,ones(L+1));
Beta_f1 = l1_beta(w,A1,Y_f1',0,L,ones(L+1));
beta_f2 = l1_beta(w,A1,Y_f2',lambda_f2,L,ones(L+1));
Beta_f2 = l1_beta(w,A1,Y_f2',0,L,ones(L+1));
beta_f3 = l1_beta(w,A1,Y_f3',lambda_f3,L,ones(L+1));
Beta_f3 = l1_beta(w,A1,Y_f3',0,L,ones(L+1));
beta_f4 = l1_beta(w,A1,Y_f4',lambda_f4,L,ones(L+1));
Beta_f4 = l1_beta(w,A1,Y_f4',0,L,ones(L+1));

for i=1:50
    t1=beta_f1*(diag(ccc))^(i);
   Lf1(i)=norm(t1*C'-GG_f1);
end
[k1 pf1]=min(Lf1);
Sbeta_f1=beta_f1*(diag(ccc))^(pf1);


for i=1:50
    t2=beta_f2*(diag(ccc))^(i);
   Lf2(i)=norm(t2*C'-GG_f2);
end
[k2 pf2]=min(Lf2);
Sbeta_f2=beta_f2*(diag(ccc))^(pf2);

for i=1:50
    t3=beta_f3*(diag(ccc))^(i);
   Lf3(i)=norm(t3*C'-GG_f3);
end
[k3 pf3]=min(Lf3);
Sbeta_f3=beta_f3*(diag(ccc))^(pf3);


for i=1:50
    t4=beta_f4*(diag(ccc))^(i);
   Lf4(i)=norm(t4*C'-GG_f4);
end
[k4 pf4]=min(Lf4);
Sbeta_f4=beta_f4*(diag(ccc))^(pf4);

%% approximation trigonometric polynomial
 p_f1=beta_f1*C'; er1=abs(p_f1-GG_f1);
 P_f1=Beta_f1*C'; Er1=abs(yy_f1-GG_f1);
 p_f2=beta_f2*C'; er2=abs(p_f2-GG_f2);
 P_f2=Beta_f2*C'; Er2=abs(yy_f2-GG_f2);
 p_f3=beta_f3*C'; er3=abs(p_f3-GG_f3);
 P_f3=Beta_f3*C'; Er3=abs(yy_f3-GG_f3);
 p_f4=beta_f4*C'; er4=abs(p_f4-GG_f4);
 P_f4=Beta_f4*C'; Er4=abs(yy_f4-GG_f4);
%% plot

fontsize_baselinet=12;
fontsize_baseline=10;
Color = [215,25,28;
0 0 128;
254,204,92;
171,217,233;
44,123,182;
255*0.8500 255*0.3250 255*0.0980;
255*0.4940 255*0.1840 255*0.5560]/255;
interspace=5;
%图第一行

    axes('position',[0.015 0.79 0.23 0.175]),
     
 plot(xx,GG_f1,'color',Color(2,:),'linewidth',1.2); set(gca, 'fontsize', fontsize_baseline),set(gca, 'XMinorGrid', 'off'),set(gca, 'YMinorGrid', 'off'), grid on;
 title('Exact $\sin$ wave' ,'interpreter','latex', 'fontsize', fontsize_baselinet)
xlabel('$x$ ','interpreter','latex', 'fontsize', fontsize_baselinet),...
% ylabel('$f(x)$', 'interpreter','latex','fontsize', fontsize_baselinet),...
 axis([-3*pi,3*pi, -1.7,1.7])
 

 axes('position',[0.265 0.79 0.23 0.175]),
 plot(xx,yy_f1,'color',Color(2,:),'linewidth',1.2); set(gca, 'fontsize', fontsize_baseline),set(gca, 'XMinorGrid', 'off'),set(gca, 'YMinorGrid', 'off'), grid on;
  title('${\rm Sin}$ wave with 10 dB noise ' ,'interpreter','latex', 'fontsize', fontsize_baselinet)
xlabel('$x$ ','interpreter','latex', 'fontsize', fontsize_baselinet),...
% ylabel('$f(x)$', 'interpreter','latex','fontsize', fontsize_baselinet),...
 axis([-3*pi,3*pi, -1.7,1.7])
 
 
 
   axes('position',[0.515 0.79 0.23 0.175]),
 plot(xx,p_f1,'color',Color(2,:),'linewidth',1.2); set(gca, 'fontsize', fontsize_baseline),set(gca, 'XMinorGrid', 'off'),set(gca, 'YMinorGrid', 'off'), grid on;
 title('Recover $\sin$ wave by LTI' ,'interpreter','latex', 'fontsize', fontsize_baselinet)
xlabel('$x$ ','interpreter','latex', 'fontsize', fontsize_baselinet),...
% ylabel('$f(x)$', 'interpreter','latex','fontsize', fontsize_baselinet),...
 axis([-3*pi,3*pi, -1.7,1.7])
 

axes('position',[0.765 0.79 0.23 0.175]),
 plot(xx(1:interspace:length(xx)),Er1(1:interspace:length(Er1)),'--','color',Color(1,:),'linewidth',0.01); 
 hold on
 plot(xx(1:interspace:length(xx)),er1(1:interspace:length(er1)),'color',Color(7,:),'linewidth',1.2); 
 title('Errors' ,'interpreter','latex', 'fontsize', fontsize_baselinet)
 set(gca, 'fontsize', fontsize_baseline),set(gca, 'XMinorGrid', 'off'),set(gca, 'YMinorGrid', 'off'), grid on; 
legend('No regularization','LTI with optimal $\lambda$','Interpreter','latex','Fontsize',12,'fontsize', fontsize_baseline)
xlabel('$x$ ','interpreter','latex', 'fontsize', fontsize_baselinet),...
% ylabel('$f(x)$', 'interpreter','latex','fontsize', fontsize_baselinet),...
axis([-3*pi,3*pi, 0,0.8])
 
 %图第二行
 axes('position',[0.015 0.545 0.23 0.175]),
 plot(xx,GG_f2,'color',Color(2,:),'linewidth',1.2); set(gca, 'fontsize', fontsize_baseline),set(gca, 'XMinorGrid', 'off'),set(gca, 'YMinorGrid', 'off'), grid on; 
 title('Exact triangular wave' ,'interpreter','latex', 'fontsize', fontsize_baselinet)
xlabel('$x$ ','interpreter','latex', 'fontsize', fontsize_baselinet),...
% ylabel('$f(x)$', 'interpreter','latex','fontsize', fontsize_baselinet),...
 axis([-3*pi,3*pi, -1.7,1.7])
 
 axes('position',[0.265 0.545 0.23 0.175]),
 plot(xx,yy_f2,'color',Color(2,:),'linewidth',1.2); set(gca, 'fontsize', fontsize_baseline),set(gca, 'XMinorGrid', 'off'),set(gca, 'YMinorGrid', 'off'), grid on;
  title('Triangular wave with 10 dB noise ' ,'interpreter','latex', 'fontsize', fontsize_baselinet)
xlabel('$x$ ','interpreter','latex', 'fontsize', fontsize_baselinet),...
% ylabel('$f(x)$', 'interpreter','latex','fontsize', fontsize_baselinet),...
 axis([-3*pi,3*pi, -1.7,1.7])
 
 
 
   axes('position',[0.515 0.545 0.23 0.175]),
 plot(xx,p_f2,'color',Color(2,:),'linewidth',1.2); set(gca, 'fontsize', fontsize_baseline),set(gca, 'XMinorGrid', 'off'),set(gca, 'YMinorGrid', 'off'), grid on;
 title('Recover triangular wave by LTI' ,'interpreter','latex', 'fontsize', fontsize_baselinet)
xlabel('$x$ ','interpreter','latex', 'fontsize', fontsize_baselinet),...
% ylabel('$f(x)$', 'interpreter','latex','fontsize', fontsize_baselinet),...
 axis([-3*pi,3*pi, -1.7,1.7])
 

axes('position',[0.765 0.545 0.23 0.175]),
 plot(xx(1:interspace:length(xx)),Er2(1:interspace:length(Er2)),'--','color',Color(1,:),'linewidth',0.01); 
 hold on
 plot(xx(1:interspace:length(xx)),er2(1:interspace:length(er2)),'color',Color(7,:),'linewidth',1.2); 
 title('Errors' ,'interpreter','latex', 'fontsize', fontsize_baselinet)
 set(gca, 'fontsize', fontsize_baseline),set(gca, 'XMinorGrid', 'off'),set(gca, 'YMinorGrid', 'off'), grid on; 
legend('No regularization','LTI with optimal $\lambda$','Interpreter','latex','Fontsize',12,'fontsize', fontsize_baseline)
xlabel('$x$ ','interpreter','latex', 'fontsize', fontsize_baselinet),...
% ylabel('$f(x)$', 'interpreter','latex','fontsize', fontsize_baselinet),...
axis([-3*pi,3*pi, 0,0.8])
 
 
  %图第三行
 axes('position',[0.015 0.3 0.23 0.175]),
 plot(xx,GG_f3,'color',Color(2,:),'linewidth',1.2); set(gca, 'fontsize', fontsize_baseline),set(gca, 'XMinorGrid', 'off'),set(gca, 'YMinorGrid', 'off'), grid on; 
 title('Exact sawtooth wave' ,'interpreter','latex', 'fontsize', fontsize_baselinet)
xlabel('$x$ ','interpreter','latex', 'fontsize', fontsize_baselinet),...
% ylabel('$f(x)$', 'interpreter','latex','fontsize', fontsize_baselinet),...
 axis([-3*pi,3*pi, -1.7,1.7])
 
 
  axes('position',[0.265 0.3 0.23 0.175]),
 plot(xx,yy_f3,'color',Color(2,:),'linewidth',1.2); set(gca, 'fontsize', fontsize_baseline),set(gca, 'XMinorGrid', 'off'),set(gca, 'YMinorGrid', 'off'), grid on;
  title('Sawtooth wave with 10 dB noise ' ,'interpreter','latex', 'fontsize', fontsize_baselinet)
xlabel('$x$ ','interpreter','latex', 'fontsize', fontsize_baselinet),...
% ylabel('$f(x)$', 'interpreter','latex','fontsize', fontsize_baselinet),...
 axis([-3*pi,3*pi, -1.7,1.7])
 
 

 
   axes('position',[0.515 0.3 0.23 0.175]),
 plot(xx,p_f3,'color',Color(2,:),'linewidth',1.2); set(gca, 'fontsize', fontsize_baseline),set(gca, 'XMinorGrid', 'off'),set(gca, 'YMinorGrid', 'off'), grid on;
 title('Recover sawtooth wave by LTI' ,'interpreter','latex', 'fontsize', fontsize_baselinet)
xlabel('$x$ ','interpreter','latex', 'fontsize', fontsize_baselinet),...
% ylabel('$f(x)$', 'interpreter','latex','fontsize', fontsize_baselinet),...
 axis([-3*pi,3*pi, -1.7,1.7])
 

axes('position',[0.765 0.3 0.23 0.175]),
 plot(xx(1:interspace:length(xx)),Er3(1:interspace:length(Er3)),'--','color',Color(1,:),'linewidth',0.01); 
 hold on
 plot(xx(1:interspace:length(xx)),er3(1:interspace:length(er3)),'color',Color(7,:),'linewidth',1.2); 
 title('Errors' ,'interpreter','latex', 'fontsize', fontsize_baselinet)
 set(gca, 'fontsize', fontsize_baseline),set(gca, 'XMinorGrid', 'off'),set(gca, 'YMinorGrid', 'off'), grid on; 
legend('No regularization','LTI with optimal $\lambda$','Interpreter','latex','Fontsize',12,'fontsize', fontsize_baseline)
xlabel('$x$ ','interpreter','latex', 'fontsize', fontsize_baselinet),...
% ylabel('$f(x)$', 'interpreter','latex','fontsize', fontsize_baselinet),...
axis([-3*pi,3*pi, 0,1])
 
 
   %图第四行
 axes('position',[0.015 0.055 0.23 0.175]),
 plot(xx,GG_f4,'color',Color(2,:),'linewidth',1.2); set(gca, 'fontsize', fontsize_baseline),set(gca, 'XMinorGrid', 'off'),set(gca, 'YMinorGrid', 'off'), grid on; 
 title('Exact square wave' ,'interpreter','latex', 'fontsize', fontsize_baselinet)
xlabel('$x$ ','interpreter','latex', 'fontsize', fontsize_baselinet),...
% ylabel('$f(x)$', 'interpreter','latex','fontsize', fontsize_baselinet),...
 axis([-3*pi,3*pi, -1.7,1.7])
 
 
  
  axes('position',[0.265 0.055 0.23 0.175]),
 plot(xx,yy_f4,'color',Color(2,:),'linewidth',1.2); set(gca, 'fontsize', fontsize_baseline),set(gca, 'XMinorGrid', 'off'),set(gca, 'YMinorGrid', 'off'), grid on;
  title('Square wave with 10 dB noise ' ,'interpreter','latex', 'fontsize', fontsize_baselinet)
xlabel('$x$ ','interpreter','latex', 'fontsize', fontsize_baselinet),...
% ylabel('$f(x)$', 'interpreter','latex','fontsize', fontsize_baselinet),...
 axis([-3*pi,3*pi, -1.7,1.7])
 

 
   axes('position',[0.515 0.055 0.23 0.175]),
 plot(xx,p_f4,'color',Color(2,:),'linewidth',1.2); set(gca, 'fontsize', fontsize_baseline),set(gca, 'XMinorGrid', 'off'),set(gca, 'YMinorGrid', 'off'), grid on;
 title('Recover square wave by LTI' ,'interpreter','latex', 'fontsize', fontsize_baselinet)
xlabel('$x$ ','interpreter','latex', 'fontsize', fontsize_baselinet),...
% ylabel('$f(x)$', 'interpreter','latex','fontsize', fontsize_baselinet),...
 axis([-3*pi,3*pi, -1.7,1.7])
 

axes('position',[0.765 0.055 0.23 0.175]),
 plot(xx(1:interspace:length(xx)),Er4(1:interspace:length(Er4)),'--','color',Color(1,:),'linewidth',0.01); 
 hold on
 plot(xx(1:interspace:length(xx)),er4(1:interspace:length(er4)),'color',Color(7,:),'linewidth',1.2); 
 title('Errors' ,'interpreter','latex', 'fontsize', fontsize_baselinet)
 set(gca, 'fontsize', fontsize_baseline),set(gca, 'XMinorGrid', 'off'),set(gca, 'YMinorGrid', 'off'), grid on; 
legend('No regularization','LTI with optimal $\lambda$','Interpreter','latex','Fontsize',12,'fontsize', fontsize_baseline)
xlabel('$x$ ','interpreter','latex', 'fontsize', fontsize_baselinet),...
% ylabel('$f(x)$', 'interpreter','latex','fontsize', fontsize_baselinet),...
axis([-3*pi,3*pi, 0,1])
 