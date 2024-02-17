clear all;close all
T=5;TT=501;
o=1;
for N=T:6:TT 
L=N-1;
K=randn(1,2);
x=(2*pi)/N:(2*pi)/N:2*pi;
xx=-pi:2*pi/2000:pi;

 W1=ones(1,N);
w=((2*pi)/N)*W1;
OO=chebop(0,2*pi)
OO.op=@(x,u) 0.001*diff(u,2)+0.001*diff(u)-cos(x).*u; 
OO.bc='periodic';
oo=OO\0.1;
%% 五个函数赋值
G_f1=sin(x);
GG_f1=sin(xx);

G_f2=sawtooth(x,1/2);
GG_f2=sawtooth(xx,1/2);



G_f3=sawtooth(x,0.99);
GG_f3=sawtooth(xx,0.99);

G_f4=tanh(50*sin(x));
GG_f4=tanh(50*sin(xx));

G_f5=cos(50*x+4*sin(5*x));
GG_f5=cos(50*xx+4*sin(5*xx));


G_f6=oo(x);
GG_f6=oo(xx);





%% 函数加噪
[Y_f1,NOISE_f1]=noisegen(G_f1,10);
[YY_f1,NOISE1_f1]=noisegen(GG_f1,10);


[Y_f2,NOISE_f2]=noisegen(G_f2,10);
[YY_f2,NOISE1_f2]=noisegen(GG_f2,10);



[Y_f3,NOISE_f3]=noisegen(G_f3,10);
[YY_f3,NOISE1_f3]=noisegen(GG_f3,10);



[Y_f4,NOISE_f4]=noisegen(G_f4,10);
[YY_f4,NOISE1_f4]=noisegen(GG_f4,10);


[Y_f5,NOISE_f5]=noisegen(G_f5,10);
[YY_f5,NOISE1_f5]=noisegen(GG_f5,10);

[Y_f6,NOISE_f6]=noisegen(G_f6,10);
[YY_f6,NOISE1_f6]=noisegen(GG_f6,10);


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
   beta_f1 = l1_beta(w,A,Y_f1',lambda(k),L,ones(L+1));
   beta_f2 = l1_beta(w,A,Y_f2',lambda(k),L,ones(L+1));
   beta_f3 = l1_beta(w,A,Y_f3',lambda(k),L,ones(L+1));
   beta_f4 = l1_beta(w,A,Y_f4',lambda(k),L,ones(L+1));
   beta_f5 = l1_beta(w,A,Y_f5',lambda(k),L,ones(L+1));
    beta_f6 = l1_beta(w,A,Y_f6',lambda(k),L,ones(L+1));
   p_f1=beta_f1*C'; p_f2=beta_f2*C'; p_f3=beta_f3*C'; p_f4=beta_f4*C'; p_f5=beta_f5*C'; p_f6=beta_f6*C';
   er(k,1)=norm(p_f1-GG_f1);er(k,2)=norm(p_f2-GG_f2);er(k,3)=norm(p_f3-GG_f3);er(k,4)=norm(p_f4-GG_f4);er(k,5)=norm(p_f5-GG_f5);er(k,6)=norm(p_f5-GG_f6);
   end
   
   [ma,t]=min(er);
   lambdaopt_f1=lambda(t(1));lambdaopt_f2=lambda(t(2));lambdaopt_f3=lambda(t(3));lambdaopt_f4=lambda(t(4));lambdaopt_f5=lambda(t(5));lambdaopt_f6=lambda(t(6));
   
   %%  三角多项式系数确定
   beta_f1= l1_beta(w,A,Y_f1',lambdaopt_f1,L,ones(L+1));
     Beta_f1= l1_beta(w,A,Y_f1',0,L,ones(L+1));

     
   beta_f2= l1_beta(w,A,Y_f2',lambdaopt_f2,L,ones(L+1));
     Beta_f2= l1_beta(w,A,Y_f2',0,L,ones(L+1));
     
   beta_f3= l1_beta(w,A,Y_f3',lambdaopt_f3,L,ones(L+1));
      Beta_f3= l1_beta(w,A,Y_f3',0,L,ones(L+1));
      
   beta_f4= l1_beta(w,A,Y_f4',lambdaopt_f4,L,ones(L+1));
      Beta_f4= l1_beta(w,A,Y_f4',0,L,ones(L+1));
      
   beta_f5= l1_beta(w,A,Y_f5',lambdaopt_f5,L,ones(L+1));
     Beta_f5= l1_beta(w,A,Y_f5',0,L,ones(L+1));
     
        beta_f6= l1_beta(w,A,Y_f6',lambdaopt_f6,L,ones(L+1));
     Beta_f6= l1_beta(w,A,Y_f6',0,L,ones(L+1));
     
   p_f1=beta_f1*C'; p_f2=beta_f2*C'; p_f3=beta_f3*C'; p_f4=beta_f4*C'; p_f5=beta_f5*C';p_f6=beta_f6*C';
  P_f1=Beta_f1*C'; P_f2=Beta_f2*C'; P_f3=Beta_f3*C'; P_f4=Beta_f4*C'; P_f5=Beta_f5*C';P_f6=Beta_f6*C';
  

  
  %% 相对二范数误差
  related_noise_f1(o)=norm(GG_f1-YY_f1)/norm(GG_f1);
  related_error_F1(o)=norm(GG_f1-P_f1)/norm(GG_f1);
  related_error_f1(o)=norm(GG_f1-p_f1)/norm(GG_f1);
  
    related_noise_f2(o)=norm(GG_f2-YY_f2)/norm(GG_f2);
  related_error_F2(o)=norm(GG_f2-P_f2)/norm(GG_f2);
  related_error_f2(o)=norm(GG_f2-p_f2)/norm(GG_f2);
  
    related_noise_f3(o)=norm(GG_f3-YY_f3)/norm(GG_f3);
  related_error_F3(o)=norm(GG_f3-P_f3)/norm(GG_f3);
  related_error_f3(o)=norm(GG_f3-p_f3)/norm(GG_f3);
  
    related_noise_f4(o)=norm(GG_f4-YY_f4)/norm(GG_f4);
  related_error_F4(o)=norm(GG_f4-P_f4)/norm(GG_f4);
  related_error_f4(o)=norm(GG_f4-p_f4)/norm(GG_f4);
  
    related_noise_f5(o)=norm(GG_f5-YY_f5)/norm(GG_f5);
  related_error_F5(o)=norm(GG_f5-P_f5)/norm(GG_f5);
  related_error_f5(o)=norm(GG_f5-p_f5)/norm(GG_f5);
  
   related_noise_f6(o)=norm(GG_f6-YY_f6)/norm(GG_f6);
  related_error_F6(o)=norm(GG_f6-P_f6)/norm(GG_f6);
  related_error_f6(o)=norm(GG_f6-p_f6)/norm(GG_f6);
  o=o+1;
end

Color = [215,25,28;
0 0 128;
254,204,92;
102, 0, 204;
255,255,255;
255*0.8500 255*0.3250 255*0.0980;
255*0.4940 255*0.1840 255*0.5560]/255;
Markersize=3;
fontsize_baseline = 15;
fontsize_baselinet = 20;
fontsize_baselinea = 15;
figure(1)
axes('position',[0.05 0.55 0.28 0.38]), 
plot(T:6:TT,related_error_F1,'o-','color',Color(2,:),'MarkerFaceColor', Color(2,:),'MarkerSize',Markersize,'linewidth',1.2), box on,...
   hold on; 
plot(T:6:TT,related_error_f1,'V-','color',  Color(1,:),'MarkerFaceColor', Color(1,:),'MarkerSize',Markersize,'linewidth',1.2),box on,...;
    plot(T:6:TT, related_noise_f1,'--','color', Color(7,:),'MarkerFaceColor', Color(2,:),'MarkerSize',Markersize,'linewidth',2),box on,...;
 legend('$\mathcal{T}_n f $','$\mathcal{T}_n^{\lambda} f $','noise level','Interpreter','latex','Fontsize',12,'Location','Northeast','fontsize', fontsize_baseline)
    set(gca, 'fontsize', fontsize_baselinea),...
    xlabel('$N$','interpreter','latex', 'fontsize', fontsize_baseline),
   ylabel(' Relative error','interpreter','latex', 'fontsize', fontsize_baseline),...
    title('Relative error for approximating $f_1(x)$','interpreter','latex', 'fontsize', fontsize_baselinet),...
     grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),axis([0,501,0,0.45]),
 
axes('position',[0.375 0.55 0.28 0.38]),
plot(T:6:TT,related_error_F2,'o-','color',Color(2,:),'MarkerFaceColor', Color(2,:),'MarkerSize',Markersize,'linewidth',1.2),box on,...
    hold on; 
plot(T:6:TT,related_error_f2,'V-','color',  Color(1,:),'MarkerFaceColor', Color(1,:),'MarkerSize',Markersize,'linewidth',1.2),box on,...;
    plot(T:6:TT, related_noise_f2,'--','color',  Color(7,:),'MarkerFaceColor', Color(6,:),'MarkerSize',Markersize,'linewidth',2),box on,...;
     legend('$\mathcal{T}_n f $','$\mathcal{T}_n^{\lambda} f $','noise level','Interpreter','latex','Fontsize',12,'fontsize', fontsize_baseline)
    set(gca, 'fontsize', fontsize_baselinea),...
    xlabel('$N$','interpreter','latex', 'fontsize', fontsize_baseline),...
    ylabel('Relative error','interpreter','latex', 'fontsize', fontsize_baseline),...
    title('Relative error for approximating $f_2(x)$','interpreter','latex', 'fontsize', fontsize_baselinet),...
     grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),axis([0,501,0,0.45]),
 
axes('position',[0.71 0.55 0.28 0.38]),
plot(T:6:TT,related_error_F3,'o-','color',Color(2,:),'MarkerFaceColor', Color(2,:),'MarkerSize',Markersize,'linewidth',1.2),box on,...
    hold on; 
plot(T:6:TT,related_error_f3,'V-','color',  Color(1,:),'MarkerFaceColor', Color(1,:),'MarkerSize',Markersize,'linewidth',1.2),box on,...;
   plot(T:6:TT, related_noise_f3,'--','color',  Color(7,:),'MarkerFaceColor', Color(6,:),'MarkerSize',Markersize,'linewidth',2),box on,...;
     legend('$\mathcal{T}_n f $','$\mathcal{T}_n^{\lambda} f $','noise level','Interpreter','latex','Fontsize',12,'Location','Northeast','fontsize', fontsize_baseline)
    set(gca, 'fontsize', fontsize_baselinea),...
    xlabel('$N$','interpreter','latex', 'fontsize', fontsize_baseline),...
    ylabel('Relative error','interpreter','latex', 'fontsize', fontsize_baseline),...
    title('Relative error for approximating $f_3(x)$','interpreter','latex', 'fontsize', fontsize_baselinet),...
     grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),axis([0,501,0,1.2]),
 
 
axes('position',[0.05 0.06 0.28 0.38]), 
plot(T:6:TT,related_error_F4,'o-','color',Color(2,:),'MarkerFaceColor', Color(2,:),'MarkerSize',Markersize,'linewidth',1.2), box on,...
   hold on; 
plot(T:6:TT,related_error_f4,'V-','color',  Color(1,:),'MarkerFaceColor', Color(1,:),'MarkerSize',Markersize,'linewidth',1.2),box on,...;
    plot(T:6:TT, related_noise_f4,'--','color',  Color(7,:),'MarkerFaceColor', Color(6,:),'MarkerSize',Markersize,'linewidth',2),box on,...;
     legend('$\mathcal{T}_n f $','$\mathcal{T}_n^{\lambda} f $','noise level','Interpreter','latex','Fontsize',12,'fontsize', fontsize_baseline)
    set(gca, 'fontsize', fontsize_baselinea),...
    xlabel('$N$','interpreter','latex', 'fontsize', fontsize_baseline),
   ylabel('Relative error','interpreter','latex', 'fontsize', fontsize_baseline),...
    title('Relative error for approximating $f_4(x)$','interpreter','latex', 'fontsize', fontsize_baselinet),...
     grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),axis([0,501,0,0.6]),
 
 axes('position',[0.375 0.06 0.28 0.38]), 
plot(T:6:TT,related_error_F5,'o-','color',Color(2,:),'MarkerFaceColor', Color(2,:),'MarkerSize',Markersize,'linewidth',1.2), box on,...
   hold on; 
plot(T:6:TT,related_error_f5,'V-','color',  Color(1,:),'MarkerFaceColor', Color(1,:),'MarkerSize',Markersize,'linewidth',1.2),box on,...;
    plot(T:6:TT, related_noise_f5,'--','color',  Color(7,:),'MarkerFaceColor', Color(6,:),'MarkerSize',Markersize,'linewidth',2),box on,...;
     legend('$\mathcal{T}_n f $','$\mathcal{T}_n^{\lambda} f $','noise level','Interpreter','latex','Fontsize',12,'fontsize', fontsize_baseline)
    set(gca, 'fontsize', fontsize_baselinea),...
    xlabel('$N$','interpreter','latex', 'fontsize', fontsize_baseline),
   ylabel('Relative error','interpreter','latex', 'fontsize', fontsize_baseline),...
    title('Relative error for approximating $f_5(x)$','interpreter','latex', 'fontsize', fontsize_baselinet),...
     grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),axis([0,501,0,2]),
 
  axes('position',[0.71 0.06 0.28 0.38]), 
plot(T:6:TT,related_error_F6,'o-','color',Color(2,:),'MarkerFaceColor', Color(2,:),'MarkerSize',Markersize,'linewidth',1.2), box on,...
   hold on; 
plot(T:6:TT,related_error_f6,'V-','color',  Color(1,:),'MarkerFaceColor', Color(1,:),'MarkerSize',Markersize,'linewidth',1.2),box on,...;
    plot(T:6:TT, related_noise_f6,'--','color',  Color(7,:),'MarkerFaceColor', Color(6,:),'MarkerSize',Markersize,'linewidth',2),box on,...;
     legend('$\mathcal{T}_n f $','$\mathcal{T}_n^{\lambda} f $','noise level','Interpreter','latex','Fontsize',12,'fontsize', fontsize_baseline)
    set(gca, 'fontsize', fontsize_baselinea),...
    xlabel('$N$','interpreter','latex', 'fontsize', fontsize_baseline),
   ylabel('Relative error','interpreter','latex', 'fontsize', fontsize_baseline),...
    title('Relative error for approximating $f_6(x)$','interpreter','latex', 'fontsize', fontsize_baselinet),...
     grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),axis([0,501,0,1.5]),
 