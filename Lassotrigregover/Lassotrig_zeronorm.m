clear all;close all
N=501; 
L=500;
K=randn(1,2);
x=-pi+(2*pi)/N:(2*pi)/N:pi;
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
noise_level=10;
[Y_f1,NOISE_f1]=noisegen(G_f1,noise_level);
[YY_f1,NOISE1_f1]=noisegen(GG_f1,noise_level);


[Y_f2,NOISE_f2]=noisegen(G_f2,noise_level);
[YY_f2,NOISE1_f2]=noisegen(GG_f2,noise_level);



[Y_f3,NOISE_f3]=noisegen(G_f3,noise_level);
[YY_f3,NOISE1_f3]=noisegen(GG_f3,noise_level);



[Y_f4,NOISE_f4]=noisegen(G_f4,noise_level);
[YY_f4,NOISE1_f4]=noisegen(GG_f4,noise_level);


[Y_f5,NOISE_f5]=noisegen(G_f5,noise_level);
[YY_f5,NOISE1_f5]=noisegen(GG_f5,noise_level);

[Y_f6,NOISE_f6]=noisegen(G_f6,noise_level);
[YY_f6,NOISE1_f6]=noisegen(GG_f6,noise_level);


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
   p_f1=beta_f1*C'; p_f2=beta_f2*C'; p_f3=beta_f3*C'; p_f4=beta_f4*C'; p_f5=beta_f5*C';p_f6=beta_f6*C';
   er(k,1)=norm(p_f1-GG_f1);er(k,2)=norm(p_f2-GG_f2);er(k,3)=norm(p_f3-GG_f3);er(k,4)=norm(p_f4-GG_f4);er(k,5)=norm(p_f5-GG_f5);er(k,6)=norm(p_f6-GG_f6);
   end
     for i=1:L+1
ccc(i)=sin(pi*(i/L))/(pi*(i/L));
end
   [ma,t]=min(er);
   lambdaopt_f1=lambda(t(1));lambdaopt_f2=lambda(t(2));lambdaopt_f3=lambda(t(3));lambdaopt_f4=lambda(t(4));lambdaopt_f5=lambda(t(5));lambdaopt_f6=lambda(t(6));
   
   %%  三角多项式系数确定
   beta_f1= l1_beta(w,A,G_f1',lambdaopt_f1,L,ones(L+1));
     Beta_f1= l1_beta(w,A,G_f1',0,L,ones(L+1));

     
   beta_f2= l1_beta(w,A,G_f2',lambdaopt_f2,L,ones(L+1));
     Beta_f2= l1_beta(w,A,G_f2',0,L,ones(L+1));
     
   beta_f3= l1_beta(w,A,G_f3',lambdaopt_f3,L,ones(L+1));
      Beta_f3= l1_beta(w,A,G_f3',0,L,ones(L+1));
      
   beta_f4= l1_beta(w,A,G_f4',lambdaopt_f4,L,ones(L+1));
      Beta_f4= l1_beta(w,A,G_f4',0,L,ones(L+1));
      
   beta_f5= l1_beta(w,A,G_f5',lambdaopt_f5,L,ones(L+1));
     Beta_f5= l1_beta(w,A,G_f5',0,L,ones(L+1));
     
       beta_f6= l1_beta(w,A,G_f6',lambdaopt_f6,L,ones(L+1));
     Beta_f6= l1_beta(w,A,G_f6',0,L,ones(L+1));
     
     
   p_f1=beta_f1*C'; p_f2=beta_f2*C'; p_f3=beta_f3*C'; p_f4=beta_f4*C'; p_f5=beta_f5*C';p_f6=beta_f6*C';
  P_f1=Beta_f1*C'; P_f2=Beta_f2*C'; P_f3=Beta_f3*C'; P_f4=Beta_f4*C'; P_f5=Beta_f5*C';P_f6=Beta_f6*C';
  
  %% 系数0范数
  zeronorm_f1=501-sum(abs(beta_f1)<10^-8)
  zeronorm_F1=501-sum(abs(Beta_f1)<10^-8)
   zeronorm_f2=501-sum(abs(beta_f2)<10^-8)
  zeronorm_F2=501-sum(abs(Beta_f2)<10^-8)
    zeronorm_f3=501-sum(abs(beta_f3)<10^-8)
  zeronorm_F3=501-sum(abs(Beta_f3)<10^-8)
    zeronorm_f4=501-sum(abs(beta_f4)<10^-8)
  zeronorm_F4=501-sum(abs(Beta_f4)<10^-8)
    zeronorm_f5=501-sum(abs(beta_f5)<10^-8)
  zeronorm_F5=501-sum(abs(Beta_f5)<10^-8)
      zeronorm_f6=501-sum(abs(beta_f6)<10^-8)
  zeronorm_F6=501-sum(abs(Beta_f6)<10^-8)
