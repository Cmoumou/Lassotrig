clear all;close all



xx=-4*pi:8*pi/2000:4*pi;
y_sin=sin(xx);
y_square=square(xx);
y_triangular=sawtooth(xx,1/2);
y_sawtooth=sawtooth(xx);
fontsize_baselinet=12;
fontsize_baseline=10;
Color = [215,25,28;
0 0 128;
254,204,92;
171,217,233;
44,123,182;
0.4940*255 0.1840*255 0.5560*255]/255;

          axes('position',[0.135 0.41 0.16 0.34]),
plot(xx,y_sin,'color',Color(6,:),'linewidth',1.2); set(gca, 'fontsize', fontsize_baseline),set(gca, 'XMinorGrid', 'off'),set(gca, 'YMinorGrid', 'off'), grid on; 
 title('Sin wave' ,'interpreter','latex', 'fontsize', fontsize_baselinet)
xlabel('$x$ ','interpreter','latex', 'fontsize', fontsize_baselinet),...
     axis([-4*pi,4*pi, -1.2,1.2])
 
        axes('position',[0.325 0.41 0.16 0.34]),
plot(xx,y_triangular,'color',Color(6,:),'linewidth',1.2); set(gca, 'fontsize', fontsize_baseline),set(gca, 'XMinorGrid', 'off'),set(gca, 'YMinorGrid', 'off'), grid on; 
 title('Triangular wave' ,'interpreter','latex', 'fontsize', fontsize_baselinet)
xlabel('$x$ ','interpreter','latex', 'fontsize', fontsize_baselinet),...
     axis([-4*pi,4*pi, -1.2,1.2])
 
 
 
  axes('position',[0.515 0.41 0.16 0.34]),
plot(xx,y_square,'color',Color(6,:),'linewidth',1.2); set(gca, 'fontsize', fontsize_baseline),set(gca, 'XMinorGrid', 'off'),set(gca, 'YMinorGrid', 'off'), grid on; 
 title('Square wave' ,'interpreter','latex', 'fontsize', fontsize_baselinet)
xlabel('$x$ ','interpreter','latex', 'fontsize', fontsize_baselinet),...
     axis([-4*pi,4*pi, -1.2,1.2])
 
  axes('position',[0.705 0.41 0.16 0.34]),
plot(xx,y_sawtooth,'color',Color(6,:),'linewidth',1.2); set(gca, 'fontsize', fontsize_baseline),set(gca, 'XMinorGrid', 'off'),set(gca, 'YMinorGrid', 'off'), grid on; 
 title('Sawtooth wave' ,'interpreter','latex', 'fontsize', fontsize_baselinet)
xlabel('$x$ ','interpreter','latex', 'fontsize', fontsize_baselinet),...
     axis([-4*pi,4*pi, -1.2,1.2])
 
