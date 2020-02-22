clear all;
close all;
clc;

dt = 0.001; 
t = 0:dt:2; 
dx = .1; 
x = -4.5:dx:4.5;
Re = .350;
Ri = .110;
Cm = 2.5;
a = 10E-4;
l = 1E-4;
gNa = 1445;
gL = 128;
Ena = 115;
El = -0.01;

I_Stim = 175; 

Ist = zeros(size(t));
Ist(100:250) = -I_Stim; 
z = 0.1;
%z = 0.2;
%z = 0.4;
%z = 0.8;

Ve = zeros(length(t), length(x));
for i = 1:length(t)
    for j = 1:length(x)
        Ve(i, j) = (Re * Ist(i)) / (4 * pi * sqrt(x(j)^2+z^2));
    end
end

im = zeros(1, length(x));

i_longitudinal = zeros(1, length(x));

V = zeros(length(t), length(x));
m = zeros(length(t), length(x));
h = zeros(length(t), length(x));


for i = 1:length(t)
    for j = 2:(length(x) - 1)
        Vm = V(i, j); 
        alpha_m = (97 + 0.363 * Vm) / (1 + exp((31 - Vm)/5.3));
        beta_m = alpha_m / exp((Vm - 23.8)/4.17);
        beta_h = 15.6 / (1 + exp((24 - Vm)/10));
        alpha_h = beta_h / exp((Vm - 5.5)/5);
        dmdt = -(alpha_m + beta_m) * m(i, j) + alpha_m;
        m(i+1, j) = m(i, j) + dmdt * dt;
        dhdt = -(alpha_h + beta_h) * h(i, j) + alpha_h;
        h(i+1, j) = h(i, j) + dhdt * dt;
        dVdt = (-gNa * m(i, j)^2 * h(i, j) * (Vm - Ena) - gL * (Vm - El) + ((2 * a * dx) / (4 * Ri * l)) * (((V(i, j-1) - 2 * Vm + V(i, j+1)) / dx^2) + ((Ve(i, j-1) - 2 * Ve(i, j) + Ve(i, j+1)) / dx^2))) / Cm;
        V(i+1, j) = V(i, j) + dVdt * dt;
    end
end

figure
surf(x, [0, t], V), title('Surface Plot for Bidirectional AP'), xlabel('Distance along axon (cm)'), ylabel('Time (ms)'), zlabel('Voltage (mV)')
shading flat
rotate3d on

% Make a video
% % video = VideoWriter('BidirectionalAP: Z=1mm.mp4');
% % open(video);
% % figure
% % for i = 1:length(t)
% %     plot(x, V(i, :));
% %     axis([-4.5, 4.5, -25, 100]);
% %     title('Bidirectional Action Potential'), xlabel('Distance along axon (cm)'), ylabel('Voltage (mV)');
% %     mov(i) = getframe(gcf);
% %     writeVideo(video, mov(i));
% % end