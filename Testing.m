% this does not account for temperature effects
clear all;
close all;
clc;
% Define constants:
dt = 0.001; % the step size for the timestep
t = 0:dt:2; % the total duration of the calculations is 2 ms
% Note: any stimulus should be removed after 1.1 ms
dx = .1; % step size for the distance along the axon.
x = -4.5:dx:4.5; % the vector of slices to be used
% 9cm total, centered around 0
Re = .350;
Ri = .110;
Cm = 2.5;
a = 10E-4;
l = 1E-4;
gNa = 1445;
gL = 128;
Ena = 115;
El = -0.01;

I_Stim = 175; %after min threshold for z = 1mm


Ist = zeros(size(t));
Ist(100:250) = -I_Stim; % 0.15ms cathodic pulse

z = 0.1;
% z = 0.2;
% z = 0.4;
% z = 0.8;

Ve = zeros(length(t), length(x));
for i = 1:length(t)
    for j = 1:length(x)
        Ve(i, j) = (Re * Ist(i)) / (4 * pi * sqrt(x(j)^2+z^2));
    end
end
% Define the matrix for membrane current which will only be used for 1
% timestep (prob 3)
im = zeros(1, length(x));
% Define the matrix for longitudinal axonal current which will only be
% used for 1 timestep (prob 3)
i_longitudinal = zeros(1, length(x));
%Define matrices to hold the membrane potential, m and h
V = zeros(length(t), length(x));
m = zeros(length(t), length(x));
h = zeros(length(t), length(x));
%boolean to keep track of the special plots, so we don't repeat them
should_plot_voltage = 1;
should_plot_voltage2 = 1;
index_half_way_to_end = 68;
%loop through time... at each point in time, consider every point along the
%axon in x with a nested inner loop
for i = 1:length(t)
    % now loop through each section of the
    for j = 2:(length(x) - 1)
        Vm = V(i, j); % define a temp variable to hold this timestep's
        % membrane voltage for the current slice in the axon
        % just to help keep things neat later.
        % define the alpha constants based off the formulas provided in the hw
        alpha_m = (97 + 0.363 * Vm) / (1 + exp((31 - Vm)/5.3));
        beta_m = alpha_m / exp((Vm - 23.8)/4.17);
        beta_h = 15.6 / (1 + exp((24 - Vm)/10));
        alpha_h = beta_h / exp((Vm - 5.5)/5);
        % define the change in m for the current timestep, then udpate the
        % next timestep's m value by adding the change
        dmdt = -(alpha_m + beta_m) * m(i, j) + alpha_m;
        m(i+1, j) = m(i, j) + dmdt * dt;
        % define the change in h for the current timestep, then update the
        % next timestep's h value by adding the change
        dhdt = -(alpha_h + beta_h) * h(i, j) + alpha_h;
        h(i+1, j) = h(i, j) + dhdt * dt;
        % define the change in membrane voltage for the current timestep,
        % then update the next timestep's membrane voltage value by adding
        % the change.
        dVdt = (-gNa * m(i, j)^2 * h(i, j) * (Vm - Ena) - gL * (Vm - El) + ((2 * a * dx) / (4 * Ri * l)) * (((V(i, j-1) - 2 * Vm + V(i, j+1)) / dx^2) + ((Ve(i, j-1) - 2 * Ve(i, j) + Ve(i, j+1)) / dx^2))) / Cm;
        V(i+1, j) = V(i, j) + dVdt * dt;
    end
end
% plot the voltage in both time and x... showing the AP.
figure
surf(x, [0, t], V), title('Surface Plot for Bidirectional AP'), xlabel('Distance along axon (cm)'), ylabel('Time (ms)'), zlabel('Voltage (mV)')
shading flat
rotate3d on

% Make a video of the AP.
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