% this does not account for temperature effects
clear all;
close all;
clc;
% Define constants:
dt = 0.001; % the step size for the timestep
t = 0:dt:2; % the total duration of the calculations is 2 ms

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
%%%% vary the values for the cathodic pulse current ampltitude until the
%%%% surface plot shows a generated AP... hone in a threshold.
% cathodic pulse current amplitude = 144; %before min threshold for z = 1mm
cathodic_pulse_current_amplitude = 175; %after min threshold for z = 1mm
% cathodic pulse current amplitude = 500; %before min threshold for z = 2mm
% cathodic pulse current amplitude = 501; %after min threshold for z = 2mm
% cathodic pulse current amplitude = 2113; %before min threshold for z = 4mm
% cathodic pulse current amplitude = 2114; %after min threshold for z = 4mm
% cathodic pulse current amplitude = 11242; %before min threshold for z = 8mm
% cathodic pulse current amplitude = 11243; %after min threshold for z = 8mm
% Define the stimulation current as a vector the same length as time. To
% make a simple square pulse, a constant value for a fixed time is sufficient.
% Include an anodic pulse from 1ms to 9ms (8ms total)
% Include a cathodic pulse from 9ms to 11ms (2ms total)
% The anodic pulse should be 1/4 the magnitude of the cathodic pulse.
Ist = zeros(size(t));
%Ist(100:900) = 0.25 * cathodic_pulse_current_amplitude; % 0.8ms anodic pulse
%Ist(900:1100) = -cathodic_pulse_current_amplitude; % 0.2ms cathodic pulse

Ist(100:250) = -cathodic_pulse_current_amplitude; % 0.15ms cathodic pulse


% z represents the perpendicular distance between the axon and the
% stimulating electrode. (make sure that the corresponding cathodic current
% ampltitude is used for the right distance "z" or the script will error out.
z = 0.1;
% z = 0.2;
% z = 0.4;
% z = 0.8;
% Define a matrix to quantify the Extracellular Potential at every point
% along the length of the axon "length(x)" for every instant in time
% "length(t)"

Ist1 = zeros(size(t));
Ist2 = zeros(size(t));
Ist1(100:200) = -300;
Ist2(100:500) = 300;
Ve = zeros(length(t), length(x));
z = 0.1;
for i = 1:length(t)
    for j = 1:length(x)
        Ve1 = (Re * Ist1(i)) / (4 * pi * sqrt(x(j)^2+z^2));
        Ve2 = (Re * Ist2(i)) / (4 * pi * sqrt((x(j) + 0.1)^2+z^2)); %only propogate to right

        %   Ve2 = (Re * Ist2(i)) / (4 * pi * sqrt((x(j) - 0.1)ˆ2+zˆ2)); %only propogate to left
        Ve(i, j) = Ve1 + Ve2;
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
    %Plot the voltage along the axon at the instant the AP propogates
    %halfway to the ends of the axon.
    if (V(i, index_half_way_to_end) > 1)
        if (V(i, index_half_way_to_end) > 90 && should_plot_voltage2 == 1)
            %plot membrane voltage
            figure;
            subplot(1, 3, 1);
            plot(x, V(i, :));
            grid on;
            xlabel('Distance');
            ylabel('Membrane Voltage');
            title('Membrane Voltage vs Distance');
            %plot membrane current
            subplot(1, 3, 2);

            for j = 2:(length(x) - 1)
                im(j) = -gNa * m(i, j)^2 * h(i, j) * (Vm - Ena) - gL * (Vm - El);
            end
            
            for j = 2:(length(x) - 1)
                i_longitudinal(j) = ((2 * a * dx) / (4 * Ri * l)) * ((V(i, j-1) - 2 * Vm + V(i, j+1)) / dx^2);
            end
            plot(x, i_longitudinal(1, :));
            grid on;
            xlabel('Distance');
            ylabel('Membrane Current');
            title('Longitudinal Current vs Distance');
            should_plot_voltage2 = 0;
        elseif (should_plot_voltage == 1)
            %plot membrane voltage
            figure;
            subplot(1, 3, 1);
            plot(x, V(i, :));
            grid on;
            xlabel('Distance');
            ylabel('Membrane Voltage');
            title('Membrane Voltage vs Distance');
            %plot membrane current
            subplot(1, 3, 2);
            for j = 2:(length(x) - 1)
                im(j) = -gNa * m(i, j)^2 * h(i, j) * (Vm - Ena) - gL * (Vm - El);
            end
            plot(x, im(1, :));
            grid on;
            xlabel('Distance');
            ylabel('Membrane Current');
            title('Membrane Current vs Distance');
            %plot longitudinal axonal current
            subplot(1, 3, 3);
            for j = 2:(length(x) - 1)
                i_longitudinal(j) = ((2 * a * dx) / (4 * Ri * l)) * ((V(i, j-1) - 2 * Vm + V(i, j+1)) / dx^2);
            end
            plot(x, i_longitudinal(1, :));
            grid on;
            xlabel('Distance');
            ylabel('Membrane Current');
            title('Longitudinal Current vs Distance');
            should_plot_voltage = 0;
        end
    end
end
% plot the voltage in both time and x... showing the AP.
figure
surf(x, [0, t], V), title('Surface Plot for Bidirectional AP'), xlabel('Distance along axon (cm)'), ylabel('Time (ms)'), zlabel('Voltage (mV)')
shading flat
rotate3d on


% Make a video of the AP.
% video = VideoWriter('MomodirectionalAP---Z=1mm.mp4');
% open(video);
% figure
% for i = 1:length(t)
%     plot(x, V(i, :));
%     axis([-4.5, 4.5, -25, 100]);
%     title('Monodirectional Action Potential'), xlabel('Distance along axon (cm)'), ylabel('Voltage (mV)');
%     mov(i) = getframe(gcf);
%     writeVideo(video, mov(i));
% end