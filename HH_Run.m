close all; clear all; clc;

T = 6.3; 
gNa = 120; 
gK = 36; 
Na_Out = 490;
Na_in = 50;
K_out = 20;
K_in = 400;
R = 8.314; 
F = 9.6485*10^4; 
E_Na = 1000*(R*T/F)*log(Na_Out/Na_in); 
E_K = 1000*(R*T/F)*log(K_out/K_in); 

V_rest = -60.34; 
m_0 = 0.05; 
n_0 = 0.32;
h_0 = 0.6; 
y_0 = [V_rest; m_0; n_0; h_0];

dt = [0,20];
options = odeset('RelTol',1e-4,'AbsTol',[1e-8,1e-8,1e-8,1e-8],'MaxStep',0.01);

[t,y] = ode45(@(t,y) HH(I),dt,y_0,options);
V = y(:,1); 
M = y(:,2); 
N = y(:,3);
H = y(:,4);

figure;
plot(t,V); title('Voltage'); xlabel('Time (ms)'); ylabel('Voltage (mV)');
figure
plot(t,M); 
hold on; plot(t,H);
hold on; plot(t,N); 
title('M, N, H Gates'); 
xlabel('Time (ms)'); 
ylabel('Probability Open');
legend('M','N','H');

I_K = gK.*N.^4.*(V-E_K);
I_Na = gK.*M.^3.*H.*(V-E_Na);

figure;
plot(t,I_K);
hold on; plot(t,I_Na); 
title('Ionic Currents'); xlabel('Time (ms)'); ylabel('Current (mA)');
legend('I_{Na}','I_{K}');