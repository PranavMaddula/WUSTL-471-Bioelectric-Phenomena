function Dydt = HH(t,y)
I_s = 0
V_rest = -65 ;
T = 6.3; 
T_k = 273.15+T; 
V = y(1,1); 
m = y(2,1); 
n = y(3,1); 
h = y(4,1); 
v_m = V-V_rest; 
gNa = 120; 
gK = 36; 
gL = 0.3; 
Na_out = 490;
Na_in = 50; 
K_out = 20; 
K_in = 400;
R = 8.314; 
F = 9.6485*10^4; 
E_Na = 1000*(R*T_k/F)*log(Na_out/Na_in); 
E_K = 1000*(R*T_k/F)*log(K_out/K_in);
E_L = -50; 
C = 1.0; 
k = 3^(0.1*T-0.63); 

alpha_m = .1*(25-v_m)./(exp((25-v_m)/10)-1); 
alpha_h = .07*exp(-v_m/20); 
alpha_n = .01*(10-v_m)./(exp((10-v_m)/10)-1); 
beta_m = 4*exp(-v_m/18);
beta_h = 1/(exp((30-v_m)/10)+1); 
beta_n = .125*exp(-v_m/80); 

dVdt = (1/C).*(I_s-gNa*m^3*h*(V-E_Na)-gK*n^4*(V-E_K)-gL*(V-E_L));
dmdt = (-(alpha_m+beta_m)*m+alpha_m)*k;
dndt = (-(alpha_n+beta_n)*n+alpha_n)*k;
dhdt = (-(alpha_h+beta_h)*h+alpha_h)*k;

Dydt = zeros(4,1);
Dydt(:,1) = [dVdt; dmdt; dndt; dhdt];
end