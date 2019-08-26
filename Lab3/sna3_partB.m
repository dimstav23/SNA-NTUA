%%Meros B
addpath(genpath('./epidemics_code'));

%%SIR
beta = 10^(-3);
gamma = [10^(-6) 10^(-5) 10^(-4) 10^(-3) 10^(-2) 10^(-1)];
S_0 = 570;
I_0 = 17;
R_0 = 0;
global i;
for i=1:6
	SIR(beta,gamma(i),S_0,I_0,R_0);
end

