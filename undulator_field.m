
me = 9.11 * 10^-31;
qe = 1.602 * 10^-19;
c = 2.998*10^8;
lambda_u = 0.05
lambda_l = 790 * 10^-9;
gamma = 538.26; 

B_u = 2*pi*me*c / (qe * lambda_u) * sqrt(2*(lambda_l/lambda_u * 2*gamma^2 - 1))

%MeV to gamma

m0c2 = 0.511; %in MeV
Ek = 1:260;
gamma = 1 + Ek/m0c2;
% plot(Ek,gamma,'.')

%compute lambda_laser from the resonance condition
B_u = 0.86514;
c = 2.998*10^8;
me = 9.11 * 10^-31;
lambda_u = 0.05;
gamma = 538.26; 
K = e*B_u*lambda_u/(2*pi*me*c);
lambda_laser = lambda_u/(2*gamma^2)*(1+0.5*K^2)

