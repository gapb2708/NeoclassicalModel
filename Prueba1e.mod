% Dynare 
% Neoclassical model

% Endogenous variables:
var y inv k A c w rk q r b_g tau l g phi;

% Exogenous variables:
varexo e_a e_tau e_g;

% Parameters:
parameters alpha, rho_a, beta, delta, sigma_a, mu, rho_tau, rho_g, sigma_tau, theta_g, sigma_g, etha, zeta, psi, tau_ee, g_ee, k_ee, w_ee, rk_ee;

alpha = 0.33;
rho_a = 0.95;
beta = 0.995;
delta = 0.02;
etha=0.3;
zeta=1;
sigma_a = 0.01;
mu = 5;
theta_g=0.05;
rho_g=0.90;
rho_tau=0.90;
sigma_tau=0.01;
sigma_g=0.01;

l_ee=0.33;
A_ee=1;
tau_ee=0.35;
b_g_ee=0;
k_ee = (((1/beta)-1+delta)/(A_ee*alpha*(l_ee^(1-alpha))))^(1/(alpha-1));
w_ee=(1-alpha)*A_ee*(k_ee^alpha)*(l_ee^(-alpha));
rk_ee=alpha*A_ee*(k_ee^(alpha-1))*(l_ee^(1-alpha));
g_ee=(w_ee*l_ee*tau_ee)-(b_g_ee/beta)+b_g_ee;
psi=(g_ee^etha)*w_ee*(1-tau_ee)/l_ee^zeta;  
y_ee = A_ee*(k_ee^(alpha))*(l_ee^(1-alpha));
inv_ee= delta*k_ee;
c_ee= y_ee-inv_ee-g_ee;
q_ee = 1;
phi_ee = delta;
r_ee= (1/beta)-1;

% The Model:
model;

%Homes

%1
exp(q)*(exp(g)^etha)/(exp(c)*(exp(g)^etha)-psi*((exp(l)^(1+zeta))/(1+zeta)))=beta*(exp(g(+1))^etha)/(exp(c(+1))*(exp(g(+1))^etha)-psi*((exp(l(+1))^(1+zeta))/(1+zeta)))*(exp(rk(+1))+exp(q(+1))*(1-delta)+exp(q(+1))*exp(phi(+1))-(exp(inv(+1))/exp(k)));
%2
exp(phi) = (exp(inv)/exp(k(-1))) - mu/2*((exp(inv)/exp(k(-1)))-delta)^2;
%3
exp(q)*(1-mu*((exp(inv)/exp(k(-1)))-delta)) = 1;
%4
exp(k) = (1-delta)*exp(k(-1))+exp(phi)*exp(k(-1));
%5
exp(w)=(psi*(exp(l)^zeta))/((exp(g)^etha)*(1-tau));
%6
(exp(g)^etha)/(exp(c)*(exp(g)^etha)-psi*((exp(l)^(1+zeta))/(1+zeta)))=(1+r)*beta*(exp(g(+1))^etha)/(exp(c(+1))*(exp(g(+1))^etha)-psi*((exp(l(+1))^(1+zeta))/(1+zeta)));
%Firms
%7
exp(y) = exp(A)*(exp(k(-1))^alpha)*(exp(l)^(1-alpha));
%8
exp(w) = (1-alpha)*exp(A)*(exp(k(-1))^(alpha))*(exp(l)^(-alpha));
%9
exp(rk) = alpha*exp(A)*(exp(k(-1))^(alpha-1))*(exp(l)^(1-alpha));
%10
A = rho_a*A(-1) + sigma_a*e_a;
%11
%Goverment
%12
exp(g)+(1+r(-1))*b_g(-1) = tau*exp(w)*exp(l)+b_g;
%13
tau = (1-rho_tau)*tau_ee+rho_tau*tau(-1)+sigma_tau*e_tau;
%14
g = (1-rho_g)*g_ee+rho_g*g(-1)-theta_g*b_g(-1)+sigma_g*e_g;
%Economy
%15 
exp(y) = exp(c) + exp(inv)+exp(g);

end;


% Declare initial Values:
initval;
k = log(k_ee);
y = log(y_ee);
A = 0;
c = log(c_ee);
l= log(l_ee);
inv = log(inv_ee);
b_g=b_g_ee;
g= g_ee;
tau=tau_ee;
w = log(w_ee);
rk = log(rk_ee);
q = log(q_ee);
phi=log(phi_ee);
r=r_ee;
end;

shocks;
var e_a = 1;
var e_tau=1;
var e_g=1;
end;

check;
resid;

stoch_simul(order = 1,irf=48,graph); 

