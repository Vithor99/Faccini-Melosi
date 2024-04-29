close all;
warning off

%%%%%%%%%%% Faccini & Melosi 

var  mu_hat %preference for consumption
z_hat %technology shock
mps_hat %monetary policy shock
pm_hat %price mark-up shock

teta_hat % thightness
u_hat %beginning of period unemployment  
phi_cap_hat %job finding rate 
omega_hat %vacancy filling rate 
s_g_hat % Firm's surplus from a good match 
s_b_hat % Firm's surplus from a bad match 
W_hat % PDV of (real) marginal revenues for service firms (= PDV of real marginal cost for price setters)
v_hat %Vacancies 
l_b_hat % matches of quality bad (beginning of period)
l_g_hat % matches of quality good (beginning of period)
C_hat % consumption
lambda_hat %MU of consumption
pi %log of inflation
varphi_hat % (spot) real marginal revenue for service firms (= spot real mc of price setters) 
r_hat % nominal interest rate 
Q_hat %Output

%Auxiliary variables 
EE_hat
UE_hat
AC_hat
mismatch_hat
;

varexo eps_mu %innovation on preference for consumption
eps_z %innovation on tech shock
eps_mps %innovation on monetary shock
eps_pm %innovation on price mark-up shock 
; 

parameters

@#ifdef BenchmarkModel
    phi_pi       %monetary rule coeff.1 
    phi_y       %monetary rule coeff.2
@#endif

@#ifdef OptimalTR 
    phi_pi      %monetary rule coeff.1 
    phi_y       %monetary rule coeff.2
@#endif

rho_mu % AR coefficient on preference shock
rho_z % AR coefficient on tech shock
rho_mps %AR coeff on monetary shock
rho_pm %AR coeff on price mark-up shock
rho_r %AR on Taylor rule

s %search intensity
delta % probability of job destruction
phi_0 %matching technology 
psi %coeff on thightness for match creations 
beta %discount factor  
k % (keriodic) cost per perod from open vacancy 
xi_g % cond. on matching, prob that match is good type 
xi_b % cond. on matching, prob that match is bad type (1-xi_g)
y_g %productivity of good match
y_b %productivity of bad match 
b %Unemployment benefit 
zeta %Calvo parameter
eta %elasticity of substitution for consumption
xi %work effort elasticity
sigma %intertemporal elasticity of consumption

%steady state variables 
mu_ss
z_ss
mps_ss
pm_ss
r_ss 
varphi_ss 
W_ss 
v_ss 
u_ss  
teta_ss 
phi_cap_ss 
omega_ss
l_b_ss 
l_g_ss
Q_ss 
C_ss
s_g_ss 
s_b_ss 
lambda_ss
%auxilary variables
EE_ss
UE_ss 
AC_ss 
mismatch_ss 

%other auxiliary parameters
kappa
A_u 
A_sg 
A_sb 
A_lb 
A_omega  

A_w 
A_lambda 
A_m
A_e
;

%%%%%%%%%%% Parametrization 
%monetary policy params
@#ifdef BenchmarkModel
    phi_pi = 1.8;%monetary rule coeff.1 
    phi_y = 0.25; %monetary rule coeff.2
@#endif

%others
s = 0.5; %search intensity 
rho_mu = 0.8; % AR coefficient on preference shock
rho_z = 0.8 ; %AR coeff on tech shock 0.8
rho_mps = 0.8; %AR coeff on monetar shock
rho_pm = 0.8; %AR coeff on price mark-up coeff
delta = 0.02;% probability of job destruction
phi_0 = 0.3284; %matching technology 
psi = 0.5;%coeff on thightness for match creations 
beta = 0.9987; %discount factor 
k = 0.49; %0.5;%(periodic) cost per perod from open vacancy 
xi_g = 0.5; % cond. on matching, prob that match is good type 
xi_b = 0.5; %cond. on matching, prob that match is bad type (1-xi_g)
y_g = 1.08; %productivity of good match
y_b = 1; %productivity of bad match 
zeta = 0.925; %Calvo parameter (portion of firms that don't reset their prices)            
eta = 6; %elasticity of substitution for consumption                                        
b = 0.8082; %unemployment benefit / Utility of leisure 
xi = 1.5; %intertemporal work effort elasticity (from MPV) 1.37
sigma = 0.33; %intertemporal consumption elasticity 

%% steady state
mu_ss= 1;
z_ss=1;
mps_ss=1;
pm_ss = 1; 
r_ss = log(1/beta);
varphi_ss = (eta-1)/eta;
W_ss = varphi_ss/(1-beta*(1-delta));

v_ss = vacancy_finder(delta, s, phi_0, psi, beta, k, xi_g, xi_b, y_g, y_b, b, W_ss, eta, xi, sigma); 
u_ss = unemp_finder(delta, s, phi_0, psi, beta, k, xi_g, xi_b, y_g, y_b, b, W_ss, eta, xi, sigma); 

teta_ss = (v_ss/(u_ss+s*(1-u_ss)));  % Example equation for variable x
phi_cap_ss = (phi_0 * ((v_ss/(u_ss+s*(1-u_ss)))^psi));  % Example equation for variable y
omega_ss = phi_cap_ss / teta_ss; 
l_b_ss =((phi_cap_ss*xi_b*u_ss)/(s*phi_cap_ss*xi_g - delta));
l_g_ss = 1 -( u_ss + l_b_ss); 
Q_ss = ((varphi_ss^(xi-1))*(y_g^xi * (l_g_ss/(1-delta)) + y_b^xi * (l_b_ss/(1-delta))))^(1/(1-sigma*(1-xi)));
C_ss = Q_ss; 
s_g_ss = y_g * W_ss - (b*C_ss/(1-beta*(1-delta)));
s_b_ss = y_b * W_ss - (b*C_ss/(1-beta*(1-delta)));
lambda_ss = (C_ss)^(-sigma);

%auxilary variables
EE_ss= (s*phi_cap_ss*(l_b_ss*(xi_g)))/(l_g_ss+l_b_ss);
UE_ss = phi_cap_ss; 
AC_ss = EE_ss/UE_ss;
mismatch_ss = l_b_ss/(l_b_ss+l_g_ss);


%% param switch
@#ifdef FastReallocation
    s=0.8;
    xi_g=0.8;
    xi_b=1-xi_g;
    psi=0.8;
@#endif

@#ifdef SlowReallocation
    s=0.3;
    xi_g=0.3;
    xi_b=1-xi_g;
    psi=0.3;
@#endif

%other auxiliary parameters
kappa = ((1-zeta)*(1-beta*zeta)*(1/zeta));
A_sg = xi_g*(u_ss+s*l_b_ss);
A_sb = (u_ss*xi_b-s*xi_g*l_b_ss);
A_w = (A_sg*y_g+A_sb*y_b);
A_lb = -(xi_g*s*(s_g_ss-s_b_ss))/A_w;

A_u = (((-(u_ss*(xi_g*s_g_ss+xi_b*s_b_ss-(k/omega_ss)*(1-s)))/A_w)- A_lb* u_ss*mismatch_ss)/varphi_ss)*kappa; %weight on unemployment
A_omega = ((-((k/omega_ss)*(u_ss*(1-s)+s))/A_w)/varphi_ss)*kappa; %weight on thightness
A_lambda = (((-((A_sg+A_sb)*(b *(lambda_ss^-1) / (1- beta*(1-delta))))/A_w)+beta*(1-delta)*W_ss)/varphi_ss)*kappa; %weight on lambda
A_m = ((A_lb * mismatch_ss *(1-u_ss))/varphi_ss)*kappa; %weight on mismatch
A_e = ((- beta*(1-delta)*W_ss)/varphi_ss)*kappa; %weight on expectation

%%%%%%%%%%%%%% Model 

model;
%Taylor rules 
@#ifdef BenchmarkModel 
    r_hat = phi_pi * pi + phi_y * Q_hat + mps_hat;       %Taylor rule 
@#endif

@#ifdef OptimalTR 
    r_hat = phi_pi * pi + phi_y * Q_hat + mps_hat;       %Taylor rule 
@#endif


%log(mu) = rho_mu * log(mu(-1)) + eps_mu;
mu_hat = rho_mu * mu_hat(-1) + eps_mu;

%log(z) = rho_z * log(z(-1)) - eps_z;
z_hat = rho_z * z_hat(-1) - eps_z;

%log(mps) = rho_mps * log(mps(-1)) + eps_mps; 
mps_hat = rho_mps * mps_hat(-1) - eps_mps;

%log(pm) = rho_pm * log(pm(-1)) + eps_pm;
pm_hat = rho_pm * pm_hat(-1) + eps_pm;

%teta = v / (u + s*(1-u));                                                  %thightness
v_hat = (1/v_ss)*(teta_ss*teta_hat*(u_ss+s-s*u_ss)+u_ss*u_hat*teta_ss*(1-s)); 

%u = (1-delta)* u(-1) * (1- phi_cap(-1)) + delta;                           %unemployment dynamics
u_hat = (1-delta)*(u_hat(-1)*(1-phi_cap_ss)-phi_cap_ss*phi_cap_hat(-1));

%phi_cap = phi_0 * (teta)^psi;                                              %matching function
phi_cap_hat = psi*teta_hat;

%lambda = mu/C;                                                             %MU of consumption
lambda_hat = mu_hat - sigma*C_hat; 

%log(lambda) = log(lambda(+1)) - log(1/beta) + r - pi(+1);                  %log linear Euler 
lambda_hat(+1)-lambda_hat+r_hat-pi(+1)=0; 

pi = kappa*(varphi_hat - z_hat) + beta * pi(+1) + pm_hat;                   %NKPC in logs

%omega = phi_cap/teta;                                                      %vacancy filling rate
omega_hat = phi_cap_hat - teta_hat;

%(k/omega) = (u/(u+s*(1-u))) *(xi_b * s_b + xi_g * s_g) + ((s*(1-u))/(u+s*(1-u))) * (xi_g * (l_b/(1-u))*(s_g - s_b)); %free entry condition
u_hat*u_ss*(xi_g*s_g_ss+xi_b*s_b_ss-(k/omega_ss)*(1-s)) + s_g_hat*s_g_ss*xi_g*(u_ss+s*l_b_ss) + s_b_hat*s_b_ss*(u_ss*xi_b-s*xi_g*l_b_ss) + l_b_hat*l_b_ss*xi_g*s*(s_g_ss-s_b_ss) + omega_hat*(k/omega_ss)*(u_ss*(1-s)+s) = 0; %log linear original

%W_ss*W_hat*((y_g^xi)*xi_g*(u_ss+s*l_b_ss)+(y_b^xi)*(u_ss*xi_b-s*xi_g*l_b_ss)) = - lambda_hat*(b *(lambda_ss^-1)/(1- beta*(1-delta)))*(xi_g*(u_ss+s*l_b_ss)+(u_ss*xi_b-s*xi_g*l_b_ss)) - u_hat*u_ss*(xi_g*s_g_ss+xi_b*s_b_ss-(k/omega_ss)*(1-s)) - l_b_hat*l_b_ss*xi_g*s*(s_g_ss-s_b_ss) - omega_hat*(k/omega_ss)*(u_ss*(1-s)+s)

% new: s_b = (y_b^xi)*W - (b *(lambda^-1) / (1- beta*(1-delta)))            Surplus function from bad match 
s_b_ss*s_b_hat=(y_b^xi)*W_ss*W_hat+ (b *(lambda_ss^-1) / (1- beta*(1-delta)))*lambda_hat;

% new: s_g = (y_g^xi)*W - (b *(lambda^-1) / (1- beta*(1-delta)));          %Surplus function from good match 
s_g_ss*s_g_hat=(y_g^xi)*W_ss*W_hat+ (b *(lambda_ss^-1) / (1- beta*(1-delta)))*lambda_hat;

%W = (1/xi)*(varphi^xi)*(lambda^(xi-1))  + beta*(1-delta)*(lambda(+1)/lambda)*W(+1)                 %Stream of real Marginal revenues 
W_ss*W_hat = ((1/xi)*(varphi_ss^xi)*(lambda_ss^(xi-1)))*(xi*varphi_hat+ (xi-1)*lambda_hat) + beta*(1-delta)*W_ss*(lambda_hat(+1)-lambda_hat+W_hat(+1));

%l_b = (1-delta)*((1- s*phi_cap(-1) * xi_g)*l_b(-1) + phi_cap(-1) * xi_b * u(-1));      %low of motion for bad matches 
l_b_ss*l_b_hat=(1-delta)*l_b_ss*l_b_hat(-1)*(1-s*phi_cap_ss*xi_g)+(1-delta)*phi_cap_ss*phi_cap_hat(-1)*(xi_b*u_ss-l_b_ss*s*xi_g)+(1-delta)*phi_cap_ss*xi_b*u_ss*u_hat(-1);

%l_g = 1 -(u + l_b);
l_g_ss*l_g_hat = -u_ss*u_hat-l_b_ss*l_b_hat; 

%Q = z*((lambda*varphi)^(xi-1))*(y_g^xi * (l_g(+1)/(1-delta)) + y_b^xi * (l_b(+1)/(1-delta)));               %Production function
Q_ss*(Q_hat-z_hat)=((y_g^xi)/(1-delta))*((lambda_ss*varphi_ss)^(xi-1))*l_g_ss*(l_g_hat(+1)+(xi-1)*(lambda_hat+varphi_hat))+((y_b^xi)/(1-delta))*((lambda_ss*varphi_ss)^(xi-1))*l_b_ss*(l_b_hat(+1)+(xi-1)*(lambda_hat+varphi_hat));

%C = Q;                                                                     %Market clearing 
C_ss*C_hat = Q_ss*Q_hat;

%auxilary variables
UE_hat = phi_cap_hat;

%EE= (s*phi_cap*(l_b*(xi_g)))/(l_b+l_g);
EE_hat*(1-u_ss*EE_ss)=s*xi_g*phi_cap_ss*l_b_ss*(phi_cap_hat+l_b_hat)+u_ss*EE_ss*u_hat; 

AC_hat = EE_hat-UE_hat;

%mismatch = l_b/(l_b+l_g);
mismatch_ss*mismatch_hat = (l_b_ss/(1-u_ss))*l_b_hat+(u_ss*mismatch_ss/(1-u_ss))*u_hat; 
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
steady_state_model;
pi=0;
end;


shocks;    
var eps_mu; stderr 0.0055; 
var eps_z; stderr 0.003; 
var eps_mps; stderr 0.001;
var eps_pm; stderr 0.0001; 
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Optimal Rules and Ramsey

@#ifdef RamseyModel
    planner_objective pi^2 + 0.25*Q_hat^2;
    ramsey_model(instruments=(r_hat),planner_discount=beta);
@#endif

@#ifdef OptimalTR
    optim_weights; 
    pi 1;
    Q_hat 0.25;
    end;
    osr_params phi_pi phi_y;
    phi_pi= 1.5;
    phi_y = 0;
    osr(irf=0, opt_algo=2);
@#endif


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
steady;
resid;
check;

options_.nograph = 1;   % Do not plot IRFs immediately

stoch_simul(order = 1,irf=30,irf_plot_threshold=0) u_hat Q_hat pi EE_hat AC_hat mismatch_hat varphi_hat pm_hat z_hat r_hat omega_hat l_b_hat lambda_hat W_hat v_hat teta_hat;

stoch_simul(order = 1
    ,irf=0
    ,periods= 10000 
    ,irf_plot_threshold=0) u_hat Q_hat pi EE_hat AC_hat mismatch_hat varphi_hat pm_hat z_hat r_hat omega_hat l_b_hat lambda_hat W_hat v_hat teta_hat;


