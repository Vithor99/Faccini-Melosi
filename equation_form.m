clear all; 

% Define the equations for each variable
v = 'x(1)';
u = 'x(2)'; 
teta = '(v/(u+s*(1-u)))';  
phi_cap = '(phi_0 * (teta^psi))';  
omega = '((phi_cap)/teta)'; 
l_b = '((phi_cap*xi_b*u)/(s*phi_cap*xi_g - delta))';
l_g = '(1 - (u + l_b))'; 
varphi = '((eta-1)/eta)';
Q = '((varphi^(xi-1))*(y_g^xi * (l_g/(1-delta)) + y_b^xi * (l_b/(1-delta))))^(1/(1-sigma*(1-xi)))';
C = 'Q';  

% Long equation with variables to be replaced
FE = '((u/(u+s*(1-u))) * ( W * (xi_g * y_g + xi_b * y_b) - ((b*C)/(1-beta*(1-delta))) ) + (s/(u+s*(1-u))) * xi_g * l_b * W * (y_g - y_b) - (c_p/omega));'; %free entry condition 
U_dyn = '((1-delta)*u*(1 - phi_cap)+delta-u)';

% Replace variables with their definitions
FE= strrep(FE, 'C', C);
FE= strrep(FE, 'Q', Q);
FE= strrep(FE, 'varphi', varphi);
FE= strrep(FE, 'l_g', l_g);
FE= strrep(FE, 'l_b', l_b);
FE= strrep(FE, 'omega', omega);
FE= strrep(FE, 'phi_cap', phi_cap);
FE = strrep(FE, 'teta', teta);
FE = strrep(FE, 'u', u);
FE = strrep(FE, 'v', v);

U_dyn = strrep(U_dyn,'phi_cap', phi_cap);
U_dyn = strrep(U_dyn, 'teta', teta);
U_dyn = strrep(U_dyn, 'u', u);
U_dyn = strrep(U_dyn, 'v', v);

%Print the modified long equation
disp('Free entry condition:');
disp(FE);
disp('Unemployment dynamics:');
disp(U_dyn)





