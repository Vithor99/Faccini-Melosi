function v = vacancy_finder(delta, s, phi_0, psi, beta, k, xi_g, xi_b, y_g, y_b, b, W, eta, xi, sigma)
options = optimoptions('fsolve','MaxIterations',100000, 'MaxFunctionEvaluations', 100000, 'FunctionTolerance',1.0000e-10, 'OptimalityTolerance', 1.0000e-10);
v0 = [0.3;0.3];
Fun = fsolve(@(x) [ ((x(2)/(x(2)+s*(1-x(2)))) * ( W * (xi_g * y_g + xi_b * y_b) - ((b*((((eta-1)/eta)^(xi-1))*(y_g^xi * ((1 - (x(2) + (((phi_0 * ((x(1)/(x(2)+s*(1-x(2))))^psi))*xi_b*x(2))/(s*(phi_0 * ((x(1)/(x(2)+s*(1-x(2))))^psi))*xi_g - delta))))/(1-delta)) + y_b^xi * ((((phi_0 * ((x(1)/(x(2)+s*(1-x(2))))^psi))*xi_b*x(2))/(s*(phi_0 * ((x(1)/(x(2)+s*(1-x(2))))^psi))*xi_g - delta))/(1-delta))))^(1/(1-sigma*(1-xi))))/(1-beta*(1-delta))) ) + (s/(x(2)+s*(1-x(2)))) * xi_g * (((phi_0 * ((x(1)/(x(2)+s*(1-x(2))))^psi))*xi_b*x(2))/(s*(phi_0 * ((x(1)/(x(2)+s*(1-x(2))))^psi))*xi_g - delta)) * W * (y_g - y_b) - (k/(((phi_0 * ((x(1)/(x(2)+s*(1-x(2))))^psi)))/(x(1)/(x(2)+s*(1-x(2)))))));
((1-delta)*x(2)*(1 - (phi_0 * ((x(1)/(x(2)+s*(1-x(2))))^psi)))+delta-x(2))], v0, options);
v = Fun(1);


