%% Housekeeping

clc
clear
set(0,'DefaultFigureWindowStyle','docked'); % I prefer docked figures

%% Additinal analysis 
dynare faccini_linear -DFastReallocation -DBenchmarkModel

SavedInfo.endo_names  = M_.endo_names;
SavedSims.BenchmarkModel.fast = oo_.endo_simul;
SavedIRF.BenchmarkModel.fast = oo_.irfs;
SavedFEVD.BenchmarkModel.fast = oo_.variance_decomposition;

dynare faccini_linear -DSlowReallocation -DBenchmarkModel

SavedInfo.BenchmarkModel.endo_names  = M_.endo_names;
SavedSims.BenchmarkModel.slow = oo_.endo_simul;
SavedIRF.BenchmarkModel.slow = oo_.irfs;
SavedFEVD.BenchmarkModel.slow = oo_.variance_decomposition;


%% fast economy - Benchmark

for i=1:size(SavedInfo.endo_names,1)
    eval([SavedInfo.endo_names{i} ' = transpose(SavedSims.BenchmarkModel.fast(find(strcmp(SavedInfo.endo_names, SavedInfo.endo_names{i})),1000:10000))']);
    clear eval
end

figure; 
scatter(mismatch_hat, teta_hat, 'filled'); 
title('Fast economy');  
xlabel('mismatch values');  
ylabel('teta values');  
grid on;  

figure; 
scatter(mismatch_hat, u_hat, 'filled'); 
title('Fast economy');  
xlabel('mismatch values');  
ylabel('U values');  
grid on; 

%% slow economy - benchmark

for i=1:size(SavedInfo.endo_names,1)
    eval([SavedInfo.endo_names{i} ' = transpose(SavedSims.BenchmarkModel.slow(find(strcmp(SavedInfo.endo_names, SavedInfo.endo_names{i})),1000:10000))']);
    clear eval
end

figure; 
scatter(mismatch_hat, teta_hat, 'filled'); 
title('Slow economy');  
xlabel('mismatch values');  
ylabel('teta values');  
grid on;  

figure; 
scatter(mismatch_hat, u_hat, 'filled'); 
title('Slow economy');  
xlabel('mismatch values');  
ylabel('U values');  
grid on; 

%%
%% decomposition the interest rate path
dynare faccini_linear -DFastReallocation -DBenchmarkModel

SavedInfo.endo_names  = M_.endo_names;
SavedSims.BenchmarkModel.fast = oo_.endo_simul;
SavedIRF.BenchmarkModel.fast = oo_.irfs;
SavedFEVD.BenchmarkModel.fast = oo_.variance_decomposition;

variables = {'r_hat','pi','Q_hat','mps_hat'};

 for i=1:size(variables,2)
     x = [variables{i},'_eps_mps'];
     eval([variables{i} ' = transpose(SavedIRF.BenchmarkModel.fast.(x))']);
     clear eval
 end


 x_ax = 1:size(pi(1:end)); % x-axis data
 y = [r_hat phi_pi*pi phi_y*Q_hat mps_hat];
 % Create stacked area chart
 %figure; % Opens a new figure window
 bar(x_ax, y(:, 2:end), 'stacked');
 hold on;
 plot(x_ax, y(:,1), '-k', 'LineWidth', 2);
 hold off;
 % Adding labels and title for clarity
 xlabel('time');
 ylabel('R');
 legend('pi', 'Q', 'mps','R');
