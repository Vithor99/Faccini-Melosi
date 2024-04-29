%% FEVD check 
%% Housekeeping

clc
clear
set(0,'DefaultFigureWindowStyle','docked'); % I prefer docked figures...
%% Benchmark
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
%% slow
% Define row names
colNames = {'u_hat', 'Q_hat', 'pi', 'EE_hat', 'AC_hat', 'mismatch_hat', 'varphi_hat', 'pm_hat', 'z_hat', 'r_hat', 'omega_hat', 'l_b_hat', 'lambda_hat', 'W_hat','v_hat','teta_hat'};
rowNames = {'eps_mu', 'eps_z', 'eps_mps', 'eps_pm'};
T = array2table(nan(size(SavedFEVD.BenchmarkModel.slow,2), size(SavedFEVD.BenchmarkModel.slow,1)), 'VariableNames', colNames, 'RowNames', rowNames);
for i = 1:length(colNames)
    T.(colNames{i}) = SavedFEVD.BenchmarkModel.slow(i,:)';
end
T_focus_slow = T(:, [1, 2, 3, 10]); 
avg =mean(T_focus_slow{:,:},2);
T_focus_slow.Avg = avg;
disp(T_focus_slow)

%% fast
% Define row names
colNames = {'u_hat', 'Q_hat', 'pi', 'EE_hat', 'AC_hat', 'mismatch_hat', 'varphi_hat', 'pm_hat', 'z_hat', 'r_hat', 'omega_hat', 'l_b_hat', 'lambda_hat', 'W_hat','v_hat','teta_hat'};
rowNames = {'eps_mu', 'eps_z', 'eps_mps', 'eps_pm'};
T = array2table(nan(size(SavedFEVD.BenchmarkModel.fast,2), size(SavedFEVD.BenchmarkModel.fast,1)), 'VariableNames', colNames, 'RowNames', rowNames);
for i = 1:length(colNames)
    T.(colNames{i}) = SavedFEVD.BenchmarkModel.fast(i,:)';
end
T_focus_fast = T(:, [1, 2, 3, 10]); 
avg =mean(T_focus_fast{:,:},2);
T_focus_fast.Avg = avg;
disp(T_focus_fast)

%% slow and fast in comparison 
colNames = {'fast','slow'};
rowNames = {'eps_mu', 'eps_z', 'eps_mps', 'eps_pm'};
T = array2table(nan(size(SavedFEVD.BenchmarkModel.fast,2),2), 'VariableNames', colNames, 'RowNames', rowNames);
T.fast = T_focus_fast.Avg;
T.slow = T_focus_slow.Avg; 

disp(T)
