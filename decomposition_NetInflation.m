%% Housekeeping

clc
clear
set(0,'DefaultFigureWindowStyle','docked'); % I prefer docked figures

%% %%%%%%%%%%%%%%%
%% Benchmark model 
models = {'fast', 'slow'};
shocks = {'mu','z','mps','pm'}; 
variables = {'u_hat','omega_hat','mismatch_hat','lambda_hat','W_hat','pi','pm_hat','z_hat','r_hat','Q_hat', 'v_hat', 'teta_hat'};

T = table(zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), 'VariableNames', {'Model', 'A_u', 'A_m', 'A_omega', 'A_lambda', 'A_e'});

for l=1:size(models,2)
    if strcmp(models{l}, 'fast')

        dynare faccini_linear -DFastReallocation -DBenchmarkModel

        SavedInfo.endo_names  = M_.endo_names;
        SavedSims.BenchmarkModel.fast = oo_.endo_simul;
        SavedIRF.BenchmarkModel.fast = oo_.irfs;
        SavedFEVD.BenchmarkModel.fast = oo_.variance_decomposition;
    else
        dynare faccini_linear -DSlowReallocation -DBenchmarkModel

        SavedInfo.BenchmarkModel.endo_names  = M_.endo_names;
        SavedSims.BenchmarkModel.slow = oo_.endo_simul;
        SavedIRF.BenchmarkModel.slow = oo_.irfs;
        SavedFEVD.BenchmarkModel.slow = oo_.variance_decomposition;
    end

    %table of coefficients
    newRow = table(models(l), A_u, A_m, A_omega, A_lambda, A_e, 'VariableNames', {'Model', 'A_u', 'A_m', 'A_omega', 'A_lambda', 'A_e'});
    T = [T; newRow];

    for j=1:size(shocks,2)
        for i=1:size(variables,2)
            x = [variables{i},'_eps_',shocks{j}];
            eval([variables{i} ' = transpose(SavedIRF.BenchmarkModel.(models{l}).(x))']);
            clear eval
        end
        
        %pi = A_u*u_hat + A_m*mismatch_hat + A_omega * omega_hat + A_lambda * lambda_hat + A_e*(lambda_hat(+1)+W_hat(+1))+beta * pi(+1) + pm_hat - kappa*z_hat;
        %normalization
        max_r = max(abs(r_hat)); 
        x_ax = 1:size(pi(1:end-1)); % x-axis data
        Net_inflation = pi(1:end-1)*(100/max_r) - beta*pi(2:end)*(100/max_r) - A_e*W_hat(2:end)*(100/max_r);
        Unemployment = A_u*u_hat(1:end-1)*(100/max_r);
        Mismatch = A_m*mismatch_hat(1:end-1)*(100/max_r);
        Omega =  A_omega*omega_hat(1:end-1)*(100/max_r);
        lambda = A_lambda*lambda_hat(1:end-1)*(100/max_r) + A_e*lambda_hat(2:end)*(100/max_r);
        shock = (pm_hat(1:end-1)-kappa*z_hat(1:end-1))*(100/max_r); 

        y = [Net_inflation Unemployment Mismatch Omega lambda shock]; 
        % Create stacked area chart
        %figure; % Opens a new figure window
        bar(x_ax, y(:, 2:end), 'stacked');
        hold on;
        plot(x_ax, y(:,1), '-k', 'LineWidth', 2);
        hold off;
        % Adding labels and title for clarity
        xlabel('time');
        ylabel('Net Inflation');
        name = ['Decomposition under shock ',shocks{j},' in ',models{l},' Economy'];
        title(name);
        legend('unemployment', 'mismatch', 'thightness','lambda', 'shock', 'net pi');
        
        %save
        short_name = [shocks{j},' ',models{l}];
        fileName = [short_name, '.png']; % You can change the file extension to .jpg, .pdf, etc.
        fullPath = fullfile('Decomp_BenchmarkModel', fileName);
        print(gcf, fullPath, '-dpng', '-r300'); % Save as a PNG with 300 dpi resolution

        %Loss function chart
        loss_Q = 0.25*(Q_hat.^2);
        loss_pi = pi.^2;
        Loss_max = max(abs(loss_Q+loss_pi));
        y = [loss_Q*(100/Loss_max) loss_pi*(100/Loss_max)]; 
        x_ax = 1:size(loss_pi); % x-axis data
        % Create stacked area chart
        %figure; % Opens a new figure window
        bar(x_ax, y(:, 1:end), 'stacked');
        hold off;
        % Adding labels and title for clarity
        xlabel('time');
        ylabel('Loss');
        name = ['Loss under shock ',shocks{j},' in ',models{l},' Economy'];
        title(name);
        legend('loss from output', 'loss from inflation');

        %save
        short_name = [shocks{j},' ',models{l}];
        fileName = [short_name, '.png']; % You can change the file extension to .jpg, .pdf, etc.
        fullPath = fullfile('Loss_BenchmarkModel', fileName);
        print(gcf, fullPath, '-dpng', '-r300'); % Save as a PNG with 300 dpi resolution

    end
end

%% %%%%%%%%%%
%% Optimal TR



for l=1:size(models,2)
    if strcmp(models{l}, 'fast')

        dynare faccini_linear -DFastReallocation -DOptimalTR

        SavedInfo.OptimalTR.endo_names  = M_.endo_names;
        SavedSims.OptimalTR.fast = oo_.endo_simul;
        SavedIRF.OptimalTR.fast = oo_.irfs;
        SavedFEVD.OptimalTR.fast = oo_.variance_decomposition;
    else
        dynare faccini_linear -DSlowReallocation -DOptimalTR

        SavedInfo.OptimalTR.endo_names  = M_.endo_names;
        SavedSims.OptimalTR.slow = oo_.endo_simul;
        SavedIRF.OptimalTR.slow = oo_.irfs;
        SavedFEVD.OptimalTR.slow = oo_.variance_decomposition;
    end

    for j=1:size(shocks,2)
        for i=1:size(variables,2)
            x = [variables{i},'_eps_',shocks{j}];
            eval([variables{i} ' = transpose(SavedIRF.OptimalTR.(models{l}).(x))']);
            clear eval
        end
        
        %pi = A_u*u_hat + A_m*mismatch_hat + A_omega * omega_hat + A_lambda * lambda_hat + A_e*(lambda_hat(+1)+W_hat(+1))+beta * pi(+1) + pm_hat - kappa*z_hat;

        max_r = max(abs(r_hat)); 
        x_ax = 1:size(pi(1:end-1)); % x-axis data
        Net_inflation = pi(1:end-1)*(100/max_r) - beta*pi(2:end)*(100/max_r) - A_e*W_hat(2:end)*(100/max_r);
        Unemployment = A_u*u_hat(1:end-1)*(100/max_r);
        Mismatch = A_m*mismatch_hat(1:end-1)*(100/max_r);
        Omega =  A_omega*omega_hat(1:end-1)*(100/max_r);
        lambda = A_lambda*lambda_hat(1:end-1)*(100/max_r) + A_e*lambda_hat(2:end)*(100/max_r);
        shock = (pm_hat(1:end-1)-kappa*z_hat(1:end-1))*(100/max_r); 

        y = [Net_inflation Unemployment Mismatch Omega lambda shock]; 
        % Create stacked area chart
        %figure; % Opens a new figure window
        bar(x_ax, y(:, 2:end), 'stacked');
        hold on;
        plot(x_ax, y(:,1), '-k', 'LineWidth', 2);
        hold off;
        % Adding labels and title for clarity
        xlabel('time');
        ylabel('Net Inflation');
        name = ['Decomposition under shock ',shocks{j},' in ',models{l},' Economy'];
        title(name);
        legend('unemployment', 'mismatch', 'thightness','lambda', 'shock', 'net pi');
        
        %save
        short_name = [shocks{j},' ',models{l}];
        fileName = [short_name, '.png']; % You can change the file extension to .jpg, .pdf, etc.
        fullPath = fullfile('Decomp_OptimalTR', fileName);
        print(gcf, fullPath, '-dpng', '-r300'); % Save as a PNG with 300 dpi resolution

        %Loss function chart
        loss_Q = 0.25*(Q_hat.^2);
        loss_pi = pi.^2;
        Loss_max = max(abs(loss_Q+loss_pi));
        y = [loss_Q*(100/Loss_max) loss_pi*(100/Loss_max)]; 
        x_ax = 1:size(loss_pi); % x-axis data
        % Create stacked area chart
        %figure; % Opens a new figure window
        bar(x_ax, y(:, 1:end), 'stacked');
        hold off;
        % Adding labels and title for clarity
        xlabel('time');
        ylabel('Loss');
        name = ['Loss under shock ',shocks{j},' in ',models{l},' Economy'];
        title(name);
        legend('loss from output', 'loss from inflation');

        %save
        short_name = [shocks{j},' ',models{l}];
        fileName = [short_name, '.png']; % You can change the file extension to .jpg, .pdf, etc.
        fullPath = fullfile('Loss_OptimalTR', fileName);
        print(gcf, fullPath, '-dpng', '-r300'); % Save as a PNG with 300 dpi resolution

    end
end


%% %%%%%%%%%%%%
%% Ramsey Model


for l=1:size(models,2)
    if strcmp(models{l}, 'fast')

        dynare faccini_linear -DFastReallocation -DRamseyModel

        SavedInfo.RamseyModel.endo_names  = M_.endo_names;
        SavedSims.RamseyModel.fast = oo_.endo_simul;
        SavedIRF.RamseyModel.fast = oo_.irfs;
        SavedFEVD.RamseyModel.fast = oo_.variance_decomposition;
    else
        dynare faccini_linear -DSlowReallocation -DRamseyModel

        SavedInfo.RamseyModel.endo_names  = M_.endo_names;
        SavedSims.RamseyModel.slow = oo_.endo_simul;
        SavedIRF.RamseyModel.slow = oo_.irfs;
        SavedFEVD.RamseyModel.slow = oo_.variance_decomposition;
    end

    for j=1:size(shocks,2)
        for i=1:size(variables,2)
            x = [variables{i},'_eps_',shocks{j}];
            eval([variables{i} ' = transpose(SavedIRF.RamseyModel.(models{l}).(x))']);
            clear eval
        end
        
        %pi = A_u*u_hat + A_m*mismatch_hat + A_omega * omega_hat + A_lambda * lambda_hat + A_e*(lambda_hat(+1)+W_hat(+1))+beta * pi(+1) + pm_hat - kappa*z_hat;

        max_r = max(abs(r_hat)); 
        x_ax = 1:size(pi(1:end-1)); % x-axis data
        Net_inflation = pi(1:end-1)*(100/max_r) - beta*pi(2:end)*(100/max_r) - A_e*W_hat(2:end)*(100/max_r);
        Unemployment = A_u*u_hat(1:end-1)*(100/max_r);
        Mismatch = A_m*mismatch_hat(1:end-1)*(100/max_r);
        Omega =  A_omega*omega_hat(1:end-1)*(100/max_r);
        lambda = A_lambda*lambda_hat(1:end-1)*(100/max_r) + A_e*lambda_hat(2:end)*(100/max_r);
        shock = (pm_hat(1:end-1)-kappa*z_hat(1:end-1))*(100/max_r); 

        y = [Net_inflation Unemployment Mismatch Omega lambda shock]; 
        % Create stacked area chart
        %figure; % Opens a new figure window
        bar(x_ax, y(:, 2:end), 'stacked');
        hold on;
        plot(x_ax, y(:,1), '-k', 'LineWidth', 2);
        hold off;
        % Adding labels and title for clarity
        xlabel('time');
        ylabel('Net Inflation');
        name = ['Decomposition under shock ',shocks{j},' in ',models{l},' Economy'];
        title(name);
        legend('unemployment', 'mismatch', 'thightness','lambda', 'shock', 'net pi');
        
        %save
        short_name = [shocks{j},' ',models{l}];
        fileName = [short_name, '.png']; % You can change the file extension to .jpg, .pdf, etc.
        fullPath = fullfile('Decomp_RamseyModel', fileName);
        print(gcf, fullPath, '-dpng', '-r300'); % Save as a PNG with 300 dpi resolution

        %Loss function chart
        loss_Q = 0.25*(Q_hat.^2);
        loss_pi = pi.^2;
        Loss_max = max(abs(loss_Q+loss_pi));
        y = [loss_Q*(100/Loss_max) loss_pi*(100/Loss_max)]; 
        x_ax = 1:size(loss_pi); % x-axis data
        % Create stacked area chart
        %figure; % Opens a new figure window
        bar(x_ax, y(:, 1:end), 'stacked');
        hold off;
        % Adding labels and title for clarity
        xlabel('time');
        ylabel('Loss');
        name = ['Loss under shock ',shocks{j},' in ',models{l},' Economy'];
        title(name);
        legend('loss from output', 'loss from inflation');

        %save
        short_name = [shocks{j},' ',models{l}];
        fileName = [short_name, '.png']; % You can change the file extension to .jpg, .pdf, etc.
        fullPath = fullfile('Loss_RamseyModel', fileName);
        print(gcf, fullPath, '-dpng', '-r300'); % Save as a PNG with 300 dpi resolution

    end
end

