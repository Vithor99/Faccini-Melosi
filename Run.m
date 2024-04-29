%% Housekeeping

clc
clear
set(0,'DefaultFigureWindowStyle','docked'); % I prefer docked figures...


%% Run models saving results

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

%% Optimized Taylor rule
dynare faccini_linear -DFastReallocation -DOptimalTR

SavedInfo.OptimalTR.endo_names  = M_.endo_names;
SavedSims.OptimalTR.fast = oo_.endo_simul;
SavedIRF.OptimalTR.fast = oo_.irfs;
SavedFEVD.OptimalTR.fast = oo_.variance_decomposition;

dynare faccini_linear -DSlowReallocation -DOptimalTR

SavedInfo.OptimalTR.endo_names  = M_.endo_names;
SavedSims.OptimalTR.slow = oo_.endo_simul;
SavedIRF.OptimalTR.slow = oo_.irfs;
SavedFEVD.OptimalTR.slow = oo_.variance_decomposition;

%% Ramsey
dynare faccini_linear -DFastReallocation -DRamseyModel

SavedInfo.RamseyModel.endo_names  = M_.endo_names;
SavedSims.RamseyModel.fast = oo_.endo_simul;
SavedIRF.RamseyModel.fast = oo_.irfs;
SavedFEVD.RamseyModel.fast = oo_.variance_decomposition;

dynare faccini_linear -DSlowReallocation -DRamseyModel

SavedInfo.RamseyModel.endo_names  = M_.endo_names;
SavedSims.RamseyModel.slow = oo_.endo_simul;
SavedIRF.RamseyModel.slow = oo_.irfs;
SavedFEVD.RamseyModel.slow = oo_.variance_decomposition;

%% Plotting IRFs

%% plot IRFs of fast and slow economy togther and saves them in IRFs folder 
models = {'BenchmarkModel','OptimalTR','RamseyModel'};

for j=1:size(models,2)
    if isfield(SavedIRF, models{j})
        fieldNames=fieldnames(SavedIRF.(models{j}).fast);
        for i =1:numel(fieldNames)
            fast = SavedIRF.(models{j}).fast.(fieldNames{i})'; 
            slow = SavedIRF.(models{j}).slow.(fieldNames{i})'; 
            
            %apply transformation
            index = strfind(fieldNames{i}, '_');
            r_hat = ['r_hat_eps_',fieldNames{i}(index(end) + 1 : end)]; 
            r_max_fast =  max(abs(SavedIRF.(models{j}).fast.(r_hat)));
            r_max_slow =  max(abs(SavedIRF.(models{j}).slow.(r_hat)));

            fast = (fast/r_max_fast)*100;
            slow = (slow/r_max_slow)*100;
            A = [fast slow];

            name = strrep(fieldNames{i}, '_', ' ');
            %plot
            plot(A(:,1), '-or', 'LineWidth', 1.5, 'MarkerSize', 6);
            hold on;
            plot(A(:,2), '-ob', 'LineWidth', 1.5, 'MarkerSize', 6);
            hold off;
            xlabel('Time Steps');
            ylabel('IRF');
            title(name);
            legend('Fast', 'Slow');

            %save
            folderName = ['IRFs_',models{j}];
            fileName = [name, '.png']; % You can change the file extension to .jpg, .pdf, etc.
            fullPath = fullfile(folderName, fileName);
            print(gcf, fullPath, '-dpng', '-r300'); % Save as a PNG with 300 dpi resolution
        end
    else
        fprintf('Field %s does not exist.\n', models{j});
    end
end
%% plot IRFs of benchmark, OptimalTR and Ramsey together 

models = {'fast','slow'};

for j=1:size(models,2)
    fieldNames=fieldnames(SavedIRF.BenchmarkModel.(models{j}));
    for i =1:numel(fieldNames)
        Benchmark = SavedIRF.BenchmarkModel.(models{j}).(fieldNames{i})';
        OptimalTR = SavedIRF.OptimalTR.(models{j}).(fieldNames{i})';
        Ramsey = SavedIRF.RamseyModel.(models{j}).(fieldNames{i})'; 

        index = strfind(fieldNames{i}, '_');
        r_hat = ['r_hat_eps_',fieldNames{i}(index(end) + 1 : end)];
        r_max_Benchmark =  max(abs(SavedIRF.BenchmarkModel.(models{j}).(r_hat)));
        r_max_OptimalTR =  max(abs(SavedIRF.OptimalTR.(models{j}).(r_hat)));
        r_max_Ramsey =  max(abs(SavedIRF.RamseyModel.(models{j}).(r_hat)));

        Benchmark = (Benchmark/r_max_Benchmark)*100;
        OptimalTR = (OptimalTR/r_max_OptimalTR)*100;
        Ramsey = (Ramsey/r_max_Ramsey)*100;
        A = [Benchmark OptimalTR Ramsey];

        name = strrep(fieldNames{i}, '_', ' ');
        %plot
        plot(A(:,1), '-or', 'LineWidth', 1.5, 'MarkerSize', 6);
        hold on;
        plot(A(:,2), '-ob', 'LineWidth', 1.5, 'MarkerSize', 6);
        hold on;
        plot(A(:,3), '-og', 'LineWidth', 1.5, 'MarkerSize', 6);
        hold off;
        xlabel('Time Steps');
        ylabel('IRF');
        title(name);
        legend('Benchmark Model','Optimal TR','Ramsey Model');

        %save
        folderName = ['IRFs_',models{j}];
        fileName = [name, '.png']; % You can change the file extension to .jpg, .pdf, etc.
        fullPath = fullfile(folderName, fileName);
        print(gcf, fullPath, '-dpng', '-r300'); % Save as a PNG with 300 dpi resolution
    end
end












