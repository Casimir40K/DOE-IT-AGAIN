function doe_visualize(design, factors)
%DOE_VISUALIZE Create visualization plots for DOE design
%
%   doe_visualize(design, factors)
%
%   Creates histograms and scatter plots for continuous factors,
%   and bar charts for categorical factors.

    T = design.runSheet;
    
    % Separate continuous and categorical factors
    cont = factors(arrayfun(@(f) strcmpi(f.type, "continuous"), factors));
    cat = factors(arrayfun(@(f) strcmpi(f.type, "categorical"), factors));
    
    %% Plot continuous factors
    if ~isempty(cont)
        contNames = string({cont.name});
        nCont = numel(contNames);
        X = zeros(height(T), nCont);
        
        for i = 1:nCont
            if ismember(contNames(i), T.Properties.VariableNames)
                X(:, i) = T.(contNames(i));
            end
        end
        
        % Histograms
        figure('Name', 'DOE: Continuous Factor Distributions', ...
               'Position', [100, 100, 1000, 600]);
        nRows = ceil(nCont / 3);
        nCols = min(nCont, 3);
        
        for i = 1:nCont
            subplot(nRows, nCols, i);
            histogram(X(:, i), 'FaceColor', [0.2 0.6 0.8], 'EdgeColor', 'k');
            title(sprintf('%s (%s)', contNames(i), cont(i).units), ...
                  'Interpreter', 'none');
            xlabel('Value');
            ylabel('Count');
            grid on;
        end
        
        sgtitle('Continuous Factor Distributions', 'FontSize', 14, 'FontWeight', 'bold');
        
        % Pairwise scatter plots
        if nCont >= 2
            figure('Name', 'DOE: Pairwise Scatter Plots', ...
                   'Position', [150, 150, 900, 900]);
            
            [~, ax] = plotmatrix(X);
            
            % Label axes
            for i = 1:nCont
                xlabel(ax(nCont, i), contNames(i), 'Interpreter', 'none');
                ylabel(ax(i, 1), contNames(i), 'Interpreter', 'none');
            end
            
            sgtitle('Pairwise Scatter (Continuous Factors)', ...
                    'FontSize', 14, 'FontWeight', 'bold');
        end
        
        % 3D scatter if 3+ factors
        if nCont >= 3
            figure('Name', 'DOE: 3D Design Space', ...
                   'Position', [200, 200, 800, 700]);
            scatter3(X(:, 1), X(:, 2), X(:, 3), 50, ...
                    'MarkerFaceColor', [0.2 0.6 0.8], ...
                    'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.6);
            xlabel(sprintf('%s (%s)', contNames(1), cont(1).units), 'Interpreter', 'none');
            ylabel(sprintf('%s (%s)', contNames(2), cont(2).units), 'Interpreter', 'none');
            zlabel(sprintf('%s (%s)', contNames(3), cont(3).units), 'Interpreter', 'none');
            title('3D Design Space (First 3 Factors)', 'FontSize', 14, 'FontWeight', 'bold');
            grid on;
            view(45, 30);
        end
    else
        disp('No continuous factors to visualize.');
    end
    
    %% Plot categorical factors
    if ~isempty(cat)
        catNames = string({cat.name});
        nCat = numel(catNames);
        
        figure('Name', 'DOE: Categorical Factor Distribution', ...
               'Position', [250, 250, 1000, 400]);
        nRows = ceil(nCat / 3);
        nCols = min(nCat, 3);
        
        for i = 1:nCat
            subplot(nRows, nCols, i);
            if ismember(catNames(i), T.Properties.VariableNames)
                catData = categorical(T.(catNames(i)));
                catCounts = countcats(catData);
                bar(catCounts, 'FaceColor', [0.8 0.4 0.2], 'EdgeColor', 'k');
                set(gca, 'XTickLabel', categories(catData));
                title(catNames(i), 'Interpreter', 'none');
                ylabel('Count');
                grid on;
                
                % Rotate labels if needed
                if max(strlength(categories(catData))) > 8
                    xtickangle(45);
                end
            end
        end
        
        sgtitle('Categorical Factor Distribution', 'FontSize', 14, 'FontWeight', 'bold');
    end
    
    %% Summary statistics
    fprintf('\n========================================\n');
    fprintf('DESIGN SUMMARY STATISTICS\n');
    fprintf('========================================\n');
    fprintf('Total runs: %d\n', height(T));
    fprintf('Blocks: %d unique values\n', numel(unique(T.Block)));
    fprintf('Replicates: %d unique values\n', numel(unique(T.Replicate)));
    
    if ~isempty(cont)
        fprintf('\nContinuous factors:\n');
        for i = 1:numel(contNames)
            if ismember(contNames(i), T.Properties.VariableNames)
                vals = T.(contNames(i));
                fprintf('  %s: min=%.3f, max=%.3f, mean=%.3f, std=%.3f\n', ...
                    contNames(i), min(vals), max(vals), mean(vals), std(vals));
            end
        end
    end
    
    if ~isempty(cat)
        fprintf('\nCategorical factors:\n');
        for i = 1:numel(catNames)
            if ismember(catNames(i), T.Properties.VariableNames)
                fprintf('  %s: %d levels\n', catNames(i), ...
                    numel(unique(T.(catNames(i)))));
            end
        end
    end
    fprintf('========================================\n\n');
end
