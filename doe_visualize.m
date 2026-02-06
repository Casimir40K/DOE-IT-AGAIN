function doe_visualize(design, factors)
%doe_visualize Quick plots for continuous factors (actual values).

T = design.runSheet;

cont = factors(arrayfun(@(f) strcmpi(f.type,"continuous"), factors));
if isempty(cont)
    disp("No continuous factors to visualize.");
    return;
end

contNames = string({cont.name});
X = zeros(height(T), numel(contNames));
for i = 1:numel(contNames)
    X(:,i) = T.(contNames(i));
end

figure("Name","DOE: Histograms");
for i = 1:numel(contNames)
    subplot(ceil(numel(contNames)/2), 2, i);
    histogram(X(:,i));
    title(contNames(i));
    grid on;
end

if size(X,2) >= 2
    figure("Name","DOE: Pairwise Scatter");
    plotmatrix(X);
    sgtitle("Pairwise scatter (continuous actual)");
end
end
