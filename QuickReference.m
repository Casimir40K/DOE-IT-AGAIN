%% ========================================================================
%  DOE TOOLKIT - QUICK REFERENCE GUIDE
%  ========================================================================
%  Common patterns and recipes for creating DOE designs
%% ========================================================================

%% PATTERN 1: Simple 2-Factor Full Factorial
% -------------------------------------------------------------------------
factors = [
    struct('name', 'Temp', 'type', 'continuous', ...
           'low', 20, 'high', 80, 'units', 'C', ...
           'levels', [], 'transform', 'none')
    struct('name', 'Time', 'type', 'continuous', ...
           'low', 10, 'high', 60, 'units', 'min', ...
           'levels', [], 'transform', 'none')
];

[design, report, ~] = createDOE(factors, 'type', 'fullfact');
% Result: 4 runs (2^2 factorial)

%% PATTERN 2: CCD for Response Surface (3 factors)
% -------------------------------------------------------------------------
factors = [
    struct('name', 'Temp', 'type', 'continuous', ...
           'low', 60, 'high', 90, 'units', 'C', ...
           'levels', [], 'transform', 'none')
    struct('name', 'pH', 'type', 'continuous', ...
           'low', 6, 'high', 9, 'units', '-', ...
           'levels', [], 'transform', 'none')
    struct('name', 'Conc', 'type', 'continuous', ...
           'low', 0.1, 'high', 1.0, 'units', 'M', ...
           'levels', [], 'transform', 'none')
];

[design, report, ~] = createDOE(factors, ...
    'type', 'ccd', ...
    'ccdAlpha', 'rotatable', ...
    'centerPoints', 6, ...
    'exportCSV', true);
% Result: ~20 runs (8 factorial + 6 axial + 6 center)

%% PATTERN 3: Mixed Continuous + Categorical
% -------------------------------------------------------------------------
factors = [
    struct('name', 'Catalyst', 'type', 'categorical', ...
           'levels', ["A", "B", "C"], ...
           'low', [], 'high', [], 'units', '-', 'transform', 'none')
    struct('name', 'Temp', 'type', 'continuous', ...
           'low', 25, 'high', 75, 'units', 'C', ...
           'levels', [], 'transform', 'none')
    struct('name', 'Time', 'type', 'continuous', ...
           'low', 30, 'high', 120, 'units', 'min', ...
           'levels', [], 'transform', 'none')
];

[design, report, ~] = createDOE(factors, ...
    'type', 'fullfact', ...
    'centerPoints', 3);
% Result: 4 continuous runs × 3 categorical = 12 base runs + centers

%% PATTERN 4: Screening Design (Many Factors)
% -------------------------------------------------------------------------
% Use fractional factorial for 5+ factors
factors = [
    struct('name', 'A', 'type', 'continuous', ...
           'low', 10, 'high', 50, 'units', '-', ...
           'levels', [], 'transform', 'none')
    struct('name', 'B', 'type', 'continuous', ...
           'low', 1, 'high', 5, 'units', '-', ...
           'levels', [], 'transform', 'none')
    struct('name', 'C', 'type', 'continuous', ...
           'low', 100, 'high', 200, 'units', '-', ...
           'levels', [], 'transform', 'none')
    struct('name', 'D', 'type', 'continuous', ...
           'low', 0.1, 'high', 1.0, 'units', '-', ...
           'levels', [], 'transform', 'none')
];

[design, report, ~] = createDOE(factors, ...
    'type', 'fracfact', ...
    'fracGenerators', 'a b c abc');  % 2^(4-1) = 8 runs
% Result: 8 runs instead of 16

%% PATTERN 5: With Constraints
% -------------------------------------------------------------------------
factors = [
    struct('name', 'Temp', 'type', 'continuous', ...
           'low', 60, 'high', 100, 'units', 'C', ...
           'levels', [], 'transform', 'none')
    struct('name', 'Pressure', 'type', 'continuous', ...
           'low', 1, 'high', 10, 'units', 'bar', ...
           'levels', [], 'transform', 'none')
];

% Define constraints
constraints = {};
constraints{1} = @(T) ~(T.Temp > 90 & T.Pressure > 8);  % Safety limit
constraints{2} = @(T) T.Temp >= 65 | T.Pressure >= 3;   % Minimum activation

[design, report, ~] = createDOE(factors, ...
    'type', 'fullfact', ...
    'constraints', constraints);

%% PATTERN 6: With Blocking and Replicates
% -------------------------------------------------------------------------
factors = [
    struct('name', 'Factor1', 'type', 'continuous', ...
           'low', 0, 'high', 100, 'units', '-', ...
           'levels', [], 'transform', 'none')
    struct('name', 'Factor2', 'type', 'continuous', ...
           'low', 0, 'high', 50, 'units', '-', ...
           'levels', [], 'transform', 'none')
];

[design, report, ~] = createDOE(factors, ...
    'type', 'ccd', ...
    'replicates', 2, ...        % Run everything twice
    'numBlocks', 3, ...         % Divide into 3 blocks
    'blockingMethod', 'sequential', ...
    'randomize', true);

%% PATTERN 7: Latin Hypercube for Computer Experiments
% -------------------------------------------------------------------------
factors = [
    struct('name', 'X1', 'type', 'continuous', ...
           'low', 0, 'high', 1, 'units', '-', ...
           'levels', [], 'transform', 'none')
    struct('name', 'X2', 'type', 'continuous', ...
           'low', 0, 'high', 1, 'units', '-', ...
           'levels', [], 'transform', 'none')
    struct('name', 'X3', 'type', 'continuous', ...
           'low', 0, 'high', 1, 'units', '-', ...
           'levels', [], 'transform', 'none')
];

[design, report, ~] = createDOE(factors, ...
    'type', 'lhs', ...
    'lhsSamples', 100, ...
    'lhsCriterion', 'maximin');
% Result: 100 well-distributed points

%% PATTERN 8: Box-Behnken (Efficient RSM)
% -------------------------------------------------------------------------
factors = [
    struct('name', 'A', 'type', 'continuous', ...
           'low', 10, 'high', 30, 'units', '-', ...
           'levels', [], 'transform', 'none')
    struct('name', 'B', 'type', 'continuous', ...
           'low', 5, 'high', 15, 'units', '-', ...
           'levels', [], 'transform', 'none')
    struct('name', 'C', 'type', 'continuous', ...
           'low', 100, 'high', 200, 'units', '-', ...
           'levels', [], 'transform', 'none')
];

[design, report, ~] = createDOE(factors, ...
    'type', 'boxbehnken', ...
    'centerPoints', 5);
% Result: 13 + 5 = 18 runs (vs 20+ for CCD)

%% TIPS FOR CHOOSING DESIGN TYPE
% -------------------------------------------------------------------------
% Full Factorial:      2-4 factors, need all interactions
% Fractional Fact:     5+ factors, screening studies
% CCD:                 Response surface, 2-5 factors, quadratic model
% Box-Behnken:         Response surface, 3+ factors, fewer runs than CCD
% LHS:                 Computer experiments, many factors, space-filling

%% COMMON PARAMETER COMBINATIONS
% -------------------------------------------------------------------------

% AGGRESSIVE SCREENING (minimize runs):
% type: 'fracfact'
% replicates: 1
% centerPoints: 0

% THOROUGH OPTIMIZATION (maximize information):
% type: 'ccd'
% ccdAlpha: 'rotatable'
% centerPoints: 6
% replicates: 2

% EXPLORATORY (balance runs vs info):
% type: 'boxbehnken'
% centerPoints: 3-5
% replicates: 1

%% LOADING SAVED DESIGNS
% -------------------------------------------------------------------------
% Load from MAT file:
load('./DOE_Output/ProjectName_bundle.mat');
disp(design.runSheet);

% Load from CSV:
runSheet = readtable('./DOE_Output/ProjectName_runsheet.csv');

%% ANALYZING RESULTS (Basic)
% -------------------------------------------------------------------------
% After running experiments, update Response columns:
% design.runSheet.Response1 = [measured_values];

% Fit linear model (requires Statistics Toolbox):
% mdl = fitlm(design.runSheet, 'Response1 ~ Temp + Time + Temp:Time');
% disp(mdl);

% Plot residuals:
% plotResiduals(mdl);

%% ========================================================================
%  END OF QUICK REFERENCE
%  ========================================================================
