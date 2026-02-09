%% ========================================================================
%  DOE TOOLKIT — QUICK REFERENCE
%  ========================================================================
%  Patterns and recipes for createDOE.
%
%  KEY CHANGE FROM THE OLD TOOLKIT:
%    When mixing categorical and continuous factors, use a CELL ARRAY
%    of structs (because MATLAB struct arrays need identical fields).
%
%    All-continuous:   f = [ struct(...); struct(...) ];   % struct array OK
%    Mixed types:      f = { struct(...); struct(...) };   % cell array
%  ========================================================================

%% PATTERN 1 — Simple Full Factorial (2 continuous factors)
% -------------------------------------------------------------------------
factors = [ struct('name','Temp','type','continuous','low',20,'high',80)
            struct('name','Time','type','continuous','low',10,'high',60) ];

[design, report] = createDOE(factors, 'type','fullfact');
%  → 4 runs (2^2)

%% PATTERN 2 — CCD for Response-Surface Modelling (3 factors)
% -------------------------------------------------------------------------
factors = [ struct('name','Temp','type','continuous','low',60,'high',90)
            struct('name','pH',  'type','continuous','low',6, 'high',9)
            struct('name','Conc','type','continuous','low',0.1,'high',1.0) ];

[design, report] = createDOE(factors, ...
    'type','ccd', ...
    'ccdAlpha','rotatable', ...
    'centerPoints',6, ...
    'exportCSV',true);
%  → 20 runs (8 factorial + 6 axial + 6 center)

%% PATTERN 3 — Mixed Categorical + Continuous
% -------------------------------------------------------------------------
factors = {
    struct('name','Catalyst','type','categorical','levels',["A","B","C"])
    struct('name','Temp','type','continuous','low',25,'high',75)
    struct('name','Time','type','continuous','low',30,'high',120)
};

[design, report] = createDOE(factors, ...
    'type','fullfact', ...
    'centerPoints',3);
%  → (4 + 3) × 3 = 21 runs

%% PATTERN 4 — Screening (Fractional Factorial)
% -------------------------------------------------------------------------
factors = [ struct('name','A','type','continuous','low',10,'high',50)
            struct('name','B','type','continuous','low',1, 'high',5)
            struct('name','C','type','continuous','low',100,'high',200)
            struct('name','D','type','continuous','low',0.1,'high',1.0) ];

[design, report] = createDOE(factors, ...
    'type','fracfact', ...
    'fracGenerators','a b c abc');
%  → 8 runs (2^3 base, D = A×B×C)

%% PATTERN 5 — Box–Behnken (efficient RSM, 3+ factors)
% -------------------------------------------------------------------------
factors = [ struct('name','A','type','continuous','low',10,'high',30)
            struct('name','B','type','continuous','low',5, 'high',15)
            struct('name','C','type','continuous','low',100,'high',200) ];

[design, report] = createDOE(factors, ...
    'type','bbd', ...
    'centerPoints',5);
%  → 18 runs  (vs 20 for CCD)

%% PATTERN 6 — Latin Hypercube for Computer Experiments
% -------------------------------------------------------------------------
factors = [ struct('name','X1','type','continuous','low',0,'high',1)
            struct('name','X2','type','continuous','low',0,'high',1)
            struct('name','X3','type','continuous','low',0,'high',1) ];

[design, report] = createDOE(factors, ...
    'type','lhs', ...
    'lhsSamples',100);
%  → 100 space-filling points

%% PATTERN 7 — Plackett–Burman (screening many factors)
% -------------------------------------------------------------------------
% Generates a very small run count (next multiple of 4 above k+1)
nFactors = 7;
factors = [];
for i = 1:nFactors
    factors = [factors; struct('name',sprintf('F%d',i), ...
        'type','continuous','low',0,'high',10)]; %#ok<AGROW>
end

[design, report] = createDOE(factors, 'type','pbdesign');
%  → 8 runs for 7 factors

%% PATTERN 8 — Constraints
% -------------------------------------------------------------------------
factors = [ struct('name','Temp','type','continuous','low',60,'high',100)
            struct('name','Press','type','continuous','low',1,'high',10) ];

constraints = {};
constraints{1} = @(T) ~(T.Temp > 90 & T.Press > 8);  % safety limit
constraints{2} = @(T) T.Temp >= 65 | T.Press >= 3;    % minimum activation

[design, report] = createDOE(factors, ...
    'type','fullfact', ...
    'constraints',constraints);

%% PATTERN 9 — Blocking + Replicates
% -------------------------------------------------------------------------
factors = [ struct('name','A','type','continuous','low',0,'high',100)
            struct('name','B','type','continuous','low',0,'high',50) ];

[design, report] = createDOE(factors, ...
    'type','ccd', ...
    'centerPoints',4, ...
    'replicates',2, ...
    'numBlocks',3, ...
    'randomize',true);

%% CHOOSING A DESIGN TYPE
% -------------------------------------------------------------------------
%   Full Factorial   2–4 factors, need all interactions
%   Frac. Factorial  5+ factors, screening — picks out main effects
%   CCD              Response surface, 2–5 factors, quadratic model
%   Box–Behnken      Response surface, 3+ factors, fewer runs than CCD
%   Plackett–Burman  Main-effects-only screening, 5+ factors, minimal runs
%   LHS              Computer experiments, many factors, space-filling

%% TYPE ALIASES (all resolve to the canonical name)
% -------------------------------------------------------------------------
%   'fullfact', 'full', 'ff', 'factorial'
%   'fracfact', 'frac', 'fractional'
%   'ccd', 'ccdesign', 'centralcomposite'
%   'bbd', 'boxbehnken', 'bb', 'box-behnken'
%   'lhs', 'latinhypercube'
%   'pb', 'pbdesign', 'plackett', 'plackettburman'

%% LOADING SAVED DESIGNS
% -------------------------------------------------------------------------
% From MAT:
%   load('./DOE_Output/MyProject_bundle.mat');
%   disp(design.runSheet);
%
% From CSV:
%   T = readtable('./DOE_Output/MyProject_runsheet.csv');

%% ANALYSING RESULTS
% -------------------------------------------------------------------------
% Fill in measured values:
%   design.runSheet.Response1 = measuredValues;
%
% Fit a model (requires Statistics Toolbox):
%   mdl = fitlm(design.runSheet, 'Response1 ~ Temp*Time + Temp^2 + Time^2');
%   disp(mdl);
%   plotResiduals(mdl);
