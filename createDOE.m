function [design, report, cfg] = createDOE(factors, varargin)
%CREATEDOE Create Design of Experiments with various design types
%
%   [design, report, cfg] = createDOE(factors, Name, Value, ...)
%
%   INPUTS:
%       factors - Struct array with fields:
%                 .name, .type ("continuous" or "categorical")
%                 For continuous: .low, .high, .units
%                 For categorical: .levels
%
%   NAME-VALUE PAIRS:
%       Design parameters:
%           type           - Design type: 'fullfact', 'fracfact', 'ccd', 
%                           'boxbehnken', 'lhs' (default: 'ccd')
%           replicates     - Number of replicates (default: 1)
%           centerPoints   - Number of center points (default: 0)
%           randomize      - Randomize run order (default: true)
%           randomSeed     - Random seed for reproducibility (default: 42)
%           numBlocks      - Number of blocks (default: 1)
%           blockingMethod - 'sequential' or 'random' (default: 'sequential')
%
%       Design-specific:
%           fracGenerators - Fractional factorial generators (default: 'a b c d ab ac')
%           ccdAlpha       - CCD alpha: 'rotatable', 'orthogonal', 'face', or numeric
%           ccdFaceType    - CCD type: 'circumscribed', 'inscribed', 'faced'
%           lhsSamples     - Number of LHS samples (default: 30)
%           lhsCriterion   - LHS criterion: 'maximin', 'correlation', etc.
%
%       Display:
%           title          - Display title for headers (default: project name or 'DOE')
%
%       Export:
%           projectName    - Project name (default: 'DOE_Project')
%           owner          - Owner name
%           notes          - Project notes
%           exportFolder   - Export folder path
%           exportBaseFilename - Base filename for exports
%           exportCSV      - Export to CSV (default: false)
%           exportMAT      - Export to MAT (default: true)
%
%   OUTPUTS:
%       design - Struct with .runSheet (table) and .meta
%       report - Struct with design summary
%       cfg    - Configuration struct (for reproducibility)

%% Parse inputs
p = inputParser;
p.FunctionName = "createDOE";

% Required
addRequired(p, "factors", @(x) isstruct(x) && ~isempty(x));

% Design
addParameter(p, "type", "ccd");
addParameter(p, "replicates", 1, @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(p, "centerPoints", 0, @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(p, "randomize", true, @islogical);
addParameter(p, "randomSeed", 42, @isnumeric);
addParameter(p, "numBlocks", 1, @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(p, "blockingMethod", "sequential");

% Fractional factorial
addParameter(p, "fracGenerators", "a b c d ab ac");

% CCD
addParameter(p, "ccdAlpha", "rotatable");
addParameter(p, "ccdFaceType", "circumscribed");

% LHS
addParameter(p, "lhsSamples", 30, @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(p, "lhsCriterion", "maximin");

% Constraints
addParameter(p, "constraints", {}, @iscell);

% Project/export
addParameter(p, "projectName", "DOE_Project");
addParameter(p, "title", "");  % Display title
addParameter(p, "owner", "");
addParameter(p, "notes", "");
addParameter(p, "exportFolder", "");
addParameter(p, "exportBaseFilename", "");
addParameter(p, "exportCSV", false, @islogical);
addParameter(p, "exportMAT", true, @islogical);

parse(p, factors, varargin{:});
S = p.Results;

%% Build configuration struct
cfg = struct();
cfg.project.name      = string(S.projectName);
cfg.project.owner     = string(S.owner);
cfg.project.createdOn = datetime("now");
cfg.project.notes     = string(S.notes);

% Set display title - use provided title, else project name, else 'DOE'
if strlength(string(S.title)) > 0
    cfg.project.displayTitle = string(S.title);
elseif strlength(string(S.projectName)) > 0 && ~strcmp(S.projectName, "DOE_Project")
    cfg.project.displayTitle = string(S.projectName);
else
    cfg.project.displayTitle = "DOE";
end

cfg.factors = factors;

cfg.design.type            = string(S.type);
cfg.design.replicates      = round(S.replicates);
cfg.design.centerPoints    = round(S.centerPoints);
cfg.design.randomize       = S.randomize;
cfg.design.randomSeed      = S.randomSeed;
cfg.design.numBlocks       = round(S.numBlocks);
cfg.design.blockingMethod  = string(S.blockingMethod);

cfg.design.frac.generators = string(S.fracGenerators);
cfg.design.ccd.alpha       = S.ccdAlpha;
cfg.design.ccd.facetype    = string(S.ccdFaceType);
cfg.design.lhs.samples     = round(S.lhsSamples);
cfg.design.lhs.criterion   = string(S.lhsCriterion);

cfg.constraints.enabled = ~isempty(S.constraints);
cfg.constraints.funs     = S.constraints;

% Export defaults
if strlength(string(S.exportFolder)) == 0
    cfg.export.folder = fullfile(pwd, "DOE_Output");
else
    cfg.export.folder = string(S.exportFolder);
end

if strlength(string(S.exportBaseFilename)) == 0
    cfg.export.baseFilename = cfg.project.name;
else
    cfg.export.baseFilename = string(S.exportBaseFilename);
end

cfg.export.writeCSV = S.exportCSV;
cfg.export.writeMAT = S.exportMAT;

%% Create DOE
fprintf('\n');
fprintf('========================================\n');
fprintf('%s\n', upper(char(cfg.project.displayTitle)));
fprintf('========================================\n');
fprintf('Creating %s design...\n', upper(char(cfg.design.type)));
[design, report] = buildDOE(cfg);

%% Optional export
if cfg.export.writeCSV || cfg.export.writeMAT
    doe_export(design, report, cfg);
end

fprintf('\nDOE creation complete!\n');
fprintf('Total runs: %d\n', report.nRuns);

end

%% =====================================================================
%  CORE DOE BUILDER
%% =====================================================================
function [design, report] = buildDOE(cfg)

    tStart = tic;
    
    % Validate and categorize factors
    [factors, idxCont, idxCat] = validateFactors(cfg.factors);
    kCont = numel(idxCont);
    kCat  = numel(idxCat);
    
    fprintf('  Continuous factors: %d\n', kCont);
    fprintf('  Categorical factors: %d\n', kCat);
    
    % Generate coded design matrix for continuous factors
    designType = lower(string(cfg.design.type));
    Xc_coded = generateDesignMatrix(designType, kCont, cfg);
    
    fprintf('  Base design runs: %d\n', size(Xc_coded, 1));
    
    % Generate categorical combinations
    C = generateCategoricalCombinations(factors, idxCat);
    
    % Cross continuous and categorical
    [Xc_coded, C] = crossDesigns(Xc_coded, C, kCont, kCat);
    
    fprintf('  After crossing: %d runs\n', size(Xc_coded, 1));
    
    % Convert coded to actual values
    Xc_actual = codedToActual(Xc_coded, factors, idxCont);
    
    % Build run sheet
    runSheet = buildRunSheet(Xc_coded, Xc_actual, C, factors, idxCont, idxCat);
    
    % Apply constraints
    runSheet = applyConstraints(runSheet, cfg);
    
    % Add replicates
    runSheet = addReplicates(runSheet, cfg.design.replicates);
    
    % Add blocking
    runSheet = addBlocking(runSheet, cfg);
    
    % Randomize
    if cfg.design.randomize
        runSheet = randomizeRuns(runSheet, cfg.design.randomSeed);
    end
    
    % Re-number runs
    runSheet.Run = (1:height(runSheet))';
    
    % Add response placeholders
    runSheet.Response1 = nan(height(runSheet), 1);
    runSheet.Response2 = nan(height(runSheet), 1);
    
    % Assemble outputs
    design = struct();
    design.runSheet = runSheet;
    design.meta     = cfg.project;
    design.cfg      = cfg;
    
    report = struct();
    report.title       = cfg.project.displayTitle;
    report.project     = cfg.project.name;
    report.designType  = upper(string(cfg.design.type));
    report.nRuns       = height(runSheet);
    report.kCont       = kCont;
    report.kCat        = kCat;
    report.blocks      = cfg.design.numBlocks;
    report.replicates  = cfg.design.replicates;
    report.randomized  = cfg.design.randomize;
    report.seconds     = toc(tStart);
    
    report.summaryText = sprintf([ ...
        '\n========================================\n' ...
        '%s - SUMMARY\n' ...
        '========================================\n' ...
        'Project: %s\n' ...
        'Design type: %s\n' ...
        'Continuous factors: %d\n' ...
        'Categorical factors: %d\n' ...
        'Total runs: %d\n' ...
        'Blocks: %d\n' ...
        'Replicates: %d\n' ...
        'Randomized: %s\n' ...
        'Build time: %.3f s\n' ...
        '========================================\n'], ...
        upper(char(cfg.project.displayTitle)), ...
        report.project, report.designType, ...
        report.kCont, report.kCat, report.nRuns, ...
        report.blocks, report.replicates, ...
        string(report.randomized), report.seconds);

end

%% =====================================================================
%  HELPER FUNCTIONS
%% =====================================================================

function [factors, idxCont, idxCat] = validateFactors(factors)
    % Validate factor definitions
    names = string({factors.name});
    if numel(unique(names)) ~= numel(names)
        error('DOE:DuplicateNames', 'Factor names must be unique.');
    end
    
    idxCont = find(arrayfun(@(f) strcmpi(f.type, "continuous"), factors));
    idxCat  = find(arrayfun(@(f) strcmpi(f.type, "categorical"), factors));
    
    for i = 1:numel(factors)
        f = factors(i);
        if strcmpi(f.type, "continuous")
            if isempty(f.low) || isempty(f.high) || ~isnumeric(f.low) || ~isnumeric(f.high)
                error('DOE:InvalidContinuous', ...
                    'Continuous factor "%s" must have numeric low/high.', f.name);
            end
            if f.high <= f.low
                error('DOE:InvalidRange', ...
                    'Continuous factor "%s" must satisfy high > low.', f.name);
            end
        elseif strcmpi(f.type, "categorical")
            if isempty(f.levels)
                error('DOE:InvalidCategorical', ...
                    'Categorical factor "%s" must define levels.', f.name);
            end
        else
            error('DOE:UnknownType', ...
                'Factor "%s" has unknown type "%s".', f.name, f.type);
        end
    end
end

function Xc = generateDesignMatrix(designType, kCont, cfg)
    % Generate coded design matrix based on design type
    
    switch designType
        case "fullfact"
            if kCont == 0
                Xc = zeros(1, 0);
            else
                Xc = 2 * fullfact(2 * ones(1, kCont)) - 3; % {-1, +1}
            end
            if cfg.design.centerPoints > 0 && kCont > 0
                Xc = [Xc; zeros(cfg.design.centerPoints, kCont)];
            end
            
        case "fracfact"
            if kCont == 0
                error('DOE:NoFactors', 'Fractional factorial requires at least 1 continuous factor.');
            end
            Xc = fracfact(char(cfg.design.frac.generators));
            if cfg.design.centerPoints > 0
                Xc = [Xc; zeros(cfg.design.centerPoints, kCont)];
            end
            
        case "ccd"
            if kCont == 0
                error('DOE:NoFactors', 'CCD requires at least 1 continuous factor.');
            end
            
            % Handle alpha parameter
            if isstring(cfg.design.ccd.alpha) || ischar(cfg.design.ccd.alpha)
                alphaStr = lower(char(string(cfg.design.ccd.alpha)));
                switch alphaStr
                    case 'rotatable'
                        alphaVal = sqrt(kCont);
                    case 'orthogonal'
                        alphaVal = (kCont * (1 + sqrt(1 + 8*kCont)) / 4)^0.25;
                    case 'face'
                        alphaVal = 1;
                    otherwise
                        alphaVal = sqrt(kCont);  % default rotatable
                end
            else
                alphaVal = cfg.design.ccd.alpha;
            end
            
            % Center points
            cp = max(0, round(cfg.design.centerPoints));
            
            % Try to use ccdesign if available, otherwise build manually
            useManualCCD = false;
            
            if exist('ccdesign', 'file') == 2
                % Try multiple ccdesign syntaxes for compatibility
                ccdSuccess = false;
                
                % Method 1: Name-value pairs (Capital - newer MATLAB)
                try
                    Xc = ccdesign(kCont, 'Alpha', alphaVal, 'Center', [cp cp], ...
                        'Type', char(cfg.design.ccd.facetype));
                    ccdSuccess = true;
                catch
                    % Method 2: Name-value pairs (lowercase - older MATLAB)
                    try
                        Xc = ccdesign(kCont, 'alpha', alphaVal, 'center', [cp cp], ...
                            'type', char(cfg.design.ccd.facetype));
                        ccdSuccess = true;
                    catch
                        % Method 3: Mixed positional and name-value
                        try
                            Xc = ccdesign(kCont, char(cfg.design.ccd.facetype), ...
                                'alpha', alphaVal, 'center', [cp cp]);
                            ccdSuccess = true;
                        catch
                            % ccdesign exists but syntax failed, use manual
                            useManualCCD = true;
                        end
                    end
                end
            else
                % ccdesign not available, use manual
                useManualCCD = true;
            end
            
            % Manual CCD builder (fallback)
            if useManualCCD
                fprintf('  Using manual CCD builder (ccdesign unavailable or incompatible)\n');
                
                % Factorial points: 2^k corners
                Xf = 2 * fullfact(2 * ones(1, kCont)) - 3;  % {-1, +1}
                
                % Axial points: ±alpha along each axis
                Xa = zeros(2 * kCont, kCont);
                for i = 1:kCont
                    Xa(2*i-1, i) = -alphaVal;
                    Xa(2*i, i) = alphaVal;
                end
                
                % Center points
                Xcp = zeros(cp, kCont);
                
                % Combine all points
                Xc = [Xf; Xa; Xcp];
                
                fprintf('  Manual CCD: %d factorial + %d axial + %d center = %d total\n', ...
                    size(Xf,1), size(Xa,1), cp, size(Xc,1));
            end
            
        case "boxbehnken"
            if kCont < 3
                error('DOE:TooFewFactors', 'Box-Behnken requires >= 3 continuous factors.');
            end
            
            if exist('bbdesign', 'file') == 2
                Xc = bbdesign(kCont);
                if cfg.design.centerPoints > 0
                    Xc = [Xc; zeros(cfg.design.centerPoints, kCont)];
                end
            else
                error('DOE:MissingToolbox', ...
                    ['Box-Behnken design requires Statistics and Machine Learning Toolbox.\n' ...
                     'The bbdesign function was not found.\n' ...
                     'Try using type=''ccd'' or type=''fullfact'' instead.']);
            end
            
        case "lhs"
            if kCont == 0
                error('DOE:NoFactors', 'LHS requires at least 1 continuous factor.');
            end
            n = cfg.design.lhs.samples;
            crit = char(cfg.design.lhs.criterion);
            U = lhsdesign(n, kCont, 'criterion', crit);
            Xc = 2 * U - 1; % Map [0,1] to [-1,1]
            
        otherwise
            error('DOE:UnknownDesign', 'Unknown design type: %s', designType);
    end
end

function C = generateCategoricalCombinations(factors, idxCat)
    % Generate all combinations of categorical levels
    kCat = numel(idxCat);
    
    if kCat == 0
        C = strings(1, 0);
        return;
    end
    
    catLevels = cell(1, kCat);
    for j = 1:kCat
        lv = factors(idxCat(j)).levels;
        % Ensure we have a string vector
        if isstring(lv) || ischar(lv)
            catLevels{j} = string(lv(:));
        elseif iscell(lv)
            catLevels{j} = string(lv(:));
        else
            catLevels{j} = string(lv(:));
        end
    end
    
    % Use allcomb to generate combinations
    C = allcomb(catLevels{:});
end

function [Xnew, Cnew] = crossDesigns(X, C, kCont, kCat)
    % Cross product of continuous and categorical designs
    
    if kCont > 0 && kCat > 0
        n = size(X, 1);
        m = size(C, 1);
        Xnew = repmat(X, m, 1);
        Cnew = repelem(C, n, 1);
    elseif kCont == 0 && kCat > 0
        Xnew = zeros(size(C, 1), 0);
        Cnew = C;
    else
        Xnew = X;
        Cnew = strings(size(X, 1), 0);
    end
end

function Xa = codedToActual(Xc, factors, idxCont)
    % Convert coded values [-1, 1] to actual values
    kCont = numel(idxCont);
    n = size(Xc, 1);
    Xa = zeros(n, kCont);
    
    for j = 1:kCont
        f = factors(idxCont(j));
        low = f.low;
        high = f.high;
        Xa(:, j) = (Xc(:, j) + 1) / 2 * (high - low) + low;
    end
end

function runSheet = buildRunSheet(Xc_coded, Xc_actual, C, factors, idxCont, idxCat)
    % Build the initial run sheet table
    
    kCont = numel(idxCont);
    kCat = numel(idxCat);
    
    runSheet = table();
    runSheet.Run = (1:size(Xc_coded, 1))';
    
    % Add continuous factor columns
    for j = 1:kCont
        f = factors(idxCont(j));
        runSheet.(f.name + "_coded") = Xc_coded(:, j);
        runSheet.(f.name) = Xc_actual(:, j);
    end
    
    % Add categorical factor columns
    for j = 1:kCat
        f = factors(idxCat(j));
        runSheet.(f.name) = string(C(:, j));
    end
end

function runSheet = applyConstraints(runSheet, cfg)
    % Apply constraint functions to filter runs
    
    if ~cfg.constraints.enabled || isempty(cfg.constraints.funs)
        return;
    end
    
    fprintf('  Applying %d constraint(s)...\n', numel(cfg.constraints.funs));
    
    maskKeep = true(height(runSheet), 1);
    for i = 1:numel(cfg.constraints.funs)
        fun = cfg.constraints.funs{i};
        ok = fun(runSheet);
        ok = logical(ok);
        if isscalar(ok)
            ok = repmat(ok, height(runSheet), 1);
        end
        maskKeep = maskKeep & ok;
    end
    
    nRemoved = sum(~maskKeep);
    if nRemoved > 0
        fprintf('  Removed %d runs due to constraints\n', nRemoved);
    end
    
    runSheet = runSheet(maskKeep, :);
end

function runSheet = addReplicates(runSheet, nRep)
    % Add replicates to the design
    
    nRep = max(1, round(nRep));
    if nRep > 1
        runSheet = repmat(runSheet, nRep, 1);
        runSheet.Replicate = repelem((1:nRep)', height(runSheet) / nRep);
    else
        runSheet.Replicate = ones(height(runSheet), 1);
    end
end

function runSheet = addBlocking(runSheet, cfg)
    % Add blocking structure to the design
    
    nb = max(1, round(cfg.design.numBlocks));
    nRuns = height(runSheet);
    
    if nb > 1
        switch lower(string(cfg.design.blockingMethod))
            case "sequential"
                edges = round(linspace(0, nRuns, nb + 1));
                B = zeros(nRuns, 1);
                for b = 1:nb
                    B((edges(b) + 1):edges(b + 1)) = b;
                end
                runSheet.Block = B;
                
            case "random"
                rng(cfg.design.randomSeed);
                runSheet.Block = randi(nb, nRuns, 1);
                
            otherwise
                error('DOE:UnknownBlocking', 'Unknown blockingMethod: %s', cfg.design.blockingMethod);
        end
    else
        runSheet.Block = ones(nRuns, 1);
    end
end

function runSheet = randomizeRuns(runSheet, seed)
    % Randomize run order within blocks
    
    rng(seed);
    blocks = unique(runSheet.Block);
    idxAll = [];
    
    for b = blocks'
        idx = find(runSheet.Block == b);
        idx = idx(randperm(numel(idx)));
        idxAll = [idxAll; idx]; %#ok<AGROW>
    end
    
    runSheet = runSheet(idxAll, :);
end
