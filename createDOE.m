function [design, report, cfg] = createDOE(factors, varargin)

%% ---------- Parse inputs ----------
p = inputParser;
p.FunctionName = "createDOE";

% Required
addRequired(p, "factors", @(x) isstruct(x) && ~isempty(x));

% Design
addParameter(p, "type", "ccd", @(s) isstring(s) || ischar(s));
addParameter(p, "replicates", 1, @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(p, "centerPoints", 0, @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(p, "randomize", true, @(x) islogical(x) && isscalar(x));
addParameter(p, "randomSeed", 42, @(x) isnumeric(x) && isscalar(x));
addParameter(p, "numBlocks", 1, @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(p, "blockingMethod", "sequential", @(s) isstring(s) || ischar(s));

% Fractional factorial
addParameter(p, "fracGenerators", "a b c d ab ac", @(s) isstring(s) || ischar(s));

% CCD
addParameter(p, "ccdAlpha", "rotatable", @(x) (isstring(x)||ischar(x)||isnumeric(x)));
addParameter(p, "ccdFaceType", "circumscribed", @(s) isstring(s) || ischar(s));

% LHS
addParameter(p, "lhsSamples", 30, @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(p, "lhsCriterion", "maximin", @(s) isstring(s) || ischar(s));

% Constraints
addParameter(p, "constraints", {}, @(c) iscell(c));

% Project/export
addParameter(p, "projectName", "DOE_Project", @(s) isstring(s) || ischar(s));
addParameter(p, "owner", "", @(s) isstring(s) || ischar(s));
addParameter(p, "notes", "", @(s) isstring(s) || ischar(s));
addParameter(p, "exportFolder", "", @(s) isstring(s) || ischar(s));
addParameter(p, "exportBaseFilename", "", @(s) isstring(s) || ischar(s));
addParameter(p, "exportCSV", false, @(x) islogical(x) && isscalar(x));
addParameter(p, "exportMAT", true, @(x) islogical(x) && isscalar(x));

parse(p, factors, varargin{:});
S = p.Results;

%% ---------- Build cfg struct (single source of truth) ----------
cfg = struct();

cfg.project.name      = string(S.projectName);
cfg.project.owner     = string(S.owner);
cfg.project.createdOn = datetime("now");
cfg.project.notes     = string(S.notes);

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

% export defaults
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

%% ---------- Create DOE ----------
[design, report] = doe_build_core(cfg);

%% ---------- Optional export ----------
if cfg.export.writeCSV || cfg.export.writeMAT
    doe_export(design, report, cfg);
end

end

%% =====================================================================
% CORE BUILDER (kept inside this file so createDOE is one-stop)
%% =====================================================================
function [design, report] = doe_build_core(cfg)

tStart = tic;

% Validate
[factors, idxCont, idxCat] = local_validateFactors(cfg.factors);
kCont = numel(idxCont);
kCat  = numel(idxCat);

designType = lower(string(cfg.design.type));
Xc_coded = [];

switch designType
    case "fullfact"
        if kCont == 0
            Xc_coded = zeros(1,0);
        else
            Xc_coded = 2*fullfact(2*ones(1,kCont)) - 3; % {-1,+1}
        end
        if cfg.design.centerPoints > 0 && kCont > 0
            Xc_coded = [Xc_coded; zeros(cfg.design.centerPoints, kCont)];
        end

    case "fracfact"
        if kCont == 0
            error("Fractional factorial requires at least 1 continuous factor.");
        end
        Xc_coded = fracfact(char(cfg.design.frac.generators));
        if cfg.design.centerPoints > 0
            Xc_coded = [Xc_coded; zeros(cfg.design.centerPoints, kCont)];
        end

    case "ccd"
        if kCont == 0
            error("CCD requires at least 1 continuous factor.");
        end
        local_requireStatsToolbox("CCD (ccdesign)");

        alpha = cfg.design.ccd.alpha;
        if isstring(alpha) || ischar(alpha)
            alphaArg = char(string(alpha));
        else
            alphaArg = alpha;
        end

        cp = max(0, round(cfg.design.centerPoints));
        centerVec = [cp cp];
        Xc_coded = ccdesign(kCont, "alpha", alphaArg, "center", centerVec, "type", char(cfg.design.ccd.facetype));

    case "boxbehnken"
        if kCont < 3
            error("Box-Behnken requires >= 3 continuous factors.");
        end
        local_requireStatsToolbox("Box-Behnken (bbdesign)");

        Xc_coded = bbdesign(kCont);
        if cfg.design.centerPoints > 0
            Xc_coded = [Xc_coded; zeros(cfg.design.centerPoints, kCont)];
        end

    case "lhs"
        if kCont == 0
            error("LHS requires at least 1 continuous factor.");
        end
        n = cfg.design.lhs.samples;
        crit = char(cfg.design.lhs.criterion);
        U = lhsdesign(n, kCont, "criterion", crit);
        Xc_coded = 2*U - 1; % map to [-1,1]

    otherwise
        error("Unknown design type: %s", cfg.design.type);
end

% Categorical combinations
if kCat > 0
    catLevels = cell(1, kCat);
    for j = 1:kCat
        lv = factors(idxCat(j)).levels;
        catLevels{j} = string(lv(:))';
    end
    C = allcomb(catLevels{:});  % cell/strings OK
else
    C = strings(size(Xc_coded,1), 0);
end

% Cross product
if kCont > 0 && kCat > 0
    [Xc_coded, C] = local_crossProduct(Xc_coded, C);
elseif kCont == 0 && kCat > 0
    Xc_coded = zeros(size(C,1), 0);
end

% Coded -> actual
Xc_actual = local_codedToActual(Xc_coded, factors, idxCont);

% Build run sheet
runSheet = table();
runSheet.Run = (1:size(Xc_coded,1))';

% Continuous columns
for j = 1:kCont
    f = factors(idxCont(j));
    runSheet.(f.name + "_coded") = Xc_coded(:,j);
    runSheet.(f.name)           = Xc_actual(:,j);
end

% Categorical columns
for j = 1:kCat
    f = factors(idxCat(j));
    runSheet.(f.name) = string(C(:,j));
end

% Constraints
maskKeep = true(height(runSheet), 1);
if cfg.constraints.enabled && ~isempty(cfg.constraints.funs)
    for i = 1:numel(cfg.constraints.funs)
        fun = cfg.constraints.funs{i};
        ok = fun(runSheet);
        ok = logical(ok);
        if isscalar(ok)
            ok = repmat(ok, height(runSheet), 1);
        end
        maskKeep = maskKeep & ok;
    end
end
runSheet = runSheet(maskKeep, :);

% Replicates
rep = max(1, round(cfg.design.replicates));
if rep > 1
    runSheet = repmat(runSheet, rep, 1);
    runSheet.Replicate = repelem((1:rep)', height(runSheet)/rep);
else
    runSheet.Replicate = ones(height(runSheet),1);
end

% Blocking
nb = max(1, round(cfg.design.numBlocks));
if nb > 1
    runSheet.Block = local_assignBlocks(height(runSheet), nb, cfg.design.blockingMethod, cfg.design.randomSeed);
else
    runSheet.Block = ones(height(runSheet),1);
end

% Randomize within blocks
if cfg.design.randomize
    rng(cfg.design.randomSeed);
    runSheet = local_randomize(runSheet);
end

% Re-number
runSheet.Run = (1:height(runSheet))';

% Response placeholders
runSheet.Response1 = nan(height(runSheet),1);
runSheet.Response2 = nan(height(runSheet),1);

% Outputs
design = struct();
design.runSheet = runSheet;
design.meta     = cfg.project;
design.cfg      = cfg; % embed cfg for reproducibility

report = struct();
report.project     = cfg.project.name;
report.designType  = upper(string(cfg.design.type));
report.nRuns       = height(runSheet);
report.kCont       = kCont;
report.kCat        = kCat;
report.blocks      = nb;
report.replicates  = rep;
report.randomized  = cfg.design.randomize;
report.seconds     = toc(tStart);

report.summaryText = sprintf([ ...
    "DOE created: %s\n" ...
    "Design type: %s\n" ...
    "Continuous: %d | Categorical: %d\n" ...
    "Runs: %d | Blocks: %d | Replicates: %d | Randomized: %d\n" ...
    "Build time: %.3f s\n"], ...
    report.project, report.designType, report.kCont, report.kCat, ...
    report.nRuns, report.blocks, report.replicates, report.randomized, report.seconds);

end

%% ===== Helpers =====
function [factors, idxCont, idxCat] = local_validateFactors(factors)
names = string({factors.name});
if numel(unique(names)) ~= numel(names)
    error("Factor names must be unique.");
end

idxCont = find(arrayfun(@(f) strcmpi(f.type,"continuous"), factors));
idxCat  = find(arrayfun(@(f) strcmpi(f.type,"categorical"), factors));

for i = 1:numel(factors)
    f = factors(i);
    if strcmpi(f.type, "continuous")
        if isempty(f.low) || isempty(f.high) || ~isnumeric(f.low) || ~isnumeric(f.high)
            error("Continuous factor '%s' must have numeric low/high.", f.name);
        end
        if ~(f.high > f.low)
            error("Continuous factor '%s' must satisfy high > low.", f.name);
        end
    elseif strcmpi(f.type, "categorical")
        if isempty(f.levels)
            error("Categorical factor '%s' must define levels.", f.name);
        end
    else
        error("Factor '%s' has unknown type '%s'.", f.name, f.type);
    end
end
end

function Xa = local_codedToActual(Xc, factors, idxCont)
kCont = numel(idxCont);
n = size(Xc,1);
Xa = zeros(n, kCont);
for j = 1:kCont
    f = factors(idxCont(j));
    low = f.low; high = f.high;
    Xa(:,j) = (Xc(:,j) + 1)/2 * (high - low) + low;
end
end

function [Xnew, Cnew] = local_crossProduct(X, C)
n = size(X,1);
m = size(C,1);
Xnew = kron(ones(m,1), X);
Cnew = repmat(C, n, 1);
end

function B = local_assignBlocks(nRuns, nBlocks, method, seed)
B = zeros(nRuns,1);
switch lower(string(method))
    case "sequential"
        edges = round(linspace(0, nRuns, nBlocks+1));
        for b = 1:nBlocks
            B((edges(b)+1):edges(b+1)) = b;
        end
    case "random"
        rng(seed);
        B = randi(nBlocks, nRuns, 1);
    otherwise
        error("Unknown blockingMethod: %s", method);
end
end

function T = local_randomize(T)
blocks = unique(T.Block);
idxAll = [];
for b = blocks'
    idx = find(T.Block == b);
    idx = idx(randperm(numel(idx)));
    idxAll = [idxAll; idx]; %#ok<AGROW>
end
T = T(idxAll,:);
end

function local_requireStatsToolbox(featureName)
if ~(exist("ccdesign","file")==2 || exist("bbdesign","file")==2)
    error("This design requires Statistics and Machine Learning Toolbox (%s).", featureName);
end
end
