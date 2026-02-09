function [design, report, cfg] = createDOE(factors, varargin)
%CREATEDOE  Create a Design of Experiments run-sheet.
%
%   [design, report, cfg] = createDOE(factors, Name, Value, ...)
%
%   Self-contained — no Statistics Toolbox required.
%
%   FACTORS (struct array or cell array of structs)
%   ───────────────────────────────────────────────
%     .name   - unique factor name (valid MATLAB variable name)
%     .type   - 'continuous' | 'categorical'
%
%     Continuous: .low, .high          (required)
%                 .units               (optional, default '-')
%     Categorical: .levels             (required, e.g. ["A","B","C"])
%
%   NAME-VALUE OPTIONS
%   ──────────────────
%   Design type
%     'type'         - 'fullfact' | 'fracfact' | 'ccd' | 'bbd' |
%                      'bbd_ccd' | 'lhs' | 'pbdesign'
%                      (default 'fullfact')
%
%   Categorical handling
%     'catMethod'    - 'cross'  : build continuous design, then cross with
%                                 every categorical combination (default)
%                      'encode' : map categorical levels to coded values
%                                 (-1, 0, +1) and include them inside the
%                                 design matrix.  Requires exactly 3 levels
%                                 per categorical factor for BBD/CCD.
%
%   General
%     'centerPoints' - # center points              (default 0)
%     'replicates'   - # complete replicates         (default 1)
%     'randomize'    - true/false                    (default true)
%     'randomSeed'   - integer seed                  (default 42)
%     'numBlocks'    - # blocks (sequential split)   (default 1)
%
%   CCD / BBD+CCD specific
%     'ccdAlpha'     - 'rotatable' | 'face' | numeric  (default 'rotatable')
%
%   Fractional-factorial
%     'fracGenerators' - generator string            (default '')
%
%   LHS
%     'lhsSamples'   - # samples                     (default 50)
%     'lhsMaxIter'   - maximin iterations             (default 20)
%
%   Constraints
%     'constraints'  - cell array of @(T) ... function handles
%
%   Export / metadata
%     'title', 'projectName', 'owner', 'notes'
%     'exportFolder', 'exportCSV', 'exportMAT'
%     'responses'    - cell array of response column names
%
%   OUTPUTS
%     design.runSheet - table (Run, Block, Replicate, factors, responses)
%     report          - summary struct
%     cfg             - full config for reproducibility
%
%   EXAMPLES
%     % 1) Full factorial
%     f = [ struct('name','Temp','type','continuous','low',60,'high',90)
%           struct('name','Time','type','continuous','low',30,'high',180) ];
%     [d,r] = createDOE(f, 'type','fullfact');
%
%     % 2) CCD with categorical CROSSED (default — full design per level)
%     f = { struct('name','Cat','type','categorical','levels',["A","B","C"])
%           struct('name','Temp','type','continuous','low',60,'high',90)
%           struct('name','Conc','type','continuous','low',0.2,'high',1) };
%     [d,r] = createDOE(f, 'type','ccd', 'centerPoints',6);
%
%     % 3) BBD with categorical ENCODED (fewer runs)
%     [d,r] = createDOE(f, 'type','bbd', 'catMethod','encode', 'centerPoints',3);
%
%     % 4) BBD core + CCD axial augmentation (hybrid)
%     [d,r] = createDOE(f, 'type','bbd_ccd', 'catMethod','encode', ...
%                        'ccdAlpha','rotatable', 'centerPoints',3);
%
%   See also: QuickReference.m

% ======================================================================
%  0. PARSE INPUTS
% ======================================================================
p = inputParser;
p.FunctionName = 'createDOE';

addRequired(p, 'factors',  @(x) ~isempty(x));

addParameter(p, 'type',           'fullfact');
addParameter(p, 'catMethod',      'cross');
addParameter(p, 'centerPoints',   0);
addParameter(p, 'replicates',     1);
addParameter(p, 'randomize',      true);
addParameter(p, 'randomSeed',     42);
addParameter(p, 'numBlocks',      1);

addParameter(p, 'ccdAlpha',       'rotatable');
addParameter(p, 'fracGenerators', '');
addParameter(p, 'lhsSamples',     50);
addParameter(p, 'lhsMaxIter',     20);

addParameter(p, 'constraints',    {});
addParameter(p, 'responses',      {});

addParameter(p, 'title',          '');
addParameter(p, 'projectName',    'DOE_Project');
addParameter(p, 'owner',          '');
addParameter(p, 'notes',          '');
addParameter(p, 'exportFolder',   './DOE_Output');
addParameter(p, 'exportCSV',      false);
addParameter(p, 'exportMAT',      false);

parse(p, factors, varargin{:});
S = p.Results;

% ======================================================================
%  1. VALIDATE & NORMALISE FACTORS
% ======================================================================
if iscell(factors)
    tmp = factors(:);
else
    tmp = num2cell(factors(:));
end
nF = numel(tmp);

allFields = {'name','type','low','high','levels','units'};
defaults  = {'',   '',    [],   [],    [],      '-'};
for i = 1:nF
    s = tmp{i};
    for fi = 1:numel(allFields)
        if ~isfield(s, allFields{fi}) || isempty(s.(allFields{fi}))
            s.(allFields{fi}) = defaults{fi};
        end
    end
    s.name = char(s.name);
    s.type = lower(char(s.type));
    tmp{i} = s;
end
factors = orderfields(tmp{1}, allFields);
for i = 2:nF
    factors(i) = orderfields(tmp{i}, allFields);
end

isCont = arrayfun(@(f) strcmp(f.type,'continuous'), factors);
isCat  = arrayfun(@(f) strcmp(f.type,'categorical'), factors);

if any(~isCont & ~isCat)
    bad = {factors(~isCont & ~isCat).name};
    error('createDOE:badType', ...
        'Unknown factor type for: %s.  Use ''continuous'' or ''categorical''.', ...
        strjoin(bad, ', '));
end

idxCont = find(isCont);
idxCat  = find(isCat);
kCont   = numel(idxCont);
kCat    = numel(idxCat);

names = {factors.name};
if numel(unique(names)) ~= nF
    error('createDOE:dupName', 'Factor names must be unique.');
end

for i = idxCont(:)'
    f = factors(i);
    assert(isnumeric(f.low) && isscalar(f.low), ...
        'createDOE:badLow',  'Factor "%s": .low must be a numeric scalar.', f.name);
    assert(isnumeric(f.high) && isscalar(f.high), ...
        'createDOE:badHigh', 'Factor "%s": .high must be a numeric scalar.', f.name);
    assert(f.high > f.low, ...
        'createDOE:badRange','Factor "%s": high (%.4g) must be > low (%.4g).', ...
        f.name, f.high, f.low);
end
for i = idxCat(:)'
    f = factors(i);
    assert(~isempty(f.levels), ...
        'createDOE:noLevels','Categorical factor "%s" must define .levels.', f.name);
    factors(i).levels = string(f.levels(:));
end

% ======================================================================
%  2. RESOLVE catMethod
% ======================================================================
catMethod = lower(char(S.catMethod));
assert(ismember(catMethod, {'cross','encode'}), ...
    'createDOE:badCatMethod', 'catMethod must be ''cross'' or ''encode''.');

if strcmp(catMethod, 'encode') && kCat > 0
    for i = idxCat(:)'
        nLev = numel(factors(i).levels);
        if nLev ~= 3
            error('createDOE:encodeLevels', ...
                ['catMethod=''encode'' requires exactly 3 levels per categorical ' ...
                 'factor.  Factor "%s" has %d levels.\n' ...
                 'Use catMethod=''cross'' for factors with != 3 levels.'], ...
                factors(i).name, nLev);
        end
    end
end

if strcmp(catMethod, 'encode')
    kDesign = kCont + kCat;
    idxDesign = [idxCont(:); idxCat(:)];
else
    kDesign = kCont;
    idxDesign = idxCont(:);
end

% ======================================================================
%  3. BUILD CONFIGURATION STRUCT
% ======================================================================
cfg = struct();
cfg.project.name        = char(S.projectName);
cfg.project.owner       = char(S.owner);
cfg.project.notes       = char(S.notes);
cfg.project.createdOn   = datestr(now, 'yyyy-mm-dd HH:MM:SS'); %#ok<TNOW1,DATST>
cfg.project.title       = char(S.title);
if isempty(cfg.project.title)
    cfg.project.title = cfg.project.name;
end

cfg.factors      = factors;
cfg.designType   = resolveTypeAlias(S.type);
cfg.catMethod    = catMethod;
cfg.centerPoints = max(0, round(S.centerPoints));
cfg.replicates   = max(1, round(S.replicates));
cfg.randomize    = logical(S.randomize);
cfg.randomSeed   = round(S.randomSeed);
cfg.numBlocks    = max(1, round(S.numBlocks));
cfg.ccdAlpha     = S.ccdAlpha;
cfg.fracGen      = char(S.fracGenerators);
cfg.lhsSamples   = max(1, round(S.lhsSamples));
cfg.lhsMaxIter   = max(1, round(S.lhsMaxIter));
cfg.constraints  = S.constraints;

if isempty(S.responses)
    cfg.responses = {'Response1','Response2'};
else
    cfg.responses = S.responses;
end

cfg.export.folder   = char(S.exportFolder);
cfg.export.csv      = logical(S.exportCSV);
cfg.export.mat      = logical(S.exportMAT);

% ======================================================================
%  4. GENERATE CODED DESIGN MATRIX
% ======================================================================
printHeader(cfg);
fprintf('  Factors: %d continuous, %d categorical  (catMethod: %s)\n', ...
    kCont, kCat, catMethod);

Xc = generateCoded(cfg.designType, kDesign, cfg);

fprintf('  Base design runs: %d  (%d columns)\n', size(Xc,1), size(Xc,2));

% ======================================================================
%  5. HANDLE CATEGORICAL FACTORS
% ======================================================================
if strcmp(catMethod, 'encode')
    Ccat = strings(0,0);
    fprintf('  Categorical factors encoded in design matrix\n');

elseif kCat > 0
    catLevels = cell(1, kCat);
    for j = 1:kCat
        catLevels{j} = factors(idxCat(j)).levels;
    end
    Ccat = localAllComb(catLevels{:});
    fprintf('  Categorical combinations: %d\n', size(Ccat,1));

    if kCont > 0
        nC = size(Xc,1);
        nL = size(Ccat,1);
        Xc   = repmat(Xc, nL, 1);
        Ccat = repelem(Ccat, nC, 1);
        fprintf('  After crossing: %d runs\n', size(Xc,1));
    else
        Xc = zeros(size(Ccat,1), 0);
    end
else
    Ccat = strings(0,0);
    if kCont == 0
        error('createDOE:noFactors', 'Need at least one factor.');
    end
end

nBase = size(Xc, 1);

% ======================================================================
%  6. DECODE TO ACTUAL VALUES & BUILD TABLE
% ======================================================================
T = table();
T.Run       = (1:nBase)';
T.Block     = ones(nBase,1);
T.Replicate = ones(nBase,1);

if strcmp(catMethod, 'encode')
    col = 0;
    for j = 1:numel(idxDesign)
        col = col + 1;
        fi = idxDesign(j);
        f  = factors(fi);
        nm = f.name;

        if strcmp(f.type, 'continuous')
            mid = (f.high + f.low) / 2;
            hwd = (f.high - f.low) / 2;
            T.(nm) = mid + Xc(:,col) * hwd;
            T.([nm '_coded']) = Xc(:,col);
        else
            % coded → label:  -1 → levels(1), 0 → levels(2), +1 → levels(3)
            coded = Xc(:,col);
            T.([nm '_coded']) = coded;
            labels = strings(nBase, 1);
            labels(coded < -0.5) = f.levels(1);
            labels(abs(coded) <= 0.5) = f.levels(2);
            labels(coded > 0.5)  = f.levels(3);
            T.(nm) = labels;
        end
    end
else
    for j = 1:kCont
        f   = factors(idxCont(j));
        nm  = f.name;
        mid = (f.high + f.low) / 2;
        hwd = (f.high - f.low) / 2;
        T.(nm) = mid + Xc(:,j) * hwd;
        T.([nm '_coded']) = Xc(:,j);
    end
    for j = 1:kCat
        nm = factors(idxCat(j)).name;
        T.(nm) = Ccat(:,j);
    end
end

% ======================================================================
%  7. CONSTRAINTS
% ======================================================================
if ~isempty(cfg.constraints)
    nBefore = height(T);
    keep = true(nBefore,1);
    for ic = 1:numel(cfg.constraints)
        fn = cfg.constraints{ic};
        ok = logical(fn(T));
        if isscalar(ok), ok = repmat(ok, nBefore, 1); end
        keep = keep & ok(:);
    end
    T = T(keep,:);
    fprintf('  Constraints removed %d / %d runs  (%d remain)\n', ...
        nBefore - height(T), nBefore, height(T));
end

% ======================================================================
%  8. REPLICATES
% ======================================================================
nRep = cfg.replicates;
if nRep > 1
    T0 = T;
    T  = repmat(T0, nRep, 1);
    T.Replicate = repelem((1:nRep)', height(T0), 1);
    fprintf('  After %d replicates: %d runs\n', nRep, height(T));
end

% ======================================================================
%  9. BLOCKING
% ======================================================================
nb = cfg.numBlocks;
if nb > 1
    nR = height(T);
    edges = round(linspace(0, nR, nb+1));
    B = zeros(nR,1);
    for b = 1:nb
        B(edges(b)+1 : edges(b+1)) = b;
    end
    T.Block = B;
end

% ======================================================================
% 10. RANDOMISE
% ======================================================================
if cfg.randomize
    rng(cfg.randomSeed, 'twister');
    blocks = unique(T.Block);
    order = [];
    for bi = 1:numel(blocks)
        idx = find(T.Block == blocks(bi));
        order = [order; idx(randperm(numel(idx)))]; %#ok<AGROW>
    end
    T = T(order,:);
end

T.Run = (1:height(T))';

% ======================================================================
% 11. RESPONSE PLACEHOLDERS
% ======================================================================
for r = 1:numel(cfg.responses)
    T.(cfg.responses{r}) = nan(height(T), 1);
end

% ======================================================================
% 12. ASSEMBLE OUTPUTS
% ======================================================================
design.runSheet = T;
design.meta     = cfg.project;

report.title      = cfg.project.title;
report.designType = upper(cfg.designType);
report.nRuns      = height(T);
report.kCont      = kCont;
report.kCat       = kCat;
report.catMethod  = catMethod;
report.blocks     = cfg.numBlocks;
report.replicates = cfg.replicates;
report.randomized = cfg.randomize;

report.summaryText = sprintf([...
    '\n========================================\n' ...
    ' %s - SUMMARY\n' ...
    '========================================\n' ...
    '  Design    : %s\n' ...
    '  Cont.     : %d factors\n' ...
    '  Categ.    : %d factors (%s)\n' ...
    '  Runs      : %d\n' ...
    '  Blocks    : %d\n' ...
    '  Replicates: %d\n' ...
    '  Randomised: %s\n' ...
    '========================================\n'], ...
    cfg.project.title, upper(cfg.designType), ...
    kCont, kCat, catMethod, ...
    height(T), cfg.numBlocks, cfg.replicates, ...
    mat2str(cfg.randomize));

fprintf('%s', report.summaryText);

% ======================================================================
% 13. EXPORT
% ======================================================================
if cfg.export.csv || cfg.export.mat
    if ~exist(cfg.export.folder, 'dir')
        mkdir(cfg.export.folder);
    end
    base = fullfile(cfg.export.folder, cfg.project.name);
    if cfg.export.csv
        csvPath = [base '_runsheet.csv'];
        writetable(T, csvPath);
        fprintf('  Saved CSV : %s\n', csvPath);
    end
    if cfg.export.mat
        matPath = [base '_bundle.mat'];
        save(matPath, 'design','report','cfg');
        fprintf('  Saved MAT : %s\n', matPath);
    end
end

fprintf('  Done.\n\n');

end  % createDOE


% ######################################################################
%  LOCAL / PRIVATE FUNCTIONS
% ######################################################################

function dtype = resolveTypeAlias(raw)
    raw = lower(char(raw));
    switch raw
        case {'fullfact','full','ff','factorial'}
            dtype = 'fullfact';
        case {'fracfact','frac','fractional'}
            dtype = 'fracfact';
        case {'ccd','ccdesign','centralcomposite'}
            dtype = 'ccd';
        case {'bbd','boxbehnken','bb','box-behnken'}
            dtype = 'bbd';
        case {'bbd_ccd','bbdccd','bbd+ccd','hybrid'}
            dtype = 'bbd_ccd';
        case {'lhs','latinhypercube'}
            dtype = 'lhs';
        case {'pb','pbdesign','plackett','plackettburman'}
            dtype = 'pbdesign';
        otherwise
            error('createDOE:badType', ...
                ['Unknown design type "%s".\n' ...
                 'Options: fullfact, fracfact, ccd, bbd, bbd_ccd, lhs, pbdesign.'], raw);
    end
end

function Xc = generateCoded(designType, k, cfg)

    switch designType

        case 'fullfact'
            if k == 0
                Xc = zeros(1, 0);
                return
            end
            Xc = localFullFact(2*ones(1,k));
            Xc = 2*Xc - 3;
            if cfg.centerPoints > 0
                Xc = [Xc; zeros(cfg.centerPoints, k)];
            end

        case 'fracfact'
            if k == 0
                error('createDOE:noFrac', 'Fractional factorial needs >=1 factor.');
            end
            gen = cfg.fracGen;
            if isempty(gen)
                error('createDOE:noGen', ...
                    'Fractional factorial requires ''fracGenerators'' (e.g. ''a b c abc'').');
            end
            Xc = localFracFact(gen);
            if size(Xc,2) ~= k
                error('createDOE:genMismatch', ...
                    'fracGenerators produced %d columns but design needs %d.', ...
                    size(Xc,2), k);
            end
            if cfg.centerPoints > 0
                Xc = [Xc; zeros(cfg.centerPoints, k)];
            end

        case 'ccd'
            if k == 0
                error('createDOE:noCCD', 'CCD needs >=1 factor.');
            end
            alpha = resolveAlpha(cfg.ccdAlpha, k);
            Xfact = localFullFact(2*ones(1,k));
            Xfact = 2*Xfact - 3;
            Xax = zeros(2*k, k);
            for j = 1:k
                Xax(2*j-1, j) = -alpha;
                Xax(2*j,   j) = +alpha;
            end
            Xcp = zeros(max(cfg.centerPoints, 0), k);
            Xc = [Xfact; Xax; Xcp];
            fprintf('  CCD: %d factorial + %d axial + %d center = %d  (alpha=%.4g)\n', ...
                size(Xfact,1), size(Xax,1), size(Xcp,1), size(Xc,1), alpha);

        case 'bbd'
            if k < 3
                error('createDOE:bbdNeed3', ...
                    'Box-Behnken requires >=3 factors (got %d).', k);
            end
            Xc = localBoxBehnken(k);
            if cfg.centerPoints > 0
                Xc = [Xc; zeros(cfg.centerPoints, k)];
            end
            nEdge = size(Xc,1) - 1 - cfg.centerPoints;
            fprintf('  BBD: %d edge-midpoint + 1 built-in center + %d extra center = %d\n', ...
                nEdge, cfg.centerPoints, size(Xc,1));

        case 'bbd_ccd'
            if k < 3
                error('createDOE:bbdccdNeed3', ...
                    'BBD+CCD hybrid requires >=3 factors (got %d).', k);
            end
            alpha = resolveAlpha(cfg.ccdAlpha, k);
            Xbbd = localBoxBehnken(k);
            Xax = zeros(2*k, k);
            for j = 1:k
                Xax(2*j-1, j) = -alpha;
                Xax(2*j,   j) = +alpha;
            end
            Xcp = zeros(max(cfg.centerPoints, 0), k);
            Xc = [Xbbd; Xax; Xcp];
            fprintf('  BBD+CCD: %d BBD + %d axial + %d extra center = %d  (alpha=%.4g)\n', ...
                size(Xbbd,1), size(Xax,1), size(Xcp,1), size(Xc,1), alpha);

        case 'lhs'
            if k == 0
                error('createDOE:noLHS', 'LHS needs >=1 factor.');
            end
            n = cfg.lhsSamples;
            Xc = localLHS(n, k, cfg.lhsMaxIter, cfg.randomSeed);
            Xc = 2*Xc - 1;
            fprintf('  LHS: %d samples, %d factors\n', n, k);

        case 'pbdesign'
            if k == 0
                error('createDOE:noPB', 'Plackett-Burman needs >=1 factor.');
            end
            Xc = localPlackettBurman(k);
            if cfg.centerPoints > 0
                Xc = [Xc; zeros(cfg.centerPoints, k)];
            end
            fprintf('  PB: %d runs for %d factors\n', size(Xc,1), k);

        otherwise
            error('createDOE:internal', 'Unhandled design type "%s".', designType);
    end
end

function alpha = resolveAlpha(spec, k)
    if isnumeric(spec)
        alpha = spec;
        return
    end
    spec = lower(char(spec));
    switch spec
        case {'rotatable','rot'}
            alpha = (2^k)^0.25;
        case {'face','faced'}
            alpha = 1;
        case 'orthogonal'
            alpha = (k*(1+sqrt(1+8*k))/4)^0.25;
        otherwise
            warning('createDOE:alphaDefault', ...
                'Unknown ccdAlpha "%s"; defaulting to rotatable.', spec);
            alpha = (2^k)^0.25;
    end
end

function printHeader(cfg)
    fprintf('\n========================================\n');
    fprintf(' %s\n', upper(cfg.project.title));
    fprintf(' Design type: %s\n', upper(cfg.designType));
    fprintf('========================================\n');
end


% ######################################################################
%  DESIGN GENERATORS
% ######################################################################

function X = localFullFact(levels)
    k = numel(levels);
    n = prod(levels);
    X = zeros(n, k);
    rep = 1;
    for j = 1:k
        nLev = levels(j);
        pat  = repelem((1:nLev)', rep, 1);
        tile = n / (rep * nLev);
        X(:,j) = repmat(pat, tile, 1);
        rep = rep * nLev;
    end
end

function X = localFracFact(generators)
    tokens = strsplit(strtrim(generators));
    nCol   = numel(tokens);
    bases = {};
    for t = 1:nCol
        for ch = tokens{t}
            if ~ismember(ch, bases)
                bases{end+1} = ch; %#ok<AGROW>
            end
        end
    end
    bases = sort(bases);
    nBase = numel(bases);
    idx = localFullFact(2*ones(1, nBase));
    baseMat = 2*idx - 3;
    baseMap = containers.Map(bases, num2cell(1:nBase));
    nRuns = size(baseMat, 1);
    X = zeros(nRuns, nCol);
    for t = 1:nCol
        tok = tokens{t};
        col = ones(nRuns, 1);
        for ch = tok
            col = col .* baseMat(:, baseMap(ch));
        end
        X(:,t) = col;
    end
end

function X = localBoxBehnken(k)
    pairs = nchoosek(1:k, 2);
    nPairs = size(pairs, 1);
    corners = [ -1 -1; -1 1; 1 -1; 1 1 ];
    nPer = size(corners, 1);
    X = zeros(nPairs * nPer, k);
    row = 0;
    for pp = 1:nPairs
        i = pairs(pp,1);
        j = pairs(pp,2);
        for c = 1:nPer
            row = row + 1;
            X(row, i) = corners(c, 1);
            X(row, j) = corners(c, 2);
        end
    end
    X = [X; zeros(1, k)];
end

function X = localLHS(n, k, maxIter, seed)
    rng(seed, 'twister');
    bestX = [];
    bestDist = -Inf;
    for it = 1:maxIter %#ok<FXUP>
        Xt = zeros(n, k);
        for j = 1:k
            perm = randperm(n);
            Xt(:,j) = (perm(:) - rand(n,1)) / n;
        end
        md = minPairDist(Xt);
        if md > bestDist
            bestDist = md;
            bestX = Xt;
        end
    end
    X = bestX;
end

function d = minPairDist(X)
    n = size(X,1);
    d = Inf;
    for i = 1:n-1
        for j = i+1:n
            dij = sqrt(sum((X(i,:)-X(j,:)).^2));
            if dij < d, d = dij; end
        end
    end
end

function X = localPlackettBurman(k)
    if k <= 3
        X = localFullFact(2*ones(1,k));
        X = 2*X - 3;
        return
    end
    n = 4 * ceil((k+1)/4);
    switch n
        case 8,  row1 = [1 1 1 -1 1 -1 -1];
        case 12, row1 = [1 1 -1 1 1 1 -1 -1 -1 1 -1];
        case 16, row1 = [1 1 1 1 -1 1 -1 1 1 -1 -1 1 -1 -1 -1];
        case 20, row1 = [1 1 -1 1 1 -1 -1 -1 -1 1 -1 1 -1 1 1 1 1 -1 -1];
        case 24, row1 = [1 1 1 1 1 -1 1 -1 1 1 -1 -1 1 1 -1 -1 1 -1 1 -1 -1 -1 -1];
        otherwise
            warning('createDOE:pbLarge', ...
                'PB generator not stored for n=%d; using full factorial.', n);
            X = localFullFact(2*ones(1,k));
            X = 2*X - 3;
            return
    end
    ncol = n - 1;
    M = zeros(n, ncol);
    M(1,:) = row1(1:ncol);
    for r = 2:ncol
        M(r,:) = circshift(M(r-1,:), [0 1]);
    end
    M(n,:) = -1;
    X = M(:, 1:k);
end

function C = localAllComb(varargin)
    nArgs = nargin;
    if nArgs == 0
        C = strings(0,0);
        return
    end
    args = varargin;
    for i = 1:nArgs
        args{i} = string(args{i}(:));
    end
    sizes = cellfun(@numel, args);
    nRows = prod(sizes);
    if nRows == 0
        C = strings(0, nArgs);
        return
    end
    C = strings(nRows, nArgs);
    rep_inner = 1;
    for j = 1:nArgs
        nj = sizes(j);
        rep_outer = nRows / (rep_inner * nj);
        idx = repmat(repelem((1:nj)', rep_inner, 1), rep_outer, 1);
        C(:,j) = args{j}(idx);
        rep_inner = rep_inner * nj;
    end
end