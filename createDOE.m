function [design, report, cfg] = createDOE(factors, varargin)
%CREATEDOE  Create a Design of Experiments run-sheet.
%
%   [design, report, cfg] = createDOE(factors, Name, Value, ...)
%
%   Self-contained — no Statistics Toolbox required.
%
%   FACTORS (struct array, one element per factor)
%   ──────────────────────────────────────────────
%     .name   - (string) unique factor name, valid MATLAB identifier
%     .type   - 'continuous' | 'categorical'
%
%     Continuous factors additionally need:
%       .low   - (scalar) lower bound
%       .high  - (scalar) upper bound
%       .units - (string) unit label, e.g. 'degC'  (optional, default '-')
%
%     Categorical factors additionally need:
%       .levels - (string/cell vector) e.g. ["A","B","C"]
%
%   NAME-VALUE OPTIONS
%   ──────────────────
%   Design type
%     'type'           - 'fullfact' | 'fracfact' | 'ccd' | 'bbd' |
%                        'boxbehnken' | 'lhs' | 'pbdesign'
%                        (default 'fullfact')
%
%   General
%     'centerPoints'   - # center points to add           (default 0)
%     'replicates'     - # complete replicates             (default 1)
%     'randomize'      - true/false                        (default true)
%     'randomSeed'     - integer seed                      (default 42)
%     'numBlocks'      - # blocks (sequential split)       (default 1)
%
%   CCD-specific
%     'ccdAlpha'       - 'rotatable' | 'face' | numeric   (default 'rotatable')
%
%   Fractional-factorial
%     'fracGenerators' - generator string for fracfact     (default '')
%
%   LHS-specific
%     'lhsSamples'     - # samples                         (default 50)
%     'lhsMaxIter'     - maximin optimisation iterations   (default 20)
%
%   Plackett-Burman
%     (no extra options — run count picked automatically)
%
%   Constraints
%     'constraints'    - cell array of function handles.
%                        Each receives a table row and returns logical true
%                        to KEEP that row.
%
%   Export / metadata
%     'title'          - display title
%     'projectName'    - used for filenames       (default 'DOE_Project')
%     'owner'          - your name
%     'notes'          - free-text notes
%     'exportFolder'   - path                     (default './DOE_Output')
%     'exportCSV'      - true/false               (default false)
%     'exportMAT'      - true/false               (default false)
%     'responses'      - cell array of response column names
%                        (default {'Response1','Response2'})
%
%   OUTPUTS
%     design.runSheet  - table with Run, Block, Replicate, factor columns,
%                        coded columns (*_coded), and response placeholders
%     report           - struct with summary info
%     cfg              - full configuration (for reproducibility)
%
%   EXAMPLES
%     % --- 1) Two-factor full factorial --------------------------------
%     f = [ struct('name','Temp','type','continuous','low',60,'high',90)
%           struct('name','Time','type','continuous','low',30,'high',180) ];
%     [d,r] = createDOE(f, 'type','fullfact');
%
%     % --- 2) CCD with a categorical factor ----------------------------
%     f = [ struct('name','Mode','type','categorical', ...
%                   'levels',["Co","Mn","Both"])
%           struct('name','Temp','type','continuous','low',60,'high',90)
%           struct('name','Conc','type','continuous','low',0.2,'high',1) ];
%     [d,r] = createDOE(f, 'type','ccd', 'centerPoints',6, ...
%                        'ccdAlpha','rotatable', 'exportCSV',true);
%
%     % --- 3) Box–Behnken with constraints -----------------------------
%     cons = { @(T) ~(T.Temp > 85 & T.Conc > 0.9) };
%     [d,r] = createDOE(f, 'type','bbd', 'constraints',cons);
%
%   See also: the companion QuickReference.m

% ======================================================================
%  0. PARSE INPUTS
% ======================================================================
p = inputParser;
p.FunctionName = 'createDOE';

addRequired(p, 'factors',  @(x) ~isempty(x));

addParameter(p, 'type',           'fullfact');
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
%  1. VALIDATE FACTORS
% ======================================================================
% Accept struct array OR cell array of structs (cell is needed when
% continuous and categorical factors have different fields).
if iscell(factors)
    factors = factors(:);
    tmp = factors;
else
    factors = factors(:);
    tmp = num2cell(factors);
end
nF = numel(tmp);

% Normalise: ensure every factor has all six canonical fields
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
% Rebuild as a struct array with uniform fields (so MATLAB allows indexing)
factors = orderfields(tmp{1}, allFields);
for i = 2:nF
    factors(i) = orderfields(tmp{i}, allFields);
end

% Separate continuous / categorical
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

% Check uniqueness
names = {factors.name};
if numel(unique(names)) ~= nF
    error('createDOE:dupName', 'Factor names must be unique.');
end

% Validate ranges / levels
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
    factors(i).levels = string(f.levels(:));   % normalise to string column
end

% ======================================================================
%  2. BUILD CONFIGURATION STRUCT
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
%  3. GENERATE CODED DESIGN MATRIX  (continuous factors only, in [-1,+1])
% ======================================================================
printHeader(cfg);
fprintf('  Factors: %d continuous, %d categorical\n', kCont, kCat);

Xc = generateCoded(cfg.designType, kCont, cfg);

fprintf('  Base continuous runs: %d\n', size(Xc,1));

% ======================================================================
%  4. GENERATE CATEGORICAL COMBINATIONS
% ======================================================================
if kCat > 0
    catLevels = cell(1, kCat);
    for j = 1:kCat
        catLevels{j} = factors(idxCat(j)).levels;
    end
    Ccat = localAllComb(catLevels{:});       % nCat-combos × kCat
    fprintf('  Categorical combinations: %d\n', size(Ccat,1));
else
    Ccat = strings(0,0);
end

% ======================================================================
%  5. CROSS continuous × categorical
% ======================================================================
if kCont > 0 && kCat > 0
    nC = size(Xc,1);
    nL = size(Ccat,1);
    % Every categorical combo gets its own copy of the continuous design
    Xc   = repmat(Xc, nL, 1);                 % (nC*nL) × kCont
    Ccat = repelem(Ccat, nC, 1);               % (nC*nL) × kCat
    fprintf('  After crossing: %d runs\n', size(Xc,1));
elseif kCont == 0 && kCat > 0
    Xc = zeros(size(Ccat,1), 0);              % no continuous columns
elseif kCont > 0 && kCat == 0
    Ccat = strings(size(Xc,1), 0);            % no categorical columns
else
    error('createDOE:noFactors', 'Need at least one factor.');
end

nBase = size(Xc, 1);    % rows after crossing, before replicates

% ======================================================================
%  6. DECODE to actual values
% ======================================================================
Xa = codedToActual(Xc, factors, idxCont);

% ======================================================================
%  7. BUILD TABLE
% ======================================================================
T = table();
T.Run = (1:nBase)';
T.Block     = ones(nBase,1);
T.Replicate = ones(nBase,1);

% Continuous columns (actual + coded)
for j = 1:kCont
    nm = factors(idxCont(j)).name;
    T.(nm) = Xa(:,j);
    T.([nm '_coded']) = Xc(:,j);
end

% Categorical columns
for j = 1:kCat
    nm = factors(idxCat(j)).name;
    T.(nm) = Ccat(:,j);
end

% ======================================================================
%  8. CONSTRAINTS
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
%  9. REPLICATES
% ======================================================================
nRep = cfg.replicates;
if nRep > 1
    T0 = T;
    T  = repmat(T0, nRep, 1);
    T.Replicate = repelem((1:nRep)', height(T0), 1);
    fprintf('  After %d replicates: %d runs\n', nRep, height(T));
end

% ======================================================================
% 10. BLOCKING
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
% 11. RANDOMISE
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

% Renumber
T.Run = (1:height(T))';

% ======================================================================
% 12. RESPONSE PLACEHOLDERS
% ======================================================================
for r = 1:numel(cfg.responses)
    T.(cfg.responses{r}) = nan(height(T), 1);
end

% ======================================================================
% 13. ASSEMBLE OUTPUTS
% ======================================================================
design.runSheet = T;
design.meta     = cfg.project;

report.title      = cfg.project.title;
report.designType = upper(cfg.designType);
report.nRuns      = height(T);
report.kCont      = kCont;
report.kCat       = kCat;
report.blocks     = cfg.numBlocks;
report.replicates = cfg.replicates;
report.randomized = cfg.randomize;

report.summaryText = sprintf([...
    '\n========================================\n' ...
    ' %s — SUMMARY\n' ...
    '========================================\n' ...
    '  Design   : %s\n' ...
    '  Cont.    : %d factors\n' ...
    '  Categ.   : %d factors\n' ...
    '  Runs     : %d\n' ...
    '  Blocks   : %d\n' ...
    '  Replicates: %d\n' ...
    '  Randomised: %s\n' ...
    '========================================\n'], ...
    cfg.project.title, upper(cfg.designType), ...
    kCont, kCat, height(T), cfg.numBlocks, cfg.replicates, ...
    mat2str(cfg.randomize));

fprintf('%s', report.summaryText);

% ======================================================================
% 14. EXPORT
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

% ======================================================================
function dtype = resolveTypeAlias(raw)
%RESOLVETYPEALIAS  Canonicalise design-type string.
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
        case {'lhs','latinhypercube'}
            dtype = 'lhs';
        case {'pb','pbdesign','plackett','plackettburman'}
            dtype = 'pbdesign';
        otherwise
            error('createDOE:badType', ...
                'Unknown design type "%s".  Options: fullfact, fracfact, ccd, bbd, lhs, pbdesign.', raw);
    end
end

% ======================================================================
function Xc = generateCoded(designType, k, cfg)
%GENERATECODED  Return n×k coded matrix in [-1,+1] (or ±alpha for CCD).

    switch designType

        % ─────────────────────────────────────────────────────────────
        case 'fullfact'
            if k == 0
                Xc = zeros(1, 0);
                return
            end
            Xc = localFullFact(2*ones(1,k));        % rows of {1,2}
            Xc = 2*Xc - 3;                          % map to {-1,+1}
            if cfg.centerPoints > 0
                Xc = [Xc; zeros(cfg.centerPoints, k)];
            end

        % ─────────────────────────────────────────────────────────────
        case 'fracfact'
            if k == 0
                error('createDOE:noFrac', 'Fractional factorial needs >=1 continuous factor.');
            end
            gen = cfg.fracGen;
            if isempty(gen)
                error('createDOE:noGen', ...
                    'Fractional factorial requires ''fracGenerators'' (e.g. ''a b c abc'').');
            end
            Xc = localFracFact(gen);
            % Verify column count matches kCont
            if size(Xc,2) ~= k
                error('createDOE:genMismatch', ...
                    'fracGenerators produced %d columns but there are %d continuous factors.', ...
                    size(Xc,2), k);
            end
            if cfg.centerPoints > 0
                Xc = [Xc; zeros(cfg.centerPoints, k)];
            end

        % ─────────────────────────────────────────────────────────────
        case 'ccd'
            if k == 0
                error('createDOE:noCCD', 'CCD needs >=1 continuous factor.');
            end
            alpha = resolveAlpha(cfg.ccdAlpha, k);

            % Factorial portion: 2^k full factorial at ±1
            Xfact = localFullFact(2*ones(1,k));
            Xfact = 2*Xfact - 3;                   % {-1,+1}

            % Axial (star) points: ±alpha on each axis
            Xax = zeros(2*k, k);
            for j = 1:k
                Xax(2*j-1, j) = -alpha;
                Xax(2*j,   j) = +alpha;
            end

            % Center points
            Xcp = zeros(max(cfg.centerPoints, 0), k);

            Xc = [Xfact; Xax; Xcp];
            fprintf('  CCD: %d factorial + %d axial + %d center = %d  (alpha=%.4g)\n', ...
                size(Xfact,1), size(Xax,1), size(Xcp,1), size(Xc,1), alpha);

        % ─────────────────────────────────────────────────────────────
        case 'bbd'
            if k < 3
                error('createDOE:bbdNeed3', ...
                    'Box–Behnken requires >=3 continuous factors (got %d).', k);
            end
            Xc = localBoxBehnken(k);
            if cfg.centerPoints > 0
                Xc = [Xc; zeros(cfg.centerPoints, k)];
            end
            fprintf('  BBD: %d edge-midpoint + %d center = %d\n', ...
                size(Xc,1)-cfg.centerPoints, cfg.centerPoints, size(Xc,1));

        % ─────────────────────────────────────────────────────────────
        case 'lhs'
            if k == 0
                error('createDOE:noLHS', 'LHS needs >=1 continuous factor.');
            end
            n = cfg.lhsSamples;
            Xc = localLHS(n, k, cfg.lhsMaxIter, cfg.randomSeed);
            Xc = 2*Xc - 1;                         % [0,1]→[-1,+1]
            fprintf('  LHS: %d samples, %d factors\n', n, k);

        % ─────────────────────────────────────────────────────────────
        case 'pbdesign'
            if k == 0
                error('createDOE:noPB', 'Plackett–Burman needs >=1 continuous factor.');
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

% ======================================================================
function alpha = resolveAlpha(spec, k)
%RESOLVEALPHA  Convert alpha specification to a numeric value.
    if isnumeric(spec)
        alpha = spec;
        return
    end
    spec = lower(char(spec));
    switch spec
        case {'rotatable','rot'}
            alpha = (2^k)^0.25;     % = 2^(k/4)
        case {'face','faced'}
            alpha = 1;
        case 'orthogonal'
            % orthogonal alpha for a CCD
            nf = 2^k;
            alpha = ( (sqrt(nf+2*k+0) - sqrt(nf)) / sqrt(2) );
            % Simplified common formula:
            alpha = sqrt( k * (sqrt(nf+2*k) - sqrt(nf)) / sqrt(2) );
            % Fallback to a reliable formula
            alpha = (k*(1+sqrt(1+8*k))/4)^0.25;
        otherwise
            warning('createDOE:alphaDefault', ...
                'Unknown ccdAlpha "%s"; defaulting to rotatable.', spec);
            alpha = (2^k)^0.25;
    end
end

% ======================================================================
function Xa = codedToActual(Xc, factors, idxCont)
%CODEDTOACTUAL  Map coded values to engineering units.
%   Mapping: coded -1 → low,  +1 → high  (linear).
%   CCD star points at ±alpha map proportionally outside [low, high].
    k  = numel(idxCont);
    n  = size(Xc,1);
    Xa = zeros(n, k);
    for j = 1:k
        f   = factors(idxCont(j));
        mid = (f.high + f.low)  / 2;
        hwd = (f.high - f.low) / 2;
        Xa(:,j) = mid + Xc(:,j) * hwd;
    end
end

% ======================================================================
function printHeader(cfg)
    fprintf('\n========================================\n');
    fprintf(' %s\n', upper(cfg.project.title));
    fprintf(' Design type: %s\n', upper(cfg.designType));
    fprintf('========================================\n');
end


% ######################################################################
%  DESIGN GENERATORS  (no toolbox needed)
% ######################################################################

% ======================================================================
function X = localFullFact(levels)
%LOCALFULLFACT  Full factorial: every combination of 1:levels(j).
%   X = localFullFact([2 3]) returns a 6×2 matrix.
    k = numel(levels);
    n = prod(levels);
    X = zeros(n, k);
    rep = 1;
    for j = 1:k
        nLev = levels(j);
        pat  = repelem((1:nLev)', rep, 1);       % repeat each element
        tile = n / (rep * nLev);                   % repeat the pattern
        X(:,j) = repmat(pat, tile, 1);
        rep = rep * nLev;
    end
end

% ======================================================================
function X = localFracFact(generators)
%LOCALFRACFACT  Two-level fractional factorial from a generator string.
%   generators = 'a b c abc' means 4 columns; the 4th is the product of
%   columns a, b, c.  Base factors a-z give independent ±1 columns.
%   Products like 'ab', 'acd' multiply the corresponding base columns.

    tokens = strsplit(strtrim(generators));
    nCol   = numel(tokens);

    % Identify base factors (single-letter tokens)
    bases = {};
    for t = 1:nCol
        if numel(tokens{t}) == 1
            if ~ismember(tokens{t}, bases)
                bases{end+1} = tokens{t}; %#ok<AGROW>
            end
        else
            for ch = tokens{t}
                if ~ismember(ch, bases)
                    bases{end+1} = ch; %#ok<AGROW>
                end
            end
        end
    end
    bases = sort(bases);
    nBase = numel(bases);

    % Full factorial on base factors: rows in {-1,+1}
    idx = localFullFact(2*ones(1, nBase));
    baseMat = 2*idx - 3;                            % {-1,+1}

    % Map base letter → column index
    baseMap = containers.Map(bases, num2cell(1:nBase));

    % Build output columns
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

% ======================================================================
function X = localBoxBehnken(k)
%LOCALBOXBEHNKEN  Box–Behnken design for k>=3 continuous factors.
%   Each run sets exactly 2 factors to ±1 and the rest to 0, plus a
%   center-point row (at least 1 is always included).

    pairs = nchoosek(1:k, 2);
    nPairs = size(pairs, 1);
    % For each pair, 4 runs: (±1, ±1)
    corners = [ -1 -1; -1 1; 1 -1; 1 1 ];
    nPer = size(corners, 1);
    X = zeros(nPairs * nPer, k);
    row = 0;
    for p = 1:nPairs
        i = pairs(p,1);
        j = pairs(p,2);
        for c = 1:nPer
            row = row + 1;
            X(row, i) = corners(c, 1);
            X(row, j) = corners(c, 2);
        end
    end
    % Always include at least one center point so the design is non-singular
    X = [X; zeros(1, k)];
end

% ======================================================================
function X = localLHS(n, k, maxIter, seed)
%LOCALLHS  Maximin Latin Hypercube Sampling on [0,1]^k.
%   Generates n points in k dimensions.  Optimises for maximin distance
%   over maxIter random starts.
    rng(seed, 'twister');
    bestX = [];
    bestDist = -Inf;
    for it = 1:maxIter
        % Random LHS: one sample per stratum in each dimension
        Xt = zeros(n, k);
        for j = 1:k
            perm = randperm(n);
            Xt(:,j) = (perm(:) - rand(n,1)) / n;
        end
        % Minimum pairwise distance
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

% ======================================================================
function X = localPlackettBurman(k)
%LOCALPLACKETTBURMAN  Plackett–Burman screening design.
%   Produces a design in n runs (n = smallest multiple of 4 >= k+1)
%   with columns in {-1,+1}.  We use the standard first-row generators
%   for n = 8, 12, 16, 20, 24.  Falls back to full factorial if k<=3.

    if k <= 3
        X = localFullFact(2*ones(1,k));
        X = 2*X - 3;
        return
    end

    % n must be a multiple of 4, >= k+1
    n = 4 * ceil((k+1)/4);

    % Known first-row generators (excluding the leading +1)
    switch n
        case 8
            row1 = [1 1 1 -1 1 -1 -1];
        case 12
            row1 = [1 1 -1 1 1 1 -1 -1 -1 1 -1];
        case 16
            row1 = [1 1 1 1 -1 1 -1 1 1 -1 -1 1 -1 -1 -1];
        case 20
            row1 = [1 1 -1 1 1 -1 -1 -1 -1 1 -1 1 -1 1 1 1 1 -1 -1];
        case 24
            row1 = [1 1 1 1 1 -1 1 -1 1 1 -1 -1 1 1 -1 -1 1 -1 1 -1 -1 -1 -1];
        otherwise
            % For very large k, fall back to a Hadamard-based approach
            % or just use a 2-level full factorial (expensive but correct)
            warning('createDOE:pbLarge', ...
                'Plackett–Burman generator not stored for n=%d; using full factorial.', n);
            X = localFullFact(2*ones(1,k));
            X = 2*X - 3;
            return
    end

    % Build the cyclic matrix  (n rows × n-1 columns)
    ncol = n - 1;
    M = zeros(n, ncol);
    M(1,:) = row1(1:ncol);
    for r = 2:ncol
        M(r,:) = circshift(M(r-1,:), [0 1]);
    end
    M(n,:) = -1;                                    % last row is all -1

    % Keep only the first k columns
    X = M(:, 1:k);
end

% ======================================================================
function C = localAllComb(varargin)
%LOCALALLCOMB  Cartesian product of string (or numeric) vectors.
%   C = localAllComb(v1, v2, ..., vn)  returns an (prod of lengths) × n
%   array with every combination of elements.

    nArgs = nargin;
    if nArgs == 0
        C = strings(0,0);
        return
    end

    % Ensure each input is a column vector of strings
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

    % Build indices via repeated tiling / element-repetition
    rep_inner = 1;                 % how many times to repeat each element
    for j = 1:nArgs
        nj    = sizes(j);
        rep_outer = nRows / (rep_inner * nj);
        idx = repmat(repelem((1:nj)', rep_inner, 1), rep_outer, 1);
        C(:,j) = args{j}(idx);
        rep_inner = rep_inner * nj;
    end
end
