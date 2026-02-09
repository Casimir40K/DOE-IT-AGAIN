%% ========================================================================
%  createDOE — COMPREHENSIVE TEST SUITE
%  ========================================================================
%  Run this to verify every design type and edge case.
%  Each test prints PASS/FAIL.  Errors are caught and reported.
%  ========================================================================
clear; clc; close all;
nPass = 0;  nFail = 0;

%% Helper: check a design
check = @(d,r, minRuns, label) localCheck(d,r,minRuns,label);

% ─────────────────────────────────────────────────────────────────────────
%% TEST 1: Two-factor full factorial (simplest case)
% ─────────────────────────────────────────────────────────────────────────
try
    f = [ struct('name','Temp','type','continuous','low',60,'high',90)
          struct('name','Time','type','continuous','low',30,'high',180) ];
    [d,r] = createDOE(f, 'type','fullfact', 'randomize',false);
    assert(r.nRuns == 4, 'Expected 4 runs for 2^2');
    assert(all(ismember(d.runSheet.Temp, [60 90])), 'Temp should be at bounds');
    fprintf('TEST 1  PASS  (2-factor fullfact: %d runs)\n', r.nRuns);
    nPass = nPass+1;
catch ME
    fprintf('TEST 1  FAIL  %s\n', ME.message); nFail = nFail+1;
end

% ─────────────────────────────────────────────────────────────────────────
%% TEST 2: Full factorial + center points
% ─────────────────────────────────────────────────────────────────────────
try
    f = [ struct('name','A','type','continuous','low',0,'high',10)
          struct('name','B','type','continuous','low',0,'high',10) ];
    [d,r] = createDOE(f, 'type','fullfact', 'centerPoints',3, 'randomize',false);
    assert(r.nRuns == 7, 'Expected 4+3=7 runs');
    % Center coded values should be 0
    coded_A = d.runSheet.A_coded;
    assert(sum(coded_A == 0) == 3, 'Should have 3 center points at coded=0');
    fprintf('TEST 2  PASS  (fullfact + 3 center: %d runs)\n', r.nRuns);
    nPass = nPass+1;
catch ME
    fprintf('TEST 2  FAIL  %s\n', ME.message); nFail = nFail+1;
end

% ─────────────────────────────────────────────────────────────────────────
%% TEST 3: Full factorial with 1 categorical × 2 continuous
% ─────────────────────────────────────────────────────────────────────────
try
    f = { struct('name','Cat','type','categorical','levels',["X","Y","Z"])
          struct('name','A','type','continuous','low',1,'high',5)
          struct('name','B','type','continuous','low',10,'high',20) };
    [d,r] = createDOE(f, 'type','fullfact', 'randomize',false);
    % 2^2 continuous × 3 categorical = 12
    assert(r.nRuns == 12, 'Expected 4*3=12 runs, got %d', r.nRuns);
    assert(all(ismember(d.runSheet.Cat, ["X","Y","Z"])), 'Cat levels wrong');
    fprintf('TEST 3  PASS  (mixed cat+cont fullfact: %d runs)\n', r.nRuns);
    nPass = nPass+1;
catch ME
    fprintf('TEST 3  FAIL  %s\n', ME.message); nFail = nFail+1;
end

% ─────────────────────────────────────────────────────────────────────────
%% TEST 4: CCD with 3 continuous factors
% ─────────────────────────────────────────────────────────────────────────
try
    f = [ struct('name','T','type','continuous','low',60,'high',90)
          struct('name','P','type','continuous','low',1,'high',10)
          struct('name','C','type','continuous','low',0.1,'high',1) ];
    [d,r] = createDOE(f, 'type','ccd', 'centerPoints',6, ...
                       'ccdAlpha','rotatable', 'randomize',false);
    % 2^3=8 factorial + 2*3=6 axial + 6 center = 20
    assert(r.nRuns == 20, 'Expected 20 CCD runs, got %d', r.nRuns);
    % Check that axial points exist (coded values > 1 or < -1)
    maxCoded = max(abs(d.runSheet.T_coded));
    assert(maxCoded > 1, 'CCD should have axial points outside [-1,1]');
    fprintf('TEST 4  PASS  (CCD 3-factor: %d runs, alpha=%.3f)\n', r.nRuns, maxCoded);
    nPass = nPass+1;
catch ME
    fprintf('TEST 4  FAIL  %s\n', ME.message); nFail = nFail+1;
end

% ─────────────────────────────────────────────────────────────────────────
%% TEST 5: CCD face-centered (alpha=1)
% ─────────────────────────────────────────────────────────────────────────
try
    f = [ struct('name','A','type','continuous','low',0,'high',100)
          struct('name','B','type','continuous','low',0,'high',100) ];
    [d,r] = createDOE(f, 'type','ccd', 'centerPoints',3, ...
                       'ccdAlpha','face', 'randomize',false);
    % 2^2=4 + 2*2=4 + 3 = 11
    assert(r.nRuns == 11, 'Expected 11 runs, got %d', r.nRuns);
    maxCoded = max(abs(d.runSheet.A_coded));
    assert(abs(maxCoded - 1) < 1e-10, 'Face-centered alpha should be 1');
    fprintf('TEST 5  PASS  (CCD face-centered: %d runs)\n', r.nRuns);
    nPass = nPass+1;
catch ME
    fprintf('TEST 5  FAIL  %s\n', ME.message); nFail = nFail+1;
end

% ─────────────────────────────────────────────────────────────────────────
%% TEST 6: Box–Behnken with 3 factors
% ─────────────────────────────────────────────────────────────────────────
try
    f = [ struct('name','A','type','continuous','low',10,'high',30)
          struct('name','B','type','continuous','low',5,'high',15)
          struct('name','C','type','continuous','low',100,'high',200) ];
    [d,r] = createDOE(f, 'type','bbd', 'centerPoints',5, 'randomize',false);
    % 3 pairs × 4 corners + 1 built-in center + 5 extra = 18
    assert(r.nRuns == 18, 'Expected 18 BBD runs, got %d', r.nRuns);
    % No corner points in BBD — no row should have all factors at ±1
    coded = [d.runSheet.A_coded, d.runSheet.B_coded, d.runSheet.C_coded];
    allCorner = all(abs(coded) == 1, 2);
    assert(~any(allCorner), 'BBD should have no corner points');
    fprintf('TEST 6  PASS  (BBD 3-factor: %d runs)\n', r.nRuns);
    nPass = nPass+1;
catch ME
    fprintf('TEST 6  FAIL  %s\n', ME.message); nFail = nFail+1;
end

% ─────────────────────────────────────────────────────────────────────────
%% TEST 7: Box–Behnken with 4 factors
% ─────────────────────────────────────────────────────────────────────────
try
    f = [ struct('name','A','type','continuous','low',0,'high',1)
          struct('name','B','type','continuous','low',0,'high',1)
          struct('name','C','type','continuous','low',0,'high',1)
          struct('name','D','type','continuous','low',0,'high',1) ];
    [d,r] = createDOE(f, 'type','bbd', 'centerPoints',3, 'randomize',false);
    % 6 pairs × 4 + 1 built-in center + 3 extra = 28
    assert(r.nRuns == 28, 'Expected 28 BBD runs, got %d', r.nRuns);
    fprintf('TEST 7  PASS  (BBD 4-factor: %d runs)\n', r.nRuns);
    nPass = nPass+1;
catch ME
    fprintf('TEST 7  FAIL  %s\n', ME.message); nFail = nFail+1;
end

% ─────────────────────────────────────────────────────────────────────────
%% TEST 8: Fractional factorial (2^(4-1))
% ─────────────────────────────────────────────────────────────────────────
try
    f = [ struct('name','A','type','continuous','low',0,'high',1)
          struct('name','B','type','continuous','low',0,'high',1)
          struct('name','C','type','continuous','low',0,'high',1)
          struct('name','D','type','continuous','low',0,'high',1) ];
    [d,r] = createDOE(f, 'type','fracfact', ...
                       'fracGenerators','a b c abc', 'randomize',false);
    assert(r.nRuns == 8, 'Expected 8 runs for 2^3 base, got %d', r.nRuns);
    % D column should be product of A,B,C (the generator)
    coded = [d.runSheet.A_coded, d.runSheet.B_coded, ...
             d.runSheet.C_coded, d.runSheet.D_coded];
    product = coded(:,1) .* coded(:,2) .* coded(:,3);
    assert(all(abs(coded(:,4) - product) < 1e-10), 'D should equal A*B*C');
    fprintf('TEST 8  PASS  (fracfact 2^(4-1): %d runs)\n', r.nRuns);
    nPass = nPass+1;
catch ME
    fprintf('TEST 8  FAIL  %s\n', ME.message); nFail = nFail+1;
end

% ─────────────────────────────────────────────────────────────────────────
%% TEST 9: LHS (50 samples, 3 factors)
% ─────────────────────────────────────────────────────────────────────────
try
    f = [ struct('name','X1','type','continuous','low',0,'high',1)
          struct('name','X2','type','continuous','low',0,'high',1)
          struct('name','X3','type','continuous','low',0,'high',1) ];
    [d,r] = createDOE(f, 'type','lhs', 'lhsSamples',50, 'randomize',false);
    assert(r.nRuns == 50, 'Expected 50 LHS runs, got %d', r.nRuns);
    % All actual values in [0,1]
    vals = [d.runSheet.X1, d.runSheet.X2, d.runSheet.X3];
    assert(all(vals(:) >= 0 & vals(:) <= 1), 'LHS values out of [0,1]');
    fprintf('TEST 9  PASS  (LHS 50×3: %d runs)\n', r.nRuns);
    nPass = nPass+1;
catch ME
    fprintf('TEST 9  FAIL  %s\n', ME.message); nFail = nFail+1;
end

% ─────────────────────────────────────────────────────────────────────────
%% TEST 10: Plackett–Burman (7 factors → 8 runs)
% ─────────────────────────────────────────────────────────────────────────
try
    f = [];
    for i = 1:7
        f = [f; struct('name',sprintf('F%d',i),'type','continuous', ...
                        'low',0,'high',10)]; %#ok<AGROW>
    end
    [d,r] = createDOE(f, 'type','pbdesign', 'randomize',false);
    assert(r.nRuns == 8, 'Expected 8 PB runs for 7 factors, got %d', r.nRuns);
    % All coded values should be ±1
    for i = 1:7
        v = d.runSheet.(sprintf('F%d_coded',i));
        assert(all(abs(v)==1), 'PB coded values should be +-1');
    end
    fprintf('TEST 10 PASS  (Plackett-Burman 7 factors: %d runs)\n', r.nRuns);
    nPass = nPass+1;
catch ME
    fprintf('TEST 10 FAIL  %s\n', ME.message); nFail = nFail+1;
end

% ─────────────────────────────────────────────────────────────────────────
%% TEST 11: Categorical-only design (no continuous factors)
% ─────────────────────────────────────────────────────────────────────────
try
    f = [ struct('name','Cat1','type','categorical','levels',["A","B"])
          struct('name','Cat2','type','categorical','levels',["X","Y","Z"]) ];
    [d,r] = createDOE(f, 'type','fullfact', 'randomize',false);
    assert(r.nRuns == 6, 'Expected 2*3=6 runs, got %d', r.nRuns);
    assert(all(ismember(d.runSheet.Cat1, ["A","B"])), 'Cat1 levels wrong');
    fprintf('TEST 11 PASS  (categorical-only: %d runs)\n', r.nRuns);
    nPass = nPass+1;
catch ME
    fprintf('TEST 11 FAIL  %s\n', ME.message); nFail = nFail+1;
end

% ─────────────────────────────────────────────────────────────────────────
%% TEST 12: Replicates
% ─────────────────────────────────────────────────────────────────────────
try
    f = [ struct('name','A','type','continuous','low',0,'high',1)
          struct('name','B','type','continuous','low',0,'high',1) ];
    [d,r] = createDOE(f, 'type','fullfact', 'replicates',3, 'randomize',false);
    assert(r.nRuns == 12, 'Expected 4*3=12 runs, got %d', r.nRuns);
    assert(max(d.runSheet.Replicate) == 3, 'Should have 3 replicate levels');
    fprintf('TEST 12 PASS  (3 replicates: %d runs)\n', r.nRuns);
    nPass = nPass+1;
catch ME
    fprintf('TEST 12 FAIL  %s\n', ME.message); nFail = nFail+1;
end

% ─────────────────────────────────────────────────────────────────────────
%% TEST 13: Blocking
% ─────────────────────────────────────────────────────────────────────────
try
    f = [ struct('name','A','type','continuous','low',0,'high',1)
          struct('name','B','type','continuous','low',0,'high',1) ];
    [d,r] = createDOE(f, 'type','fullfact', 'numBlocks',2, 'randomize',false);
    assert(numel(unique(d.runSheet.Block)) == 2, 'Should have 2 blocks');
    fprintf('TEST 13 PASS  (2 blocks)\n');
    nPass = nPass+1;
catch ME
    fprintf('TEST 13 FAIL  %s\n', ME.message); nFail = nFail+1;
end

% ─────────────────────────────────────────────────────────────────────────
%% TEST 14: Constraints
% ─────────────────────────────────────────────────────────────────────────
try
    f = [ struct('name','Temp','type','continuous','low',60,'high',100)
          struct('name','Press','type','continuous','low',1,'high',10) ];
    cons = { @(T) ~(T.Temp > 90 & T.Press > 8) };
    [d,r] = createDOE(f, 'type','fullfact', 'constraints',cons, 'randomize',false);
    % Original: 4 runs.  The (100,10) corner should be removed → 3 runs
    assert(r.nRuns == 3, 'Expected 3 runs after constraint, got %d', r.nRuns);
    fprintf('TEST 14 PASS  (constraint removed 1 run: %d remain)\n', r.nRuns);
    nPass = nPass+1;
catch ME
    fprintf('TEST 14 FAIL  %s\n', ME.message); nFail = nFail+1;
end

% ─────────────────────────────────────────────────────────────────────────
%% TEST 15: CCD + categorical (the RUNSCRIPT use case)
% ─────────────────────────────────────────────────────────────────────────
try
    f = { struct('name','Mode','type','categorical', ...
                 'levels',["Co-only","Mn-only","Together"])
          struct('name','Temperature','type','continuous','low',60,'high',90)
          struct('name','Time','type','continuous','low',30,'high',180)
          struct('name','Concentration','type','continuous','low',0.2,'high',1.0)
          struct('name','SLRatio','type','continuous','low',0.02,'high',0.2) };
    [d,r] = createDOE(f, 'type','ccd', 'centerPoints',6, ...
                       'ccdAlpha','rotatable', 'randomize',true, ...
                       'randomSeed',42);
    % 2^4=16 + 2*4=8 + 6 = 30 per category × 3 = 90
    assert(r.nRuns == 90, 'Expected 90 runs, got %d', r.nRuns);
    assert(all(ismember(d.runSheet.Mode, ["Co-only","Mn-only","Together"])));
    fprintf('TEST 15 PASS  (CCD+categorical lithium formate: %d runs)\n', r.nRuns);
    nPass = nPass+1;
catch ME
    fprintf('TEST 15 FAIL  %s\n', ME.message); nFail = nFail+1;
end

% ─────────────────────────────────────────────────────────────────────────
%% TEST 16: BBD + categorical (the RUNSCRIPT 'bbd' alias)
% ─────────────────────────────────────────────────────────────────────────
try
    f = { struct('name','Mode','type','categorical', ...
                 'levels',["Co-only","Mn-only","Together"])
          struct('name','Temperature','type','continuous','low',60,'high',90)
          struct('name','Time','type','continuous','low',30,'high',180)
          struct('name','Concentration','type','continuous','low',0.2,'high',1.0)
          struct('name','SLRatio','type','continuous','low',0.02,'high',0.2) };
    [d,r] = createDOE(f, 'type','bbd', 'centerPoints',6, 'randomize',false);
    % BBD(4): 6 pairs × 4 + 1 built-in center + 6 = 31 per cat × 3 = 93
    assert(r.nRuns == 93, 'Expected 93 runs, got %d', r.nRuns);
    fprintf('TEST 16 PASS  (BBD+categorical: %d runs)\n', r.nRuns);
    nPass = nPass+1;
catch ME
    fprintf('TEST 16 FAIL  %s\n', ME.message); nFail = nFail+1;
end

% ─────────────────────────────────────────────────────────────────────────
%% TEST 17: Randomisation is reproducible
% ─────────────────────────────────────────────────────────────────────────
try
    f = [ struct('name','A','type','continuous','low',0,'high',1)
          struct('name','B','type','continuous','low',0,'high',1) ];
    [d1,~] = createDOE(f, 'type','fullfact', 'randomize',true, 'randomSeed',99);
    [d2,~] = createDOE(f, 'type','fullfact', 'randomize',true, 'randomSeed',99);
    assert(isequal(d1.runSheet, d2.runSheet), 'Same seed should give same order');
    fprintf('TEST 17 PASS  (reproducible randomisation)\n');
    nPass = nPass+1;
catch ME
    fprintf('TEST 17 FAIL  %s\n', ME.message); nFail = nFail+1;
end

% ─────────────────────────────────────────────────────────────────────────
%% TEST 18: Custom response column names
% ─────────────────────────────────────────────────────────────────────────
try
    f = [ struct('name','A','type','continuous','low',0,'high',1) ];
    [d,~] = createDOE(f, 'type','fullfact', ...
                       'responses',{'Yield','Purity','Color'});
    assert(ismember('Yield',  d.runSheet.Properties.VariableNames));
    assert(ismember('Purity', d.runSheet.Properties.VariableNames));
    assert(ismember('Color',  d.runSheet.Properties.VariableNames));
    fprintf('TEST 18 PASS  (custom response columns)\n');
    nPass = nPass+1;
catch ME
    fprintf('TEST 18 FAIL  %s\n', ME.message); nFail = nFail+1;
end

% ─────────────────────────────────────────────────────────────────────────
%% TEST 19: Single continuous factor CCD
% ─────────────────────────────────────────────────────────────────────────
try
    f = struct('name','X','type','continuous','low',0,'high',10);
    [d,r] = createDOE(f, 'type','ccd', 'centerPoints',3, 'randomize',false);
    % 2^1=2 factorial + 2*1=2 axial + 3 center = 7
    assert(r.nRuns == 7, 'Expected 7 runs, got %d', r.nRuns);
    fprintf('TEST 19 PASS  (1-factor CCD: %d runs)\n', r.nRuns);
    nPass = nPass+1;
catch ME
    fprintf('TEST 19 FAIL  %s\n', ME.message); nFail = nFail+1;
end

% ─────────────────────────────────────────────────────────────────────────
%% TEST 20: Type aliases ('ff', 'bb', 'box-behnken', etc.)
% ─────────────────────────────────────────────────────────────────────────
try
    f = [ struct('name','A','type','continuous','low',0,'high',1)
          struct('name','B','type','continuous','low',0,'high',1) ];
    [~,r1] = createDOE(f, 'type','ff',        'randomize',false);
    [~,r2] = createDOE(f, 'type','fullfact',   'randomize',false);
    [~,r3] = createDOE(f, 'type','factorial',   'randomize',false);
    assert(r1.nRuns == r2.nRuns && r2.nRuns == r3.nRuns, 'Aliases should match');
    fprintf('TEST 20 PASS  (type aliases work)\n');
    nPass = nPass+1;
catch ME
    fprintf('TEST 20 FAIL  %s\n', ME.message); nFail = nFail+1;
end

% ─────────────────────────────────────────────────────────────────────────
%% TEST 21: Error on BBD with <3 factors
% ─────────────────────────────────────────────────────────────────────────
try
    f = [ struct('name','A','type','continuous','low',0,'high',1)
          struct('name','B','type','continuous','low',0,'high',1) ];
    threw = false;
    try
        createDOE(f, 'type','bbd');
    catch
        threw = true;
    end
    assert(threw, 'Should have thrown error for BBD with 2 factors');
    fprintf('TEST 21 PASS  (BBD <3 factors rejected)\n');
    nPass = nPass+1;
catch ME
    fprintf('TEST 21 FAIL  %s\n', ME.message); nFail = nFail+1;
end

% ─────────────────────────────────────────────────────────────────────────
%% TEST 22: Actual values map correctly
% ─────────────────────────────────────────────────────────────────────────
try
    f = struct('name','X','type','continuous','low',100,'high',200);
    [d,~] = createDOE(f, 'type','fullfact', 'randomize',false);
    coded = d.runSheet.X_coded;
    actual = d.runSheet.X;
    % coded=-1 → actual=100, coded=+1 → actual=200
    for i = 1:height(d.runSheet)
        expected = 150 + coded(i)*50;
        assert(abs(actual(i) - expected) < 1e-10, 'Actual value mismatch');
    end
    fprintf('TEST 22 PASS  (coded→actual mapping correct)\n');
    nPass = nPass+1;
catch ME
    fprintf('TEST 22 FAIL  %s\n', ME.message); nFail = nFail+1;
end


%% ====================================================================
%  SUMMARY
% =====================================================================
fprintf('\n========================================\n');
fprintf('  RESULTS: %d passed, %d failed\n', nPass, nFail);
fprintf('========================================\n');
if nFail == 0
    fprintf('  ALL TESTS PASSED\n');
else
    fprintf('  *** %d TEST(S) FAILED ***\n', nFail);
end
fprintf('========================================\n\n');
