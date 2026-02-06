%% CCDESIGN SYNTAX TESTER
% This script tests which ccdesign syntax works in your MATLAB version
% Run this to diagnose the issue, then we'll fix createDOE accordingly

fprintf('Testing ccdesign syntax for your MATLAB version...\n\n');

k = 3;  % Test with 3 factors
alpha = sqrt(3);
cp = 3;

%% Test 1: Name-Value pairs with Capital letters (R2020a+)
fprintf('Test 1: Name-Value pairs (Capital)... ');
try
    X = ccdesign(k, 'Alpha', alpha, 'Center', [cp cp], 'Type', 'circumscribed');
    fprintf('SUCCESS! ✓\n');
    fprintf('  Your MATLAB uses: ccdesign(k, ''Alpha'', val, ''Center'', [n n], ''Type'', type)\n\n');
    disp('First few rows:');
    disp(X(1:5,:));
    return;
catch ME
    fprintf('Failed: %s\n', ME.message);
end

%% Test 2: Name-value pairs with lowercase (older style)
fprintf('Test 2: Name-Value pairs (lowercase)... ');
try
    X = ccdesign(k, 'alpha', alpha, 'center', [cp cp], 'type', 'circumscribed');
    fprintf('SUCCESS! ✓\n');
    fprintf('  Your MATLAB uses: ccdesign(k, ''alpha'', val, ''center'', [n n], ''type'', type)\n\n');
    disp('First few rows:');
    disp(X(1:5,:));
    return;
catch ME
    fprintf('Failed: %s\n', ME.message);
end

%% Test 3: Positional with name-value (mixed)
fprintf('Test 3: Positional + name-value (mixed)... ');
try
    X = ccdesign(k, 'circumscribed', 'alpha', alpha, 'center', [cp cp]);
    fprintf('SUCCESS! ✓\n');
    fprintf('  Your MATLAB uses: ccdesign(k, type, ''alpha'', val, ''center'', [n n])\n\n');
    disp('First few rows:');
    disp(X(1:5,:));
    return;
catch ME
    fprintf('Failed: %s\n', ME.message);
end

%% Test 4: Pure positional arguments
fprintf('Test 4: Pure positional arguments... ');
try
    X = ccdesign(k, alpha, 'circumscribed', [cp cp]);
    fprintf('SUCCESS! ✓\n');
    fprintf('  Your MATLAB uses: ccdesign(k, alpha, type, [n n])\n\n');
    disp('First few rows:');
    disp(X(1:5,:));
    return;
catch ME
    fprintf('Failed: %s\n', ME.message);
end

%% Test 5: Minimal syntax (just k)
fprintf('Test 5: Minimal syntax ccdesign(k)... ');
try
    X = ccdesign(k);
    fprintf('SUCCESS! ✓\n');
    fprintf('  Your MATLAB uses minimal: ccdesign(k)\n');
    fprintf('  Default design created (may need manual adjustment for alpha/centers)\n\n');
    disp('First few rows:');
    disp(X(1:5,:));
    fprintf('\n  Note: This may not support custom alpha or center points directly.\n');
    return;
catch ME
    fprintf('Failed: %s\n', ME.message);
end

%% If we get here, nothing worked
fprintf('\n');
fprintf('========================================\n');
fprintf('ERROR: None of the ccdesign syntaxes worked!\n');
fprintf('========================================\n');
fprintf('This could mean:\n');
fprintf('1. Statistics Toolbox is not installed\n');
fprintf('2. Your MATLAB version is very old/new with different syntax\n');
fprintf('3. ccdesign function is not available\n\n');

fprintf('Checking if ccdesign exists...\n');
if exist('ccdesign', 'file') == 2
    fprintf('  ✓ ccdesign function found\n');
    fprintf('  Try: help ccdesign\n');
    fprintf('  to see the correct syntax for your version\n\n');
else
    fprintf('  ✗ ccdesign function NOT found\n');
    fprintf('  You need to install Statistics and Machine Learning Toolbox\n\n');
end

fprintf('ALTERNATIVE: Use type=''fullfact'' or type=''boxbehnken'' instead of CCD\n');
