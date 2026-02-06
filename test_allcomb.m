%% Test allcomb function
fprintf('Testing allcomb function...\n\n');

%% Test 1: String arrays
fprintf('Test 1: String arrays\n');
modes = ["Co-only", "Mn-only", "Together"];
fprintf('Input: ');
disp(modes);

try
    result = allcomb(modes);
    fprintf('Output (%dx%d):\n', size(result, 1), size(result, 2));
    disp(result);
    fprintf('✓ SUCCESS\n\n');
catch ME
    fprintf('✗ FAILED: %s\n\n', ME.message);
end

%% Test 2: Two string arrays
fprintf('Test 2: Two string arrays (cross product)\n');
cat1 = ["A", "B"];
cat2 = ["X", "Y", "Z"];
fprintf('Input 1: ');
disp(cat1);
fprintf('Input 2: ');
disp(cat2);

try
    result = allcomb(cat1, cat2);
    fprintf('Output (%dx%d):\n', size(result, 1), size(result, 2));
    disp(result);
    fprintf('✓ SUCCESS\n\n');
catch ME
    fprintf('✗ FAILED: %s\n\n', ME.message);
end

%% Test 3: Numeric arrays
fprintf('Test 3: Numeric arrays\n');
num1 = [1, 2];
num2 = [10, 20, 30];
fprintf('Input 1: ');
disp(num1);
fprintf('Input 2: ');
disp(num2);

try
    result = allcomb(num1, num2);
    fprintf('Output (%dx%d):\n', size(result, 1), size(result, 2));
    disp(result);
    fprintf('✓ SUCCESS\n\n');
catch ME
    fprintf('✗ FAILED: %s\n\n', ME.message);
end

%% Test 4: Single element
fprintf('Test 4: Single string array (edge case)\n');
single = ["OnlyOne"];
fprintf('Input: ');
disp(single);

try
    result = allcomb(single);
    fprintf('Output (%dx%d):\n', size(result, 1), size(result, 2));
    disp(result);
    fprintf('✓ SUCCESS\n\n');
catch ME
    fprintf('✗ FAILED: %s\n\n', ME.message);
end

fprintf('All tests complete!\n');
