%% ========================================================================
%  LITHIUM FORMATE DOE - RUN SCRIPT
%  ========================================================================
%  This script creates a Design of Experiments for lithium formate synthesis
%  using a Central Composite Design (CCD) with categorical factors.
%
%  FACTORS:
%    - Mode (Categorical): Co-only, Mn-only, Together
%    - Temperature (Continuous): 60-90°C
%    - Time (Continuous): 30-180 min
%    - Concentration (Continuous): 0.20-1.00 M
%    - S/L Ratio (Continuous): 0.02-0.20 g/mL
%
%  DESIGN: CCD with alpha=2 (rotatable design)
%  ========================================================================

clear; clc; close all;

%% 1. DEFINE FACTOR RANGES
% -------------------------------------------------------------------------
% Temperature - explicitly defined from report
T_low  = 60;    % °C
T_high = 90;    % °C

% Time - adjust these based on your experimental knowledge
t_min  = 30;    % min
t_max  = 180;   % min

% Concentration - adjust these based on your experimental knowledge
C_min  = 0.20;  % M
C_max  = 1.00;  % M

% Solid/Liquid Ratio - adjust these based on your experimental knowledge
R_min  = 0.02;  % g/mL
R_max  = 0.20;  % g/mL

% Categorical levels
modes = ["Co-only", "Mn-only", "Together"];

%% 2. CREATE FACTORS STRUCTURE
% -------------------------------------------------------------------------
% Each factor needs: name, type, and type-specific parameters
% Continuous: low, high, units
% Categorical: levels

factors = [
    struct('name', 'Mode', ...
           'type', 'categorical', ...
           'levels', modes, ...
           'low', [], ...
           'high', [], ...
           'units', '-', ...
           'transform', 'none')
    
    struct('name', 'Temperature', ...
           'type', 'continuous', ...
           'low', T_low, ...
           'high', T_high, ...
           'units', 'degC', ...
           'levels', [], ...
           'transform', 'none')
    
    struct('name', 'Time', ...
           'type', 'continuous', ...
           'low', t_min, ...
           'high', t_max, ...
           'units', 'min', ...
           'levels', [], ...
           'transform', 'none')
    
    struct('name', 'Concentration', ...
           'type', 'continuous', ...
           'low', C_min, ...
           'high', C_max, ...
           'units', 'M', ...
           'levels', [], ...
           'transform', 'none')
    
    struct('name', 'SLRatio', ...
           'type', 'continuous', ...
           'low', R_min, ...
           'high', R_max, ...
           'units', 'g/mL', ...
           'levels', [], ...
           'transform', 'none')
];

%% 3. DEFINE CONSTRAINTS (OPTIONAL)
% -------------------------------------------------------------------------
% Add constraint functions to exclude infeasible combinations
% Each function receives the runSheet table and returns logical array

constraints = {};

% Example constraint: Exclude high temperature + high concentration
% Uncomment and modify as needed:
% constraints{1} = @(T) ~(T.Temperature > 85 & T.Concentration > 0.8);

% Example constraint: Specific mode combinations
% constraints{2} = @(T) ~(strcmp(T.Mode, "Together") & T.Time < 60);

%% 4. CREATE THE DOE
% -------------------------------------------------------------------------
fprintf('\n');
fprintf('========================================\n');
fprintf('CREATING LITHIUM FORMATE DOE\n');
fprintf('========================================\n\n');

[design, report, cfg] = createDOE(factors, ...
    'type', 'ccd', ...                    % Central Composite Design
    'ccdAlpha', 2, ...                    % Rotatable design (alpha = 2)
    'ccdFaceType', 'circumscribed', ...   % Face-centered, inscribed, or circumscribed
    'centerPoints', 6, ...                % Number of center points (recommended: 3-6)
    'replicates', 1, ...                  % Number of complete replicates
    'numBlocks', 1, ...                   % Number of blocks
    'randomize', true, ...                % Randomize run order
    'randomSeed', 42, ...                 % For reproducibility
    'constraints', constraints, ...       % Apply any constraints
    'projectName', 'LithiumFormate_CCD', ...
    'owner', 'YourName', ...              % Add your name here
    'notes', 'CCD design for lithium formate synthesis optimization', ...
    'exportFolder', './DOE_Output', ...   % Output folder
    'exportCSV', true, ...                % Export to CSV
    'exportMAT', true);                   % Export to MAT

%% 5. DISPLAY SUMMARY
% -------------------------------------------------------------------------
fprintf('%s\n', report.summaryText);

%% 6. PREVIEW RUN SHEET
% -------------------------------------------------------------------------
fprintf('========================================\n');
fprintf('RUN SHEET PREVIEW (first 10 runs)\n');
fprintf('========================================\n\n');
disp(design.runSheet(1:min(10, height(design.runSheet)), :));

%% 7. VISUALIZE DESIGN
% -------------------------------------------------------------------------
fprintf('\nGenerating visualization plots...\n');
doe_visualize(design, factors);

%% 8. SAVE WORKSPACE
% -------------------------------------------------------------------------
fprintf('\nSaving workspace...\n');
save('./DOE_Output/workspace_LithiumFormate.mat', ...
    'design', 'report', 'cfg', 'factors');
fprintf('✓ Workspace saved\n\n');

%% 9. INSTRUCTIONS FOR NEXT STEPS
% -------------------------------------------------------------------------
fprintf('========================================\n');
fprintf('NEXT STEPS\n');
fprintf('========================================\n');
fprintf('1. Review the run sheet in: ./DOE_Output/\n');
fprintf('2. Run experiments according to randomized order\n');
fprintf('3. Record responses in Response1, Response2 columns\n');
fprintf('4. Use statistical analysis tools to build models\n');
fprintf('========================================\n\n');

%% NOTES
% -------------------------------------------------------------------------
% CCD with 4 continuous factors creates:
%   - 2^4 = 16 factorial points (corners of hypercube)
%   - 2*4 = 8 axial points (star points at ±alpha)
%   - n center points (your choice, typically 3-6)
%   Total: 24 + center points
%
% With 3 categorical levels (Mode), total runs = (24 + center) * 3
%
% Alpha = 2 creates a "rotatable" design where prediction variance
% depends only on distance from center, not direction
%
% TIPS:
%   - Start with 1 replicate to minimize experiments
%   - Add center points to check for curvature
%   - Use constraints to exclude infeasible regions
%   - Consider blocking if experiments span multiple days/batches
