# Design of Experiments (DOE) Builder for MATLAB

A comprehensive MATLAB toolkit for creating, exporting, and visualizing Design of Experiments (DOE) for process optimization and research.

## 📋 Overview

This toolkit provides a streamlined way to create various experimental designs including:
- Full Factorial
- Fractional Factorial
- Central Composite Design (CCD)
- Box-Behnken Design
- Latin Hypercube Sampling (LHS)

**Key Features:**
- Mix continuous and categorical factors
- Apply custom constraints
- Automatic randomization and blocking
- Export to CSV and MAT formats
- Built-in visualization tools
- Fully reproducible designs

## 📁 Files Included

1. **createDOE.m** - Main function to create DOE designs
2. **allcomb.m** - Helper function for Cartesian products
3. **doe_export.m** - Export designs to CSV/MAT files
4. **doe_visualize.m** - Visualization tools
5. **RUNSCRIPT.m** - Example script for lithium formate synthesis DOE
6. **README.md** - This file

## 🚀 Quick Start

### Basic Usage

```matlab
% 1. Define your factors
factors = [
    struct('name', 'Temperature', ...
           'type', 'continuous', ...
           'low', 60, ...
           'high', 90, ...
           'units', 'degC', ...
           'levels', [], ...
           'transform', 'none')
    
    struct('name', 'Catalyst', ...
           'type', 'categorical', ...
           'levels', ["A", "B", "C"], ...
           'low', [], ...
           'high', [], ...
           'units', '-', ...
           'transform', 'none')
];

% 2. Create the DOE
[design, report, cfg] = createDOE(factors, ...
    'type', 'ccd', ...
    'title', 'Temperature Optimization', ...  % Optional: custom display title
    'exportCSV', true);

% 3. View the results
disp(report.summaryText);
disp(design.runSheet);

% 4. Visualize
doe_visualize(design, factors);
```

## 📖 Detailed Usage

### Factor Definition

Each factor requires a struct with these fields:

**Continuous Factors:**
```matlab
struct('name', 'FactorName', ...      % String identifier
       'type', 'continuous', ...       % Factor type
       'low', 10, ...                  % Lower bound
       'high', 50, ...                 % Upper bound
       'units', 'degC', ...            % Units (for display)
       'levels', [], ...               % Leave empty
       'transform', 'none')            % 'none', 'log', 'sqrt'
```

**Categorical Factors:**
```matlab
struct('name', 'FactorName', ...      % String identifier
       'type', 'categorical', ...      % Factor type
       'levels', ["A", "B", "C"], ... % Level names
       'low', [], ...                  % Leave empty
       'high', [], ...                 % Leave empty
       'units', '-', ...               % Units (typically '-')
       'transform', 'none')            % Leave as 'none'
```

### Design Types

#### 1. Full Factorial (`'fullfact'`)
Tests all combinations of factor levels.

```matlab
[design, report, cfg] = createDOE(factors, ...
    'type', 'fullfact', ...
    'centerPoints', 3);
```

**Best for:** 
- 2-4 factors
- Small experiments
- Complete coverage

#### 2. Fractional Factorial (`'fracfact'`)
Tests a fraction of all combinations using confounding.

```matlab
[design, report, cfg] = createDOE(factors, ...
    'type', 'fracfact', ...
    'fracGenerators', 'a b c abc');
```

**Best for:**
- Screening many factors
- Limited resources
- Main effects + some interactions

#### 3. Central Composite Design (`'ccd'`)
Factorial + axial points + center points for response surface modeling.

```matlab
[design, report, cfg] = createDOE(factors, ...
    'type', 'ccd', ...
    'ccdAlpha', 'rotatable', ...  % or 'orthogonal', 'face', or numeric
    'ccdFaceType', 'circumscribed', ...
    'centerPoints', 6);
```

**Best for:** 
- Response surface methodology
- Quadratic models
- Optimization studies

**Alpha options:**
- `'rotatable'` - Equal prediction variance at equal distances from center (alpha = sqrt(k))
- `'orthogonal'` - Minimizes correlation between coefficients
- `'face'` - Axial points on cube faces (alpha = 1)
- Numeric value - Custom alpha distance

**Note:** The code automatically handles different MATLAB versions' ccdesign syntax.

#### 4. Box-Behnken Design (`'boxbehnken'`)
Efficient design for quadratic models (requires ≥3 factors).

```matlab
[design, report, cfg] = createDOE(factors, ...
    'type', 'boxbehnken', ...
    'centerPoints', 5);
```

**Best for:**
- 3+ factors
- Avoiding extreme combinations
- Quadratic models with fewer runs than CCD

#### 5. Latin Hypercube Sampling (`'lhs'`)
Space-filling design for computer experiments.

```matlab
[design, report, cfg] = createDOE(factors, ...
    'type', 'lhs', ...
    'lhsSamples', 50, ...
    'lhsCriterion', 'maximin');
```

**Best for:**
- Computer simulations
- Many factors
- Uniform space coverage

### Advanced Features

#### Constraints
Exclude infeasible combinations:

```matlab
constraints = {};

% Temperature constraint
constraints{1} = @(T) T.Temperature <= 85 | T.Pressure <= 10;

% Categorical constraint
constraints{2} = @(T) ~(strcmp(T.Catalyst, "Expensive") & T.Time > 120);

[design, report, cfg] = createDOE(factors, ...
    'constraints', constraints);
```

#### Blocking
Organize runs into blocks (e.g., different days/batches):

```matlab
[design, report, cfg] = createDOE(factors, ...
    'numBlocks', 3, ...
    'blockingMethod', 'sequential');  % or 'random'
```

#### Replicates
Add complete replicates for better error estimation:

```matlab
[design, report, cfg] = createDOE(factors, ...
    'replicates', 3);
```

#### Randomization
Control run order randomization:

```matlab
[design, report, cfg] = createDOE(factors, ...
    'randomize', true, ...
    'randomSeed', 42);  % For reproducibility
```

## 📊 Outputs

### design struct
```matlab
design.runSheet  % Table with all runs
design.meta      % Project metadata
design.cfg       % Configuration (for reproducibility)
```

### runSheet columns
- `Run` - Run number (randomized order)
- `FactorName_coded` - Coded values (-1 to +1) for continuous factors
- `FactorName` - Actual values for continuous factors
- `CategoricalName` - Levels for categorical factors
- `Block` - Block assignment
- `Replicate` - Replicate number
- `Response1`, `Response2` - Placeholders for experimental results

### report struct
```matlab
report.project      % Project name
report.designType   % Design type used
report.nRuns        % Total number of runs
report.kCont        % Number of continuous factors
report.kCat         % Number of categorical factors
report.blocks       % Number of blocks
report.replicates   % Number of replicates
report.summaryText  % Formatted summary
```

## 📁 Export Options

```matlab
[design, report, cfg] = createDOE(factors, ...
    'exportFolder', './MyResults', ...
    'exportBaseFilename', 'MyExperiment', ...
    'exportCSV', true, ...  % Creates MyExperiment_runsheet.csv
    'exportMAT', true);     % Creates MyExperiment_bundle.mat
```

**Output files:**
- `*_runsheet.csv` - Run sheet for lab use
- `*_bundle.mat` - Complete design + config for MATLAB

## 📈 Visualization

```matlab
doe_visualize(design, factors);
```

**Creates:**
- Histograms of continuous factor distributions
- Pairwise scatter plots
- 3D design space visualization (if 3+ factors)
- Categorical factor level counts
- Summary statistics

## 💡 Example: Lithium Formate Synthesis

See `RUNSCRIPT.m` for a complete example with:
- 1 categorical factor (Mode: Co-only, Mn-only, Together)
- 4 continuous factors (Temperature, Time, Concentration, S/L Ratio)
- CCD with rotatable alpha
- Constraint examples
- Full visualization

**To run:**
```matlab
RUNSCRIPT
```

## ⚙️ Requirements

- MATLAB R2019b or later
- **For CCD and Box-Behnken:** Statistics and Machine Learning Toolbox
- **For other designs:** Base MATLAB only

## 🔍 Tips & Best Practices

1. **Start small:** Begin with 1 replicate and add more if needed
2. **Center points:** 3-6 center points help detect curvature
3. **Blocking:** Use blocks if experiments span multiple days/batches
4. **Constraints:** Add constraints for safety or feasibility
5. **Randomization:** Always randomize unless there's a specific reason not to
6. **Coded values:** Use for balanced model fitting, actual values for interpretation

## 🐛 Troubleshooting

**Error: "Statistics and Machine Learning Toolbox required"**
- Solution: Use 'fullfact', 'fracfact', or 'lhs' designs, or install the toolbox

**Error: "Factor names must be unique"**
- Solution: Check that all factor names are different

**Error: "high > low required"**
- Solution: Verify continuous factor ranges are valid

**Too many runs generated:**
- Solution: Use fractional factorial, add constraints, or reduce factor levels

**Design doesn't match expectations:**
- Solution: Check `design.runSheet` and verify factor definitions

## 📚 References

- Box, G. E., Hunter, J. S., & Hunter, W. G. (2005). Statistics for Experimenters
- Montgomery, D. C. (2017). Design and Analysis of Experiments
- Myers, R. H., Montgomery, D. C., & Anderson-Cook, C. M. (2016). Response Surface Methodology

## 📄 License

MIT License - Feel free to use and modify for your research.

## 🤝 Contributing

Suggestions and improvements welcome! Key areas:
- Additional design types
- More visualization options
- Export to other formats
- Integration with analysis tools

## 📧 Support

For issues or questions, please check:
1. This README
2. Function documentation (`help createDOE`)
3. Example script (RUNSCRIPT.m)

---

**Version:** 2.0  
**Last Updated:** February 2026
