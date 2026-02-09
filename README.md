# DOE Builder for MATLAB

A single-file MATLAB function for generating experimental designs. No Statistics Toolbox required — every design generator is built from scratch.

Written with substantial help from Claude (Anthropic). If something is elegant, that's probably Claude. If something is oddly specific to lithium formate leaching, that's me.

## What it does

`createDOE` takes a list of factors (continuous, categorical, or both), a design type, and some options, and returns a randomised run-sheet as a MATLAB table. It handles the crossing of categorical and continuous factors, coded-to-actual value mapping, blocking, replication, constraints, and CSV/MAT export.

Supported designs:

| Design | Alias(es) | Minimum factors | Use case |
|--------|-----------|-----------------|----------|
| Full factorial | `fullfact`, `ff`, `factorial` | 1 continuous or 1 categorical | Complete coverage, few factors |
| Fractional factorial | `fracfact`, `frac` | 1 continuous | Screening, many factors |
| Central Composite (CCD) | `ccd`, `ccdesign` | 1 continuous | Response surface / quadratic models |
| Box–Behnken | `bbd`, `boxbehnken`, `bb` | 3 continuous | Response surface, fewer runs than CCD |
| Plackett–Burman | `pbdesign`, `pb`, `plackett` | 1 continuous | Main-effects-only screening |
| Latin Hypercube | `lhs`, `latinhypercube` | 1 continuous | Computer experiments, space-filling |

## Files

| File | What it is |
|------|------------|
| `createDOE.m` | The function. This is the only file you need. |
| `RUNSCRIPT.m` | Example: lithium formate synthesis DOE (CCD + categorical). |
| `QuickReference.m` | Copy-pasteable recipes for common design patterns. |
| `test_createDOE.m` | 22 automated tests. Run this first if you're sceptical. |
| `README.md` | You're here. |

## Quick start

```matlab
% Two continuous factors, full factorial
factors = [ struct('name','Temp','type','continuous','low',60,'high',90)
            struct('name','Time','type','continuous','low',30,'high',180) ];

[design, report] = createDOE(factors, 'type','fullfact');
disp(design.runSheet);
```

## Defining factors

Continuous factors need `name`, `type`, `low`, and `high`. Categorical factors need `name`, `type`, and `levels`. The `units` field is optional (defaults to `'-'`).

```matlab
% All-continuous — struct array is fine
factors = [ struct('name','Temp', 'type','continuous', 'low',60, 'high',90)
            struct('name','Conc', 'type','continuous', 'low',0.1,'high',1.0) ];
```

When mixing categorical and continuous factors, use a **cell array** instead of a struct array. MATLAB won't let you concatenate structs with different fields, and categorical factors don't have `low`/`high` while continuous ones don't have `levels`. The function normalises everything internally.

```matlab
% Mixed — use a cell array
factors = {
    struct('name','Catalyst', 'type','categorical', 'levels',["A","B","C"])
    struct('name','Temp',     'type','continuous',   'low',60, 'high',90)
    struct('name','Time',     'type','continuous',   'low',30, 'high',180)
};
```

## Options

Pass these as name-value pairs after the factors argument.

**Design control:**

| Parameter | Default | Notes |
|-----------|---------|-------|
| `'type'` | `'fullfact'` | See the table above for aliases |
| `'centerPoints'` | `0` | Added to the base design before crossing with categoricals |
| `'replicates'` | `1` | Complete replicates of the whole design |
| `'numBlocks'` | `1` | Sequential split into blocks |
| `'randomize'` | `true` | Shuffles within each block |
| `'randomSeed'` | `42` | For reproducibility |

**Design-specific:**

| Parameter | Applies to | Default | Notes |
|-----------|-----------|---------|-------|
| `'ccdAlpha'` | CCD | `'rotatable'` | Also accepts `'face'`, `'orthogonal'`, or a number |
| `'fracGenerators'` | Fractional factorial | `''` | e.g. `'a b c abc'` for a 2^(4-1) design |
| `'lhsSamples'` | LHS | `50` | Number of sample points |
| `'lhsMaxIter'` | LHS | `20` | Maximin optimisation iterations |

**Constraints:**

```matlab
% Cell array of functions. Each receives the run-sheet table,
% returns a logical vector: true = keep, false = discard.
constraints = { @(T) ~(T.Temp > 85 & T.Conc > 0.9) };

[d,r] = createDOE(factors, 'type','ccd', 'constraints',constraints);
```

**Export and metadata:**

| Parameter | Default | Notes |
|-----------|---------|-------|
| `'projectName'` | `'DOE_Project'` | Used for filenames |
| `'title'` | same as projectName | Display title in console output |
| `'owner'` | `''` | Stored in metadata |
| `'notes'` | `''` | Stored in metadata |
| `'exportFolder'` | `'./DOE_Output'` | Created if it doesn't exist |
| `'exportCSV'` | `false` | Writes `<projectName>_runsheet.csv` |
| `'exportMAT'` | `false` | Writes `<projectName>_bundle.mat` |
| `'responses'` | `{'Response1','Response2'}` | Column names for response placeholders |

## What you get back

```matlab
[design, report, cfg] = createDOE(factors, ...);
```

`design.runSheet` is a MATLAB table with columns: `Run`, `Block`, `Replicate`, then actual-value columns for each factor, coded-value columns (`*_coded`) for continuous factors, categorical columns, and NaN-filled response placeholders.

`report` is a struct with `nRuns`, `designType`, `kCont`, `kCat`, `summaryText`, etc.

`cfg` is the full configuration — enough to reproduce the design exactly.

## CCD alpha values

For the record, since this tripped me up and it'll trip you up too:

| Alpha setting | Formula | k=2 | k=3 | k=4 |
|---------------|---------|-----|-----|-----|
| `'rotatable'` | (2^k)^(1/4) | 1.414 | 1.682 | 2.000 |
| `'face'` | 1 | 1 | 1 | 1 |
| `'orthogonal'` | (k(1+√(1+8k))/4)^(1/4) | varies | varies | varies |

A lot of references (and the previous version of this code) use `sqrt(k)` for rotatable alpha, which only happens to be correct for k=4. The proper formula is `(2^k)^(1/4)`.

## Requirements

MATLAB R2019b or later. No toolboxes.

The old version of this toolkit depended on `ccdesign`, `bbdesign`, and `lhsdesign` from the Statistics Toolbox. This version implements everything from scratch, so those functions are not needed. If you have the toolbox installed, it won't interfere — nothing from it is called.

## Running the tests

```matlab
test_createDOE
```

This runs 22 tests covering every design type, edge cases (single-factor CCD, categorical-only designs), constraint filtering, replication, blocking, randomisation reproducibility, and coded-to-actual mapping. It prints a pass/fail count at the end.

## References

- Montgomery, D.C. (2017). *Design and Analysis of Experiments*, 9th ed. Wiley.
- Box, G.E.P., Hunter, J.S., & Hunter, W.G. (2005). *Statistics for Experimenters*, 2nd ed. Wiley.
- Myers, R.H., Montgomery, D.C., & Anderson-Cook, C.M. (2016). *Response Surface Methodology*, 4th ed. Wiley.

## Licence

MIT. Do whatever you like with it.
