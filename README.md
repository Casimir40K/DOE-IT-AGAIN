# DOE-IT-AGAIN
DOE solver for matlab, because I didn't want to go to the PhD room to use MiniTab.

% =========================================================================
%  createDOE
% =========================================================================
%  General-purpose Design of Experiments (DOE) generator.
%
%  This function provides a reusable, configurable framework for generating
%  experimental designs across a wide range of applications (chemical,
%  mechanical, process, and data-driven experimentation).
%
%  The solver supports common DOE methodologies and produces a standardized,
%  reproducible run sheet including coded and actual factor values, optional
%  blocking, randomization, replication, and constraint handling.
%
% =========================================================================
%  FUNCTION SIGNATURE
% =========================================================================
%  [design, report, cfg] = createDOE(factors, Name,Value,...)
%
% =========================================================================
%  REQUIRED INPUT
% =========================================================================
%  factors
%    Struct array defining all experimental factors.
%
%    Each factor MUST include:
%      - name   : string
%                Unique factor name. Should be a valid MATLAB identifier
%                for clean table variable naming.
%
%      - type   : "continuous" or "categorical"
%
%    Continuous factors REQUIRE:
%      - low    : numeric scalar
%                Lower bound of the factor (actual units).
%
%      - high   : numeric scalar
%                Upper bound of the factor (actual units).
%
%      - units  : string
%                Documentation only (not used in calculations).
%
%    Categorical factors REQUIRE:
%      - levels : string array
%                Allowed discrete values.
%
%    Example:
%      factors = [
%        struct("name","Temperature","type","continuous", ...
%               "low",40,"high",90,"units","degC", ...
%               "levels",string.empty,"transform","none")
%
%        struct("name","Solvent","type","categorical", ...
%               "low",[],"high",[],"units","-", ...
%               "levels",["Water","MeOH","EtOH"], ...
%               "transform","none")
%      ];
%
% =========================================================================
%  DESIGN OPTIONS (NAME–VALUE PAIRS)
% =========================================================================
%  General design settings:
%
%    "type"
%       DOE type. One of:
%         "fullfact"    - 2-level full factorial
%         "fracfact"    - fractional factorial
%         "ccd"         - central composite design
%         "boxbehnken"  - Box–Behnken design
%         "lhs"         - Latin hypercube sampling
%       Default: "ccd"
%
%    "centerPoints"
%       Number of center points (where applicable).
%       Default: 0
%
%    "replicates"
%       Number of whole-design replications.
%       Default: 1
%
%    "randomize"
%       Logical flag to randomize run order.
%       Default: true
%
%    "randomSeed"
%       Random number generator seed for reproducibility.
%       Default: 42
%
%    "numBlocks"
%       Number of experimental blocks.
%       Default: 1 (no blocking)
%
%    "blockingMethod"
%       Method for assigning runs to blocks:
%         "sequential" | "random"
%       Default: "sequential"
%
% -------------------------------------------------------------------------
%  Method-specific options
% -------------------------------------------------------------------------
%
%  Fractional factorial:
%    "fracGenerators"
%       Generator string passed to fracfact().
%
%  Central composite design (CCD):
%    "ccdAlpha"
%       Axial distance:
%         "rotatable" | "orthogonal" | numeric
%
%    "ccdFaceType"
%       CCD geometry:
%         "circumscribed" | "inscribed" | "faced"
%
%  Latin hypercube sampling:
%    "lhsSamples"
%       Number of LHS points.
%
%    "lhsCriterion"
%       LHS optimization criterion (e.g. "maximin").
%
% =========================================================================
%  CONSTRAINT HANDLING
% =========================================================================
%  "constraints"
%    Cell array of function handles used to remove invalid runs.
%
%    Each constraint function:
%      - Accepts the RUNSHEET TABLE (actual units)
%      - Returns:
%           * logical scalar, or
%           * logical vector of length height(runSheet)
%      - TRUE  → keep run
%      - FALSE → discard run
%
%    Constraints are applied after DOE generation and before replication
%    and randomization.
%
%    Example:
%      constraints = {
%        @(T) T.Temperature <= 85
%        @(T) ~(T.Solvent=="MeOH" & T.pH > 8)
%      };
%
% =========================================================================
%  OUTPUTS
% =========================================================================
%  design
%    Struct containing:
%      - runSheet : table of experimental runs including:
%          * Run index
%          * Coded continuous factors
%          * Actual continuous factors
%          * Categorical factors
%          * Block number
%          * Replicate number
%          * Placeholder response columns
%
%      - meta     : project metadata
%      - cfg      : full configuration used to generate the DOE
%
%  report
%    Struct summarizing the design:
%      - Design type
%      - Number of runs
%      - Number of continuous and categorical factors
%      - Blocking and replication settings
%      - Build time
%
%  cfg
%    Complete configuration struct used internally. Returned to ensure
%    full reproducibility of the design.
%
% =========================================================================
%  EXPORT OPTIONS
% =========================================================================
%  "exportCSV"
%     Write run sheet to CSV file.
%
%  "exportMAT"
%     Save full DOE bundle (design, report, cfg) to MAT-file.
%
%  "exportFolder"
%     Output directory for exported files.
%
%  "exportBaseFilename"
%     Base filename used for exports.
%
% =========================================================================
%  TOOLBOX REQUIREMENTS
% =========================================================================
%  - Full factorial, fractional factorial, LHS:
%       Base MATLAB
%
%  - Central composite and Box–Behnken designs:
%       Statistics and Machine Learning Toolbox
%       (functions: ccdesign, bbdesign)
%
% =========================================================================
%  DESIGN PHILOSOPHY
% =========================================================================
%  - This function is intended to be reused without modification.
%  - Experimental definition should be isolated to factor definitions and
%    Name–Value options in the calling script.
%  - All outputs are reproducible when randomSeed is fixed.
%
% =========================================================================
%  AUTHOR / VERSION
% =========================================================================
%  Author  : Casimir Solomon
%  Version : 1.0
%  Date    : 06/02/2026
%
% =========================================================================
