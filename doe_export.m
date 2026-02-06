function doe_export(design, report, cfg)
%DOE_EXPORT Save run sheet and design bundle to disk
%
%   doe_export(design, report, cfg)
%
%   Saves the DOE design to CSV and/or MAT files based on cfg settings.

    outDir = char(cfg.export.folder);
    
    % Create output directory if it doesn't exist
    if ~exist(outDir, "dir")
        mkdir(outDir);
        fprintf('Created output directory: %s\n', outDir);
    end
    
    base = char(cfg.export.baseFilename);
    T = design.runSheet;
    
    % Export CSV
    if cfg.export.writeCSV
        csvPath = fullfile(outDir, base + "_runsheet.csv");
        writetable(T, csvPath);
        fprintf('✓ Exported CSV: %s\n', csvPath);
    end
    
    % Export MAT bundle
    if cfg.export.writeMAT
        matPath = fullfile(outDir, base + "_bundle.mat");
        save(matPath, "design", "report", "cfg", '-v7.3');
        fprintf('✓ Exported MAT: %s\n', matPath);
    end
    
    fprintf('\n');
end
