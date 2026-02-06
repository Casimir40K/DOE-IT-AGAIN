function doe_export(design, report, cfg)
%doe_export Save run sheet + bundle.

outDir = char(cfg.export.folder);
if ~exist(outDir, "dir")
    mkdir(outDir);
end

base = char(cfg.export.baseFilename);
T = design.runSheet;

if cfg.export.writeCSV
    csvPath = fullfile(outDir, base + "_runsheet.csv");
    writetable(T, csvPath);
    fprintf("Wrote CSV: %s\n", csvPath);
end

if cfg.export.writeMAT
    matPath = fullfile(outDir, base + "_bundle.mat");
    save(matPath, "design", "report", "cfg");
    fprintf("Wrote MAT: %s\n", matPath);
end
end
