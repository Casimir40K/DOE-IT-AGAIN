function C = allcomb(varargin)
%ALLCOMB Cartesian product of input vectors/arrays
%   C = allcomb(V1, V2, ..., Vn) returns all combinations of elements from
%   the input vectors. Works with numeric arrays, strings, and cells.
%
%   Example:
%       C = allcomb(["A","B"], ["X","Y","Z"]);
%       % Returns 6x2 string array with all combinations

    if nargin == 0
        C = [];
        return;
    end
    
    args = varargin;
    n = numel(args);
    
    % Get sizes
    sizes = zeros(1, n);
    for i = 1:n
        if isstring(args{i}) || ischar(args{i})
            args{i} = string(args{i}(:));  % Ensure column vector of strings
        end
        sizes(i) = numel(args{i});
    end
    
    nRows = prod(sizes);
    
    % Handle edge case
    if nRows == 0
        if isstring(args{1})
            C = strings(0, n);
        else
            C = zeros(0, n);
        end
        return;
    end
    
    % Generate combinations using ndgrid
    grids = cell(1, n);
    [grids{:}] = ndgrid(1:sizes(1), 1:sizes(2:end));
    
    % Build output based on first argument type
    if isstring(args{1})
        % String output
        C = strings(nRows, n);
        for i = 1:n
            idx = grids{i}(:);
            C(:, i) = args{i}(idx);
        end
    elseif isnumeric(args{1})
        % Numeric output
        C = zeros(nRows, n);
        for i = 1:n
            idx = grids{i}(:);
            C(:, i) = args{i}(idx);
        end
    else
        % Cell output
        C = cell(nRows, n);
        for i = 1:n
            idx = grids{i}(:);
            C(:, i) = args{i}(idx);
        end
    end
end
