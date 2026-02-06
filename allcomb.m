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
    
    % Convert all inputs to cells for uniform handling
    for i = 1:n
        if ~iscell(args{i})
            args{i} = num2cell(args{i}(:)');
        end
    end
    
    % Get sizes
    sizes = cellfun(@numel, args);
    nRows = prod(sizes);
    
    % Generate combinations using ndgrid
    [F{1:n}] = ndgrid(1:sizes(1), 1:sizes(2:end));
    
    % Build output
    if isstring(varargin{1})
        % String output
        C = strings(nRows, n);
        for i = 1:n
            idx = F{i}(:);
            C(:,i) = string(args{i}(idx));
        end
    elseif isnumeric(varargin{1})
        % Numeric output
        C = zeros(nRows, n);
        for i = 1:n
            idx = F{i}(:);
            C(:,i) = [args{i}{idx}];
        end
    else
        % Cell output
        C = cell(nRows, n);
        for i = 1:n
            idx = F{i}(:);
            C(:,i) = args{i}(idx)';
        end
    end
end
