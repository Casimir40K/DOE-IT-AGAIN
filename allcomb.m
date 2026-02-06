function C = allcomb(varargin)
%allcomb Cartesian product of input vectors/arrays (returns array of strings if inputs are strings)
% Example:
%   C = allcomb(["A","B"], ["X","Y","Z"]);
%   -> 6x2

args = varargin;
n = numel(args);

[F{1:n}] = ndgrid(args{:});
for i = 1:n
    F{i} = F{i}(:);
end
C = [F{:}];
end
