% Copyright 2019 Lukas F. Lang and Sebastian Neumayer
%
% This file is part of TBIR.
%
%    TBIR is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    TBIR is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with TBIR.  If not, see <http://www.gnu.org/licenses/>.
function y = flatten3d(x, cols)
%FLATTEN3D Takes a 3D array and returns a flattened 2D version.

m = size(x);
maxval = max(x(:));
c = squeeze(mat2cell(x, m(1), m(2), ones(m(3), 1)));
rows = ceil(m(3) / cols);

% Pad missing images.
for k=m(3)+1:rows*cols
    c{k} = zeros(m(1:2));
end
c = reshape(c, cols, rows)';

% Add white frame.
c = cellfun(@(x) padarray(x, [1, 1], maxval), c, 'UniformOutput', false);
y = cell2mat(c);

end