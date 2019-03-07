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
function [K, Kadj, cleanup, ndet] = createRadon3d(m, theta, gpuIdx)
%RADON3D Creates function handles for 3D Radon transform using
%ASTRA Toolbox.
%
%   [K, Kadj, cleanup, ndet] = RADON3D(m, theta, gpuIdx) takes size of 3D
%   geometry m, and angles of projections theta, and returns function 
%   handles and the number of detectors. gpuIdx is the index of the GPU.
%
%   The setup uses parallel beam projections and the detector geometry is
%   set up as [0, 1] and the number of detectors is
%   ndet = ceil(hypot(m(1), m(2))).
%
%   Note that both K and Kadj accept and return vectors! The function 
%   requires GPU as ASTRA only provides a GPU implementation of 3D 
%   operations.
%
% Input:
%   m       size of 3D geometry, vector of length three.
%   theta   angles of projections in degrees (0-180), vector of length k.
%   
% Output:
%   K       function handle, K maps a vector of length prod(m) to a matrix
%           of size [ndet(1), length(theta), ndet(2)].
%   Kadj    function handle, Kadj maps a matrix of size
%           [ndet(1), length(theta), ndet(2)] to a vector of length prod(m).
%   cleanup function handle that can be called to free resources occupied
%           by ASTRA.
%   ndet    The number of detectors. Vector of length two.
%
% Note that ASTRA seems to require m(1) == m(2)!

% Set up detector size and geometries.
ndet = [m(3), ceil(sqrt(sum(m(1:2).^2)))];
vol_geom = astra_create_vol_geom(m);
proj_geom = astra_create_proj_geom('parallel3d', 1.0, 1.0, ndet(1), ndet(2), theta * pi / 180);

% Create operator function handles.
K = @(x) radon3d(reshape(x, m), proj_geom, vol_geom, gpuIdx);
Kadj = @(x) reshape(radon3dadj(x, proj_geom, vol_geom, gpuIdx), [], 1);

% Create cleanup function handle.
cleanup = @() astra_cleanup(proj_geom, vol_geom);

end

function sino = radon3d(data, proj_geom, vol_geom, gpuIdx)
    astra_mex('set_gpu_index', gpuIdx);
    [sino_id, sino] = astra_create_sino3d_cuda(data, proj_geom, vol_geom);
    astra_mex_data3d('delete', sino_id);
end

function vol = radon3dadj(sino, proj_geom, vol_geom, gpuIdx)
    astra_mex('set_gpu_index', gpuIdx);
    [vol_id, vol] = astra_create_backprojection3d_cuda(sino, proj_geom, vol_geom);
    astra_mex_data3d('delete', vol_id);
end

function astra_cleanup(proj_geom, vol_geom)
    astra_mex_data3d('delete', proj_geom, vol_geom);
end
