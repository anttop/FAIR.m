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
function [K, Kadj, cleanup, ndet] = createRadon2d(m, theta, scale)
%RADON2D Creates function handles for 2D Radon transform using
%ASTRA Toolbox.
%
%   [K, Kadj, cleanup, ndet] = RADON2D(m, theta, scale) takes size of 2D
%   geometry m, and angles of projections theta, and a scaling factor,
%   and returns function handles and the number of detectors.
%
%   The setup uses parallel beam projections and the detector geometry is
%   set up as [0, 1] and the number of detectors is
%   ndet = ceil(hypot(m(1), m(2))).
%
%   Note that both K and Kadj accept and return vectors!
%
% Input:
%   m       size of 2D geometry, vector of length two.
%   theta   angles of projections in degrees (0-180), vector of length k.
%   scale   a scalar > 0.
%   
% Output:
%   K       function handle, K maps a vector of length prod(m) to a matrix
%           of size [length(theta), ndet].
%   Kadj    function handle, Kadj maps a matrix of size
%           [length(theta), ndet] to a vector of length prod(m).
%   cleanup function handle that can be called to free resources occupied
%           by ASTRA.
%   ndet    The number of detectors.

% Set up detector size and geometries.
ndet = 1.5 * m(1);
vol_geom = astra_create_vol_geom(m);
proj_geom = astra_create_proj_geom('parallel', 1.0, ndet, deg2rad(theta));

% Create operator function handles.
K = @(x) scale * radon2d(reshape(x, m), proj_geom, vol_geom);
Kadj = @(x) reshape(radon2dadj(scale * x, proj_geom, vol_geom), [], 1);

% Create cleanup function handle.
cleanup = @() astra_cleanup(proj_geom, vol_geom);

end

function sino = radon2d(data, proj_geom, vol_geom)
    proj_id = astra_create_projector('linear', proj_geom, vol_geom);
    [sino_id, sino] = astra_create_sino(data, proj_id);
    astra_mex_data2d('delete', sino_id);
    astra_mex_projector('delete', proj_id);
end

function vol = radon2dadj(sino, proj_geom, vol_geom)
    recon_id = astra_mex_data2d('create', '-vol', vol_geom, 0);
    sinogram_id = astra_mex_data2d('create', '-sino', proj_geom, sino);
    proj_id = astra_create_projector('linear', proj_geom, vol_geom);
    
    cfg = astra_struct('BP');
    cfg.ProjectorId = proj_id;
    cfg.ProjectionDataId = sinogram_id;
    cfg.ReconstructionDataId = recon_id;

    alg_id = astra_mex_algorithm('create', cfg);
    astra_mex_algorithm('run', alg_id);
    vol = astra_mex_data2d('get', recon_id);

    astra_mex_data2d('delete', sinogram_id, recon_id);
    astra_mex_projector('delete', proj_id);
    astra_mex_algorithm('delete', alg_id);
end

function astra_cleanup(proj_geom, vol_geom)
    astra_mex_data2d('delete', proj_geom, vol_geom);
end
