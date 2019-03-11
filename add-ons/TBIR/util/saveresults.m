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
function [resfile, paramfile] = saveresults(name, outputfolder, image1,...
    image2, ~, rec, dist, reg, objfun, imageModel, N, nt, alpha, theta,...
    sigma, elapsed)
%SAVERESULTS A script that saves results and parameters.

% Create string for result filename.
str = sprintf('%s_params_%s_%s_%s_%.2f', name, dist, reg, objfun, sigma);
resfile = fullfile(outputfolder, sprintf('%s.png', str));
str = sprintf('%s_result_%s_%s_%s_%.2f', name, dist, reg, objfun, sigma);
paramfile = fullfile(outputfolder, sprintf('%s.txt', str));

% Save result.
imwrite(rec / max(max(rec(:)), 255), resfile);

% Save parameters.
fid = fopen(paramfile, 'wt+');
fprintf(fid, sprintf('Distance:\t %s\n',dist));
fprintf(fid, sprintf('Regulariser:\t %s\n', reg));
fprintf(fid, sprintf('Objective:\t %s\n', objfun));
fprintf(fid, sprintf('Image model:\t %s\n', imageModel));
fprintf(fid, sprintf('RK steps:\t %i\n', N));
fprintf(fid, sprintf('Time steps:\t %i\n', nt));
for k=1:length(alpha)
    fprintf(fid, sprintf('Alpha %i:\t %g\n', k, alpha(k)));
end
fprintf(fid, sprintf('Angles (deg.):\t %s\n', num2str(theta)));
fprintf(fid, sprintf('Noise level:\t %g\n', sigma));
fprintf(fid, sprintf('SSIM (recon.):\t %g\n', ssim(rec, image2)));
fprintf(fid, sprintf('SSIM (ref.):\t %g\n', ssim(image1, image2)));
fprintf(fid, sprintf('PSNR (recon.):\t %g\n', psnr(rec, image2)));
fprintf(fid, sprintf('PSNR (ref.):\t %g\n', psnr(image1, image2)));
fprintf(fid, sprintf('Elapsed time:\t %g s\n', elapsed));
fclose(fid);

end
