function [resfile,resfile1,resfile2,resfile3,resfile4,paramfile] = saveres(name, outputfolder, image1,...
    image2,image3, def,source, dist, reg, objfun, imageModel, N, nt, alpha, theta,...
    sigma, elapsed, scaleimg, lambda)
%SAVERESULTS A script that saves results and parameters.
% Create string for result filename.
str = sprintf('%s_result%.2f', name, sigma);
str1=sprintf('%s_source%.2f', name, sigma);
str2=sprintf('%s_deformation%.2f', name, sigma);
str3=sprintf('%s_error%.2f', name, sigma);
str4=sprintf('%s_deformation_error%.2f', name, sigma);
resfile1=fullfile(outputfolder, sprintf('%s.png', str1));
resfile2=fullfile(outputfolder, sprintf('%s.png', str2));
resfile3=fullfile(outputfolder, sprintf('%s.png', str3));
resfile4=fullfile(outputfolder, sprintf('%s.png', str4));
resfile = fullfile(outputfolder, sprintf('%s.png', str));
str = sprintf('%s_params%.2f', name, sigma);
paramfile = fullfile(outputfolder, sprintf('%s.txt', str));

rec=def+source;
error=abs(image2-(rec));
error_def=abs(image3-def);

% Save result.
if(scaleimg)
    imwrite(uint8(255 * rec / max(rec(:))), resfile); %save result
    imwrite(uint8(255 * source / max(rec(:))), resfile1); %save source part
    imwrite(uint8(255 * def / max(rec(:))), resfile2); %save source part
    imwrite(uint8(255 * error / max(rec(:))), resfile3); %save error
    imwrite(uint8(255 * error_def / max(rec(:))), resfile4); %save deformation error
else
    imwrite(uint8(rec), resfile); 
    imwrite(uint8(source), resfile1);
    imwrite(uint8(def), resfile2);
    imwrite(uint8(error), resfile3);
    imwrite(uint8(error_def), resfile4);
end


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
fprintf(fid, sprintf('TV-parameter 1:\t %g\n', lambda));
fprintf(fid, sprintf('Angles (deg.):\t %s\n', num2str(theta)));
fprintf(fid, sprintf('Noise level:\t %g\n', sigma));
fprintf(fid, sprintf('SSIM (recon.):\t %g\n', ssim(rec, image2)));
fprintf(fid, sprintf('SSIM (recon. normed):\t %g\n', ssim(rec/max(rec(:)), image2/max(image2(:))) ) );
fprintf(fid, sprintf('SSIM (ref.):\t %g\n', ssim(image1, image2)));
fprintf(fid, sprintf('SSIM (def.):\t %g\n', ssim(def, image3)));
fprintf(fid, sprintf('PSNR (recon.):\t %g\n', psnr(rec, image2)));
fprintf(fid, sprintf('PSNR (ref.):\t %g\n', psnr(image1, image2)));
fprintf(fid, sprintf('PSNR (def.):\t %g\n', psnr(def, image3)));
fprintf(fid, sprintf('max-error:\t %g\n', max(error(:)/255)));
fprintf(fid, sprintf('SSD-error:\t %g\n', norm(error(:)/255)^2/prod(size(error))));
fprintf(fid, sprintf('max-error(def.):\t %g\n', max(error_def(:)/255)));
fprintf(fid, sprintf('SSD-error(def.):\t %g\n', norm(error_def(:)/255)^2/prod(size(error))));
fprintf(fid, sprintf('Elapsed time:\t %g s\n', elapsed));
fclose(fid);
