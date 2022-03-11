%This script compares several approaches to compute prox_{lam*TV}.
%load image
path = fullfile(FAIRpath, 'add-ons','experiments', 'data');
file1 = 'deform2.png';
image1 = double(imresize(imread(fullfile(path, file1)), [128, 128]));
N=size(image1);
name='test';
lambda=10;
noise=0.1;

image2=addnoise(image1,noise);

rec1=proxTV_toolbox(image2,lambda);
rec2=proximalTV(image2,lambda);
rec3=proximalTV1(image2,lambda);

fprintf('toolbox:  SSIM=%.3f.\n', ssim(rec1, image1));
fprintf('indirect: SSIM=%.3f.\n', ssim(rec2, image1));
fprintf('direct: SSIM=%.3f.\n', ssim(rec3, image1));
outputfolder = fullfile(FAIRpath, 'add-ons', 'experiments', 'results','testings','compare denoising');
mkdir(outputfolder);
imwrite(image1 / 255, fullfile(outputfolder, sprintf('%s_source.png', name)));
imwrite(rec1/max(max(rec1(:)), 255) ,  fullfile(outputfolder,sprintf('%s_result_toolbox_%2g_noise_%2g.png', name,lambda,noise)));
imwrite(rec2/max(max(rec2(:)), 255) ,  fullfile(outputfolder,sprintf('%s_result_indirect_%2g_noise_%2g.png', name,lambda,noise)));
imwrite(rec3/max(max(rec3(:)), 255) ,  fullfile(outputfolder,sprintf('%s_result_direct_%2g_noise_%2g.png', name,lambda,noise)));