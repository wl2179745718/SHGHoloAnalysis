%%
% This code-file patches together the grid-positioned 3D refractive-index
% patches. Each patch has variables 'rowInd', 'colInd', and 'coordSt',
% which together denote position of patch within a total volume of lateral
% size 'totalFOV_size' x 'totalFOV_size'. All these variables are contained
% within the individual .MAT files.
%
% NOTE: this code is only beta-tested. The user is encouraged to go
% cell-by-cell in this code .M file to confirm that code works as intended.
%
% On a Windows computer, 'cntrl'+'enter' runs a specific cell at a time
%
% Author: Shwetadwip Chowdhury; July 24, 2020

path = 'D:\Shwetadwip\matlab_code\Data\Dataset_08\reconstructions_v1\mat\';
load([path 'recon_910_v1_patch' sprintf('%0.2d',7) '.mat']);
%%
gridNum = 5;
totalFOV =  single(zeros([totalFOV_size,totalFOV_size,size(reconObj,3)]).*NaN);

for stripNum = 1:gridNum
    horStrip =  single(zeros([totalFOV_size,totalFOV_size,size(reconObj,3)]).*NaN);
    for patchIndx = (stripNum-1)*gridNum+[1:gridNum]
        fileName = [path 'recon_910_v1_patch' sprintf('%0.2d',patchIndx) '.mat'];
        if isfile(fileName)
            load(fileName);
            rows = coordSt(rowInd,:);
            cols = coordSt(colInd,:);
            tot1 = single(zeros([totalFOV_size,totalFOV_size,size(reconObj,3)]).*NaN);
            tot1(rows,cols,:) = reconObj(pdar+1:end-pdar,pdar+1:end-pdar,:);
            horStrip = mergePatch(tot1,horStrip);
        end
        fprintf('loaded patch %0.2d \n',patchIndx);
    end
    fprintf('.......loaded strip %0.2d ....... \n',stripNum);
    totalFOV = mergePatch(totalFOV,horStrip);
end
disp('done!');
%%

imagesc(totalFOV(:,:,end/2)); axis equal; axis off; colormap gray; caxis([-0.06,0.14]); colorbar;
%%
% save([folderPath1 'FOV_01_reconstruction_CElegan_totalFOV_notuning.mat'], ...
%     'totalFOV','regParam','prox_iter','t_k','NA', 'lambda', 'ps', 'psz', ...
%     'step_size', 'n_imm', 'n_m', 'pdar','z_plane');