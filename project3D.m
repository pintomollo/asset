function [projection] = project3D(mymovie,planes,opts)
%PROJECT3D Summary of this function goes here
%   Detailed explanation goes here

% Make sure the data have been copied to the current field
if (~isfield(mymovie.data, 'centers') | isempty(mymovie.data.centers) | opts.recompute)
   mymovie = duplicate_segmentation(mymovie, 'data', opts);
   save(mymovie.experiment, 'mymovie', 'trackings','opts');
end

timesteps=size(mymovie.dic.centers,2)/planes

for i = 1:timesteps
   stack(:,:,1:5)=load_data(mymovie.data.fname,i:97:485);
   projection(:,:,i)=max(stack,[],3);
end

