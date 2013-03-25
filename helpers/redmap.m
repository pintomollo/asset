function cMap = redmap(N)

  if nargin < 1
     N = size(get(gcf,'colormap'),1);
  end

  refMap = redbluecmap;
  refMap = refMap(((end-1)/2)+1:end,:);
  pos = [0:size(refMap,1)-1];
  pos = pos/pos(end);

  new_pos = [0:N-1];
  new_pos = new_pos/new_pos(end);

  cMap = interp1(pos(:), refMap, new_pos(:));

  return;
end
