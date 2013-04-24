function cMap = blueredmap(N)

  if nargin < 1
     N = size(get(gcf,'colormap'),1);
  end

  refMap = redbluecmap;
  pos = [0:size(refMap,1)-1];
  pos = pos/pos(end);

  new_pos = [0:N-1];
  new_pos = new_pos/new_pos(end);

  cMap = interp1(pos(:), refMap, new_pos(:));
  cMap = cMap(end:-1:1,:);

  return;
end
