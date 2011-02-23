function display_tracking(imgs, spots, links)

  % Adapting to possibly different types of images
  % If it's a string, we assume it's a filename, if it's a structure we assume it's a 'mymovie'
  if (ischar(imgs) |isstruct(imgs))
    % In this case, we store the filename and load the number of frames from it
    fname = imgs;
    nframes = size_data(fname);

  % Otherwise we initialize the filename and treat imgs as a stack
  else
    fname = '';
    nframes = size(imgs, 3);
  end

  for i=1:nframes
    % If we have a file, we load the image
    if (~isempty(fname))
      img = load_data(fname, i);

    % Otherwise we get the plane
    else
      img = imgs(:, :, i);
    end
    
    imshow(img);
    hold on;

    if (i > 1)
      link = links{i};
      %indxs = (link > 0);
      %indxs = find(indxs);

      plot([spots{i-1}(link(:,2), 2) spots{i}(link(:,1), 2)].', [spots{i-1}(link(:,2), 1) spots{i}(link(:,1), 1)].', 'g');

      scatter(spots{i-1}(:,2), spots{i-1}(:,1), 'b');
    end
    scatter(spots{i}(:,2), spots{i}(:,1), 'r');

    hold off;
    input('');
  end

  return;
end
