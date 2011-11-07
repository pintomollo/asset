function display_trackings(trackings, links)

  if (nargin == 1)
    links = {{}};
  end

  nframes = length(trackings);

  if (isstruct(trackings))
    mytracks = trackings;
    trackings = cell(nframes, 1);
    links = cell(nframes, 1);

    for i=1:nframes-1
      trackings{i} = mytracks(i).carth;
      links{i} = mytracks(i).cluster;
    end
    trackings{end} = mytracks(end).carth;
  end

  figure;
  hold on;

  prev_pts = [];
  assigns = [];

  for i=1:nframes
    pts = trackings{i};
    assigns = links{i};

    scatter3(pts(:,1), pts(:,2), i*ones(size(pts(:,1))), 'b');

    if (~isempty(assigns))
      tmp_frames = unique(assigns(:,3)).';

      for j = tmp_frames
        current = (assigns(:, 3) == j);

        plot3([pts(assigns(current,1), 1) trackings{j}(assigns(current,2),1)].', [pts(assigns(current,1), 2) trackings{j}(assigns(current,2),2)].', [i*ones(size(assigns(current, :), 1), 1) assigns(current, 3)].', 'k');
      end
    end
  end
  
  return;
end
