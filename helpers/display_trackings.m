function display_trackings(trackings)

  figure;
  hold on;

  prev_pts = [];
  assigns = [];

  nframes = length(trackings);
  for i=nframes:-1:1
    pts = trackings(i).carth;
    scatter3(pts(:,1), pts(:,2), i*ones(size(pts(:,1))), 'b');

    if (~isempty(assigns))
      plot3([pts(assigns, 1) prev_pts(:,1)].', [pts(assigns, 2) prev_pts(:,2)].', [i*ones(numel(assigns), 1) (i+1)*ones(numel(assigns), 1)].', 'k');
    end

    prev_pts = pts;
    assigns = trackings(i).cluster(:,1);
    prev_pts = prev_pts(assigns ~= 0,:);
    assigns = assigns(assigns ~= 0);
  end
  
  return;
end
