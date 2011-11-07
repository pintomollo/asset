function [spots, links] = close_tracking_gaps(spots, gap_links)

  nframes = length(gap_links);
  links = cell(nframes, 1);

  for i=nframes:-1:1
    tmp_links = gap_links{i};

    if (isempty(tmp_links))
      continue;
    end

    prev_spots = (tmp_links(:,3) == i-1);
    links{i} = tmp_links(prev_spots, :);
    tmp_links = tmp_links(~prev_spots, :);

    for j=1:size(tmp_links, 1)
      first_spot = tmp_links(j, :);
      last_spot = gap_links{first_spot(3)}(first_spot(2), :);

      dist = first_spot(3) - last_spot(2);

      dmov = diff([first_spot(1,1:2); last_spot(1,1:2)]);
      interm_pos = first_spot(ones(1,dist-1),1:2) +  ([1:dist-1] / dist).' * dmov;

      for k = 1:dist-1
        frame_index = k+first_spot(3);
        new_pt = ...
      end
    end
  end

  return;
end
