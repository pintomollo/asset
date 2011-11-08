function [spots, links] = close_tracking_gaps(spots, gap_links)

  nframes = length(gap_links);
  links = cell(nframes, 1);
  new_pt = zeros(1, size(cat(1, spots{:}),2));

  for i=1:nframes
    tmp_links = gap_links{i};

    if (isempty(tmp_links))
      continue;
    end

    prev_spots = (tmp_links(:,3) == i-1);
    links{i} = tmp_links(prev_spots, :);
    tmp_links = tmp_links(~prev_spots, :);

    for j=1:size(tmp_links, 1)
      last_spot = spots{i}(tmp_links(j,1), :);
      first_spot = spots{tmp_links(j,3)}(tmp_links(j,2), :);

      prev_index = tmp_links(j,2);
      dist = i - tmp_links(j,3);

      dmov = diff([first_spot(1,1:2); last_spot(1,1:2)]);
      interm_pos = first_spot(ones(1,dist-1),1:2) +  ([1:dist-1] / dist).' * dmov;

      for k = 1:dist-1
        frame_index = k + tmp_links(j,3);
        new_pt(1:2) = interm_pos(k, :);
        new_index = size(spots{frame_index}, 1) + 1;

        spots{frame_index}(new_index, :) = new_pt;
        links{frame_index}(end+1, :) = [new_index, prev_index, frame_index-1];

        prev_index = new_index;
      end

      links{i}(end+1, :) = [tmp_links(j, 1), prev_index, frame_index];
    end
    links{i} = sortrows(links{i});
  end

  return;
end
