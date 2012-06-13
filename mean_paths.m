function [mean_path] = mean_paths(paths, varargin)

  [outside, type] = parse_input(varargin{:});

  [ntypes, nframes, npaths] = size(paths);

  if (npaths == 1)
    if (nframes == 1)
      paths = shiftdim(paths, -2);
    else
      paths = shiftdim(paths, -1);
    end

    [ntypes, nframes, npaths] = size(paths);
  end
  
  nbreaks = 0;
  maxs = ones(ntypes, nframes, npaths);
  mean_path = cell([ntypes, nframes]);

  switch type
    case 'radial'
      center = zeros(ntypes, nframes, npaths, 2);
      axes_length = zeros(ntypes, nframes, npaths, 2);
      orientation = NaN(ntypes, nframes, npaths, 1);

      maxs = maxs * 2*pi;
    case 'linear'
      dist = cell([ntypes, nframes, npaths]);
  end

  if (outside)
    paths = path_periphery(paths);
  end

  for t=1:ntypes
    for f=1:nframes
      for i=1:npaths
        tmp_path = paths{t,f,i};

        if (numel(tmp_path) == 0)
          continue;
        end

        n = size(tmp_path, 1);
        if (n > nbreaks)
          nbreaks = n;
        end

        switch type
          case 'radial'
            [center(t,f,i,:), axes_length(t,f,i,:), orientation(t,f,i,1)] = fit_ellipse(tmp_path);
          case 'linear'
            tmp_dist = [0; cumsum(sqrt(sum(diff(tmp_path(:,1:2), [], 1).^2, 2)), 1)];
            maxs(t,f,i) = tmp_dist(end);
            dist(t,f,i) = {tmp_dist};
        end
      end
    end
  end

  if (nbreaks == 0) 
    return;
  end

  npts = nbreaks * 4;
  nbins = nbreaks;

  pos = [0:npts-1]' / npts;
  full_pos = [0:nbins-1]' / (nbins-1);
  
  all_maxs = max(maxs, [], 3);
  
  switch type
    case 'radial'
      oks = ~isnan(orientation);
      noks = sum(oks, 3);

      center = sum(center, 3) ./ cat(4, noks, noks);
      axes_length = sum(axes_length, 3) ./ cat(4, noks, noks);

      orientation = align_orientations(orientation, [], 3);
      orientation(~oks) = 0;
      orientation = sum(orientation, 3) ./ noks;
  end

  for t=1:ntypes
    for f=1:nframes
        
      all_pts = [];
      for i=1:npaths
        tmp_path = paths{t,f,i};

        if (isempty(tmp_path))
          continue;
        elseif (noks(t,f) == 1)
          mean_path(t,f) = {tmp_path};
          continue;
        end

        pts_pos = maxs(t,f,i)*pos;

        switch type
          case 'radial'
            tmp_path = carth2elliptic(tmp_path, squeeze(center(t,f,1,:)), squeeze(axes_length(t,f,1,:)), orientation(t,f,1,1));

            try
            pts = interp_elliptic(tmp_path, pts_pos);
            catch
              keyboard
            end
            pts = pts(:,1:2);
          case 'linear'
            pts = interp1q(dist{t,f,i}, tmp_path, pts_pos);
            pts = pts(:,1:2);
            pts = [pts_pos pts];
        end
        all_pts = cat(3,all_pts, pts);
      end

      if (isempty(all_pts))
        continue;
      else
        nall = size(all_pts, 3);
      end

      tmp_pos = full_pos * max(all_maxs(t,f));
      [junk,ibin] = histc(all_pts(:,1,:),tmp_pos);
      mean_pts = zeros(nbins,size(pts,2));
      std_pts = zeros(nbins,size(pts,2));

      for i=1:nbins
        bin = [];
        for j=1:nall
          bin = [bin; all_pts(ibin(:,:,j)==i,:,j)];
        end
        if (isempty(bin))
          mean_pts(i,:) = NaN;
        else
          mean_pts(i,:) = mean(bin, 1);
          std_pts(i,:) = std(bin, [], 1);
        end
      end

      indxs = (~any(isnan(mean_pts),2));
      mean_pts = mean_pts(indxs,:);
      std_pts = std_pts(indxs,:);

      switch type
        case 'radial'
          mean_pts = elliptic2carth(mean_pts, center(t,f,1,:), axes_length(t,f,1,:), orientation(t,f,1,1));
          std_pts = elliptic2carth(std_pts, [0; 0], axes_length(t,f,:), orientation(t,f,1,1));

        case 'linear'
          mean_pts = mean_pts(:,2:end);
          std_pts = std_pts(:,2:end);
      end

      mean_path(t,f) = {[mean_pts std_pts]};
    end
  end
    
  return;
end

function [outside, mytype] = parse_input(varargin)

  outside = false;
  mytype = 'radial';

  if (nargin > 0)
    for i=1:length(varargin)
      type = get_type(varargin{i});
      switch type
        case 'bool'
          outside = varargin{i};
        case 'char'
          mytype = varargin{i};
        case 'struct'
          if (isfield(varargin{i}, 'follow_periphery'))
            outside = varargin{i}.follow_periphery;
          end
      end
    end
  end

  return;
end
