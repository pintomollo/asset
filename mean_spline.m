function [mean_spline, splines] = mean_spline(splines, varargin)
  
  [outside, type] = parse_input(varargin{:});

  [ntypes, nframes, nsplines] = size(splines);

  if (nsplines == 1)
    if (nframes == 1)
      splines = shiftdim(splines, -2);
    else
      splines = shiftdim(splines, -1);
    end

    [ntypes, nframes, nsplines] = size(splines);
  end

  nbreaks = 0;
  
  mins = zeros(ntypes, nframes, nsplines);
  maxs = ones(ntypes, nframes, nsplines);

  center = zeros(ntypes, nframes, nsplines, 2);
  axes_length = zeros(ntypes, nframes, nsplines, 2);
  orientation = NaN(ntypes, nframes, nsplines, 1);

  for t=1:ntypes
    for f=1:nframes
      for i=1:nsplines
        if (isempty(splines(t,f,i).breaks))
          continue;
        end

        if (outside)
          splines(t,f,i) = spline_periphery(splines(t,f,i));
        end

        n = length(splines(t,f,i).breaks);

        if (n > nbreaks)
          nbreaks = n;
        end
        switch type
          case 'radial'
            [center(t,f,i,:), axes_length(t,f,i,:), orientation(t,f,i,1)] = fit_ellipse(splines(t,f,i));
        end

        if (n > 0)  
          a = splines(t,f,i).breaks(1);
          b = splines(t,f,i).breaks(end);

          if (b < a)
            [a,b] = deal(b,a);
          end

          mins(t,f,i) = a;
          maxs(t,f,i) = b;
        end
        
      end
    end
  end

  npts = nbreaks * 4;
  nbins = nbreaks;

  pos = [0:1/(npts):(npts-1)/npts];
  full_pos = [0:1/(nbins-1):1]';
  
  all_mins = min(mins, [], 3);
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
      
      all_mins = zeros(size(all_mins));
      all_maxs = 2 * pi * ones(size(all_maxs));
  end

  %figure;hold on

  for t=1:ntypes
    for f=1:nframes

      all_pts = [];
      %draw_ellipse(center,axes_length,orientation,'r')
      for i=1:nsplines
        %splines(t,f,i) = realign_spline(splines(t,f,i), center, orientation);

        %fnplt(splines(i),'og')
        %fnplt(splines(i),[0 10],'k')
        if (isempty(splines(t,f,i).breaks))
          continue;
        elseif (noks(t,f) == 1)
          mean_spline(t,f) = splines(t,f,i);
          continue;
        end

        pts_pos = mins(t,f,i) + (maxs(t,f,i) - mins(t,f,i))*pos;
        pts = fnval(splines(t,f,i), pts_pos);
        pts = pts(1:2,:);
      %end

      %figure;hold on;
      %for i=1:nsplines
        switch type
          case 'radial'
            pts = carth2elliptic(pts, squeeze(center(t,f,1,:)), squeeze(axes_length(t,f,1,:)), orientation(t,f,1,1));
          case 'linear'
            pts = [pts_pos; pts].';
        end
        %plot(polar(:,1),polar(:,2),'g')
        all_pts = cat(3,all_pts, pts);
      end

      if (isempty(all_pts))
        if (noks(t,f) == 0)
          mean_spline(t,f) = create_spline(all_pts);
        end

        continue;
      else
        nall = size(all_pts, 3);
      end

      tmp_pos = min(all_mins(t,f)) + full_pos * (max(all_maxs(t,f)) - min(all_mins(t,f)));
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
          %std_pts(:,1) = std_pts(:,1) + mean_pts(:,1);

          mean_pts = elliptic2carth(mean_pts, center(t,f,1,:), axes_length(t,f,1,:), orientation(t,f,1,1));
          std_pts = elliptic2carth(std_pts, [0; 0], axes_length(t,f,:), orientation(t,f,1,1));
          %mean_spline(t,f) = create_spline([mean_pts orth]);
          spline_pos = [];

        case 'linear'
          spline_pos = mean_pts(:,1);
          mean_pts = mean_pts(:,2:end);
          std_pts = std_pts(:,2:end);
      end

      %%%% Useful to compute the norm  !!

      %norm = hypot(std_pts(:,1), std_pts(:,2));
      %added = false;
      %if (all(mean_pts(end,:) ~= mean_pts(1,:)))
      %  mean_pts = mean_pts([1:end 1],:);

      %  added = true;
      %end

      %if (size(mean_pts, 1) > 1)
      
      %  dt = sum(diff(mean_pts).^2, 2);
      %  deriv = (mean_pts(2:end,:) - mean_pts(1:end-1,:)) ./ dt(:, [1 1]);
      %  orth = [deriv(:,2) -deriv(:,1)] ./ repmat(hypot(deriv(:,1), deriv(:,2)), 1, 2);

      %  if (added)
      %    mean_pts = mean_pts(1:end-1,:);
      %  else
      %    deriv = deriv([1:end 1],:);
      %    orth = orth([1:end 1],:);
      %  end

      %  orth = orth .* norm(:, [1 1]);
      %else
      %  orth = std_pts;
      %end

      mean_spline(t,f) = create_spline([mean_pts std_pts], spline_pos);
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
