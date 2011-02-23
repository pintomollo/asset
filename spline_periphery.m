function outside_splines = spline_periphery(mysplines, center, axes_length, orientation, opts)

  update = [];
  if (~isstruct(mysplines))
    mysplines = create_spline(mysplines);
  elseif (~isfield(mysplines,'breaks'))
    if (isfield(mysplines,'update'))
      update = mysplines.update;
    end
    mysplines = mysplines.splines;
  end

  if (nargin == 2)
    opts = center;
  elseif (nargin ~= 5)
    opts = get_struct('RECOS',1);
  end

  [n,m,o] = size(mysplines);
  outside_splines = get_struct('spline',[n,m,o]);
  if (isempty(update))
    update = true([n,m,o]);
  end

  for i=1:n
    for j=1:m
      if (~update(i,j) & ~opts.recompute)
        outside_splines(i,j,:) = mysplines(i,j,:);

        continue;
      end
      for k=1:o
        if (isempty(mysplines(i,j,k).breaks))
          outside_splines(i,j,k) = mysplines(i,j,k);

          continue;
        end

        npts = 4 * length(mysplines(i,j,k).breaks);
        nbins = round(npts / 4);

        if (nargin < 4)
          [center, axes_length, orientation] = fit_ellipse(mysplines(i,j,k));
        end

        full_theta = [0:2*pi/nbins:2*pi]';
        outside_pts = zeros(nbins,2);

        pos = [0:1/(npts):(npts-1)/npts]*mysplines(i,j,k).breaks(end);

        pts = fnval(mysplines(i,j,k), pos);
        pts = pts(1:2,:);

        polar = carth2elliptic(pts,center,axes_length,orientation);

        if (all(diff(polar(:,1))>=0))
          % No detectable invagination
          
          outside_splines(i,j,k) = mysplines(i,j,k);

          continue;
        end

        [junk,ibin] = histc(polar(:,1),full_theta);

        for b=1:nbins
          bin = polar(ibin==b,:);

          indx = find(bin(:,2)==max(bin(:,2)));
          if (length(indx)>1)
            dist = abs(bin(:,1) - (full_theta(b) + pi/nbins));
            indx = find(bin(:,2)==max(bin(:,2)) & dist==min(dist),1);
          end
          if (isempty(bin) | isempty(indx))
            outside_pts(b,:) = NaN;
          else
            outside_pts(b,:) = bin(indx,:);
          end
        end

        outside_pts = outside_pts(~any(isnan(outside_pts),2),:);
        outside_pts = elliptic2carth(outside_pts, center, axes_length, orientation);
        outside_splines(i,j,k) = create_spline(outside_pts);
      end
    end
  end

  %figure;
  %fnplt(myspline);
  %hold on;
  %fnplt(outside_spline,'r')

  return;
end
