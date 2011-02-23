function [ellpts, radius] = carth2elliptic(ptsx, ptsy, center, axes_length, orient, align)

  if (isempty(ptsx))
    ellpts = [];
    if (nargout == 2)
      radius = [];
    end

    return;
  end

  thresh = 5;
  if (nargin==5 & islogical(orient))
    align = orient;
    thresh = thresh + 1;
  elseif (nargin < 6)
    align = false;
  end

  if (nargin<thresh)
    if (nargin == 1)
      [center, axes_length, orient] = fit_ellipse(ptsx);
    else
      [center, axes_length, orient] = deal(ptsy, center, axes_length);
    end

    if (isstruct(ptsx))
      if (isempty(ptsx.breaks))
        ellpts = ptsx;

        return;
      end
      pts = fnval(ptsx, ptsx.breaks);

      if (nargout == 2)
        [ellpts, radius] = carth2elliptic(pts, center, axes_length, orient)
      else
        if (~ispolycw(pts(1,:), pts(2,:)))
          pts = pts(:,[end:-1:1]);
        end

        ellpts = carth2elliptic(pts, center, axes_length, orient);
        if (isempty(ellpts))
          if (nargout == 2)
            radius = [];
          end

          return;
        end

        min_indx = find(ellpts(:,1) == min(ellpts(:,1)),1);

        if (~isempty(min_indx) && min_indx ~= 1)
          ellpts = ellpts([min_indx:end 1:min_indx-1],:);
        end
        % Add pts on both sides to ensure circularity
        nadd = 10;
        npts = size(ellpts, 1);
        if (nadd > npts)
          nadd = npts;
        end

        ellpts = ellpts([end-nadd+1:end 1:end 1:nadd],:);
        ellpts(1:nadd,1) = ellpts(1:nadd,1) - 2*pi;
        ellpts(end-nadd+1:end,1) = ellpts(end-nadd+1:end,1) + 2*pi;
        ellpts = create_spline(ellpts(:,2:end), ellpts(:,1));
      end

      return;
    end

    if (size(ptsx,2) > 4)
      ptsx = ptsx.';
    end

    if (size(ptsx, 2) == 4)
      [pts_o, pts_r] = carth2elliptic(ptsx(:,1), ptsx(:,2), center, axes_length, orient);
      [stds_o, stds_r] = carth2elliptic(ptsx(:,3), ptsx(:,4), [0; 0], axes_length, orient);

      if (nargout == 2)
        ellpts = [pts_o stds_o];
        radius = [pts_r stds_r];
      else
        ellpts = [pts_o pts_r stds_o stds_r];
      end

      return;
    end

    ptsy = ptsx(:,2);
    ptsx = ptsx(:,1);
  else
    ptsx = ptsx(:);
    ptsy = ptsy(:);
  end

  ellpts = zeros(length(ptsx),2);

  x = ptsx - center(1);
  y = -(ptsy - center(2));

  ratio = axes_length(1) / axes_length(2);

  corient = cos(orient); 
  sorient = sin(orient); 

  long_axe = (x*corient + y*sorient)/axes_length(1);
  short_axe = (y*corient - x*sorient)/axes_length(2);
  
  ellpts(:,1) = atan2(short_axe, long_axe);
  r = hypot(long_axe, short_axe);
  %r = sqrt((long_axe).^2 + (short_axe).^2);

  negs = ellpts(:,1)<0;
  ellpts(negs,1) = ellpts(negs,1) + 2*pi;
  ellpts(:,2) = r;

  % Correct a rounding error
  if (ellpts(1,1)==2*pi)
    ellpts(1,1) = 0;
  end

  if (align)
    indx = find(ellpts(:,1)==min(ellpts(:,1)));

    if (all(ellpts(1,:)) == ellpts(end,:))
      ellpts = [ellpts(indx:end-1,:); ellpts(1:indx-1,:); ellpts(indx,:)];
    else
      ellpts = [ellpts(indx:end,:); ellpts(1:indx-1,:)];
    end
  end

  if(nargout==2)
    radius = ellpts(:,2);
    ellpts(:,2) = [];
  end

  return;
end
