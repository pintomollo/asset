function [warped, ell_pts] = carth2normalized(pts, warper, opts, varargin)

  %keyboard

  if (isstruct(pts))
    if (nargin == 1)
      opts = get_struct('ASSET', 1);
    elseif (~isempty(warper))
      opts = warper;
    end
    warped = convert_struct(pts, opts);

    return;
  end

  switch opts.warp_type
    case {'radial', 'hyperbolic'}
      if (isempty(warper))
        warped = get_struct('warper',1);

        if (length(varargin) < 3 |isempty(pts))
          return;
        end

        warped.original.centers = varargin{1};
        warped.original.axes_length = varargin{2};
        warped.original.orientations = varargin{3};

        if (length(varargin) >= 6)
          warped.reference.centers = varargin{4};
          warped.reference.axes_length = varargin{5};
          warped.reference.orientations = varargin{6};
        end

        warp = carth2elliptic(pts, warped.original.centers, warped.original.axes_length, warped.original.orientations, opts.warp_type);
        warped.warp = linearize_egg(warp);
      else
        if (isempty(pts))
          warped = [];
          ell_pts = pts;
          return;
        end

        ell_pts = carth2elliptic(pts, warper.original.centers, warper.original.axes_length, warper.original.orientations, opts.warp_type);

        if (isempty(ell_pts))
          warped = [];
        else
          [ordered, index] = sort(ell_pts(:,1));
          [junk, index] = sort(index); 

          %keyboard

          try
          corr = interp_elliptic(warper.warp, ordered);
          ell_pts(:,2) = ell_pts(:,2) ./ corr(index,2);
          catch
            beep;keyboard;
          end

          warped = elliptic2carth(ell_pts, warper.reference.centers, warper.reference.axes_length, warper.reference.orientations, opts.warp_type);
        end
      end
  end

  return;
end

function mymovie = convert_struct(mymovie, opts)

  %keyboard

  fields = fieldnames(mymovie);
  for i=1:length(fields)
    if (~isempty(mymovie.(fields{i})) & isfield(mymovie.(fields{i}), 'centers'))
      field = fields{i};

      dv_inversion = false;
      if (isfield(mymovie.(field), 'dv_inverted') & mymovie.(field).dv_inverted)
        dv_inversion = true;
      end

      if (isfield(mymovie.(field), 'fname'))
        nframes = size_data(mymovie.(field));
      elseif (strncmp(field, 'markers', 7) && isfield(mymovie, 'cortex') && isfield(mymovie.cortex, 'fname'))
        nframes = size_data(mymovie.cortex); 
      else
        nframes = size(mymovie.(field).centers,2);
      end
      warpers = get_struct('warper',[1 nframes]);
      
      subfields = fieldnames(mymovie.(field));
      %first = true;
      for k=1:length(subfields)
        if (~(isempty(mymovie.(field).(subfields{k}))) & isfield(mymovie.(field).(subfields{k}), 'carth') & ~strncmp(subfields{k}, 'eggshell', 8))
          subfield = subfields{k};
          for j = 1:nframes
            if (isempty(warpers(j).warp) | opts.recompute)
              warpers(j) = carth2normalized(mymovie.(field).eggshell(j).carth, [], opts, mymovie.(field).centers(:,j), mymovie.(field).axes_length(:,j), mymovie.(field).orientations(1,j));
              %first = false;
            end

            if (any(~isfinite(warpers(j).original.axes_length)))
              mymovie.(field).(subfield)(j).warped = [];
            else
              mymovie.(field).(subfield)(j).warped = carth2normalized(mymovie.(field).(subfield)(j).carth, warpers(j), opts); 

              if (dv_inversion)
                mymovie.(field).(subfield)(j).warped(:,2) = -mymovie.(field).(subfield)(j).warped(:,2);

                if (strncmp(subfield, 'cortex', 6))
                  mymovie.(field).(subfield)(j).warped = mymovie.(field).(subfield)(j).warped(end:-1:1,:);
                end
              end
            end
          end
        end
      end
      mymovie.(field).warpers = warpers;
    end
  end

  return;
end

function pts = linearize_egg(pts)

  theta = pts(:,1);
  dt = diff(theta);
  step = mean(abs(dt));
  
  negs = (dt < 0);
  smalls = (abs(dt) < 2*step) & negs;
  others = ~smalls & negs;
  
  pts = pts(~smalls, :);
  others = others(~smalls);

  if (sum(others) > 1)
    error('Invalid eggshell as it is not linearly increasing')
  else
    indx = find(others);
    if (~isempty(indx))
      pts = [pts(indx+1:end,:); pts(1:indx,:)];
    end
  end

  return;
end
