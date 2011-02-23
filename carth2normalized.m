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
    case 'radial'
      if (isempty(warper))
        warped = get_struct('warper',1);

        if (length(varargin) < 3)
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

        warped.warp = carth2elliptic(pts, warped.original.centers, warped.original.axes_length, warped.original.orientations);
      else
        ell_pts = carth2elliptic(pts, warper.original.centers, warper.original.axes_length, warper.original.orientations);

        if (isempty(ell_pts))
          warped = [];
        else
          %keyboard

          corr = interp_elliptic(warper.warp, ell_pts(:,1));
          ell_pts(:,2) = ell_pts(:,2) ./ corr(:,2);

          warped = elliptic2carth(ell_pts, warper.reference.centers, warper.reference.axes_length, warper.reference.orientations);
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

      nframes = size(mymovie.(field).centers,2);
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
      mymovie.(field).warpers = warpers;
    end
  end

  return;
end
