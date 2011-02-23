function [warped, ell_pts] = carth2RECOS(pts, warper, opts, varargin)

  if (isstruct(pts))
    if (isfield(pts, 'breaks'))
      pts = fnval(pts, pts.breaks);
    else
      if (nargin == 1)
        opts = get_struct('RECOS', 1);
      elseif (~isempty(warper))
        opts = warper;
      end
      warped = convert_struct(pts, opts);

      return;
    end
  end

  switch opts.warp_type
    case 'radial'
      if (isempty(warper))
        warped = get_struct('warper',1);

        if (length(varargin) < 3)
          return;
        end

        warped.original.center = varargin{1};
        warped.original.axes_length = varargin{2};
        warped.original.orientation = varargin{3};

        if (length(varargin) >= 6)
          warped.reference.center = varargin{4};
          warped.reference.axes_length = varargin{5};
          warped.reference.orientation = varargin{6};
        end

        warped.warp = elliptic_spline(pts, warped.original.center, warped.original.axes_length, warped.original.orientation);
      else
        ell_pts = carth2elliptic(pts, warper.original.center, warper.original.axes_length, warper.original.orientation);

        if (isempty(ell_pts))
          warped = [];
        else
          corr = fnval(warper.warp, ell_pts(:,1));
          ell_pts(:,2) = ell_pts(:,2) ./ corr;

          warped = elliptic2carth(ell_pts, warper.reference.center, warper.reference.axes_length, warper.reference.orientation);
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
            if (isempty(warpers(j).warp))
              warpers(j) = carth2RECOS(mymovie.(field).eggshell(j).carth, [], opts, mymovie.(field).centers(:,j), mymovie.(field).axes_length(:,j), mymovie.(field).orientations(1,j));
              %first = false;
            end
          
            mymovie.(field).(subfield)(j).warped = carth2RECOS(mymovie.(field).(subfield)(j).carth, warpers(j), opts); 

            if (dv_inversion)
              mymovie.(field).(subfield)(j).warped(:,2) = -mymovie.(field).(subfield)(j).warped(:,2);

              if (strncmp(subfield, 'cortex', 6))
                mymovie.(field).(subfield)(j).warped = mymovie.(field).(subfield)(j).warped(end:-1:1,:);
              end
            end
          end
        end
%        if (isfield(mymovie.(field), 'ruffles'))
%          mymovie.(field).ruffles(j).warped = carth2RECOS(mymovie.(field).ruffles(j).carth, warpers(j), opts); 
%          if (dv_inversion & ~isempty(mymovie.(field).ruffles(j).warped))
%            mymovie.(field).ruffles(j).warped(:,2) = -mymovie.(field).ruffles(j).warped(:,2);
%          end
%        end
%        if (isfield(mymovie.(field), 'centrosomes'))
%          mymovie.(field).centrosomes(j).warped = carth2RECOS(mymovie.(field).centrosomes(j).carth, warpers(j), opts); 
%          if (dv_inversion)
%            mymovie.(field).centrosomes(j).warped(:,2) = -mymovie.(field).centrosomes(i).warped(:,2);
%          end
%        end
      end
      mymovie.(field).warpers = warpers;
    end
  end

  return;
end
