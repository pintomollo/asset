function errors = spline_error(ref, splines, varargin)

  if (numel(splines) == 0)
    errors = NaN;

    return;
  end

  [nbins, polar, opts] = parse_inputs(varargin{:});

  update = [];
  orig_errors = [];
  if (isfield(splines, 'splines'))
    if (isfield(splines, 'update'))
      update = splines.update;
    end
    if (isfield(splines, 'errors'))
      orig_errors = splines.errors;
    end

    splines = splines.splines;
  end

  [orig_t, orig_f, junk, junk] = size(splines);

  reference = [];
  if (polar && isfield(ref, 'elliptic') && isfield(ref, 'reference'))
    reference = ref.reference;
    ref = ref.elliptic; 
  elseif (isfield(ref, 'mean'))
    ref = ref.mean;
  end

  [ntypes, nframes] = size(ref);
  [n, m, o] = size(splines);
  dims = [n m o];

  permdim = perms([3;2;1]);
  dimensions = dims(permdim);

  indx = find(dimensions(:,1) == ntypes & dimensions(:,2) == nframes, 1, 'first');
  splines = permute(splines, permdim(indx,:));

  nsplines = dimensions(indx,3);

  errors = zeros([ntypes nframes nsplines nbins]);
  if (isempty(update))
    update = true([ntypes nframes]);
  end

  if (polar)
    if (mod(splines(1,1,1).dim,2) == 0 && numel(reference) == 0)
      errors = NaN(size(errors));

      return;
    end

    theta = (2*pi) / nbins;
    bins = [0 0.5:nbins nbins] * theta;

    ell_spline = get_struct('spline', 1);

    for t=1:ntypes
      for f=1:nframes
        
        if (~update(t,f) & ~opts.recompute)

          if (t <= orig_t & f <= orig_f)
            errors(t,f,:,:) = orig_errors(t,f,:,:);
          end

          continue;
        end
          
        if (isempty(ref(t,f).breaks))
          errors(t,f,:,:) = NaN;
        else
          if (ref(t,f).dim > 1)
            [junk, ref_knots] = fnplt(ref(t,f));
          else
            ref_knots = fnplt(ref(t,f));
            ref_knots = ref_knots(1,:);
          end

          for i=1:nsplines
            if (isempty(splines(t,f,i).breaks))
              errors(t,f,i,:) = NaN;
            else
              if (numel(reference) ~= 0 & mod(splines(t,f,i).dim,2) == 0)
                ell_spline = carth2elliptic(splines(t,f,i), reference.center(:,f), reference.axes_length(:, f), reference.orientation(1,f));
              else
                ell_spline = splines(t,f,i);
              end

              if (ell_spline.dim > 1)
                [junk, knots] = fnplt(ell_spline);
              else
                knots = fnplt(ell_spline);
                knots = knots(1,:);
              end

              vals = unique([bins, ref_knots, knots]);
              vals = vals(vals >= 0 & vals <= 2*pi);

              ref_pts = fnval(ref(t,f), vals);
              ref_pts = ref_pts(1,:);
              pts = fnval(ell_spline, vals);
              pts = pts(1,:);

              dt = diff(vals);
              dr = ref_pts - pts;
              ds = abs(diff(sign(dr)));
              dr = abs(dr);

              areas = (dr(1:end-1) + dr(2:end)) .* dt ./ 2;

              if (any(ds ~= 0))
                targets = (ds ~= 0);
                denum = dt(targets) ./ (2* (dr([targets false]) + dr([false targets])));
                indxs = ~isinf(denum);

                denum = denum(indxs);
                targets(~indxs) = false;
                
                areas(targets) = (dr([targets false]).^2 + dr([false targets]).^2) .* denum;
              end

              [junk, map] = histc(vals, bins);

              for k=1:nbins
                errors(t,f,i,k) = sum(areas(map == k));
              end
              errors(t,f,i,1) = errors(t,f,i,1) + sum(areas(map(1:end-1) == (nbins + 1)));
            end
          end
          errors(t,f,:,:) = errors(t,f,:,:) / (2*pi);
        end
      end
    end

  else

    thresh = 0.1;

    if (nbins > 1)
      [center , axes_length, orientation] = fit_ellipse(ref(1,1));
      theta = 2*pi / nbins;
      half = theta / 2;
      pie = [0; 0];
      for i=0:pi/8:theta
        pie = [pie [cos(i - half); sin(i - half)]];
      end
      pie = [pie [cos(half) 0; sin(half) 0]];
      pie = pie * 1.5;
      pie = pie(:,[end:-1:1]);
    
      rotmat = [cos(theta) sin(theta); -sin(theta) cos(theta)];
    else
      pie = [];
      theta = 0;
    end

    for t=1:ntypes
      for f=1:nframes
       
        if (isempty(ref(t,f).breaks))
          errors(t,f,:,:) = NaN;
        else
          refpts = fnplt(ref(t,f));
          refarea = polyarea(refpts(1,:),refpts(2,:));

          if (~isempty(pie))
            [c, a, o] = fit_ellipse(ref(t,f));
            ref_pie = ((pie * a(1,1)).' * [cos(o) sin(o); -sin(o) cos(o)]).';
          end

          for i=1:nsplines
            %figure;
            %fnplt(ref(t,f));
            %hold on;
            %fnplt(splines(t,f,i),'g');

            if (isempty(splines(t,f,i).breaks))
              errors(t,f,i) = NaN;
            else

              pts = fnplt(splines(t,f,i));
              for b=1:nbins

                if (~isempty(pie))
                  if (b > 1)
                    ref_pie = (ref_pie.' * rotmat).';
                  end

                  [rx, ry] = polybool('and', refpts(1,:), refpts(2,:), ref_pie(1,:) + c(1), ref_pie(2,:) + c(2));
                  [px, py] = polybool('and', pts(1,:), pts(2,:), ref_pie(1,:) + c(1), ref_pie(2,:) + c(2));
                else
                  rx = refpts(1,:);
                  ry = refpts(2,:);
                  
                  px = pts(1,:);
                  py = pts(2,:);
                end
                refarea = polyarea(rx,ry);
                spline_area = polyarea(px,py);
                [x,y] = polybool('xor',rx,ry,px,py);

                %figure;plot(x,y)
                %[th,r] = cart2pol(x-mean(x),y-mean(y));
                %figure;plot(th)
                %keyboard
                [x,y] = polysplit(x,y);

                refpatch = 0;
                splinepatch = 0;

                %figure;
                %myplot(refpts);
                %hold on;
                %myplot(pts,'Color',[0 1 0]);
                %myplot(ref_pie, 'Color', [ 0 0 0]);
                for j=1:length(x)
                  %patch(x{j},y{j},'r')
                  tmparea = polyarea(x{j},y{j});
                  if (abs(1 - tmparea/spline_area) < thresh & splinepatch == 0)
                    splinepatch = tmparea;
                    %tmparea = -tmparea;
                  elseif (abs(1 - tmparea/refarea) < thresh)
                    refpatch = tmparea;
                  else
                    errors(t,f,i,b) = errors(t,f,i,b) + tmparea;
                  end
                end
                %end
                errors(t,f,i,b) = errors(t,f,i,b) + abs(refpatch - splinepatch);
              end
              errors(t,f,i,:) = errors(t,f,i,:) / splines(t,f,i).breaks(end);
            end
          end
        end
      end
    end
  end

  return;
end

function spline = extract_dim(spline, indx)
  
  pts = fnval(spline, spline.breaks);
  spline = create_spline(pts(indx,:), spline.breaks);

  return;
end

function [nbins, polar, opts] = parse_inputs(varargin)

  nbins = 100;
  polar = true;
  opts = get_struct('RECOS', 1);

  if (nargin > 0)
    for i=1:length(varargin)
      type = get_type(varargin{i});
      switch type
        case 'bool'
          polar = varargin{i};
        case 'num'
          nbins = varargin{i};
        case 'struct'
          opts = varargin{i};
      end
    end
  end

  return;
end
