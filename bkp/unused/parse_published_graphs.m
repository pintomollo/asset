function parse_published_graphs(pweight)

  if (nargin == 0)
    pweight = 0.6;
  end

  vals = textread('quantify_graphs.txt', '%s');
  nvals = length(vals);

  all_data = cell(0,3);
  xscale = zeros(1,2);
  yscale = zeros(1,2);
  origin = zeros(1,2);
  lims = zeros(2,2);

  pts = zeros(0,2);
  stds = zeros(0,2);
  pts_indx = 0;

  data_type = 0;

  for i=1:nvals
    if (isempty(vals{i}))
      continue;
    end

    data = regexp(vals{i},',','split');
    values = str2double(data);

    switch data{1}
      case 'x'
        if (values(4)==0 && values(5)==0)
          xscale = xscale + values(2:3);
        else
          origin(2) = values(3);
          lims(:,1) = [values(2); values(2)+values(4)];
        end
      case 'y'
        if (values(4)==0 && values(5)==0)
          yscale = yscale + values(2:3);
        else
          origin(1) = values(2);
          lims(:,2) = [values(3); values(3)+values(5)];
        end
      case 's'
        pts_indx = pts_indx + 1;
        width = values(4:5)/2;
        center = values(2:3) + width;
        switch stds(pts_indx, 2)
          case -1
            stds(pts_indx,1) = stds(pts_indx,1) + width(2);
            pts(pts_indx,2) = pts(pts_indx,2) + width(2);
          case 1
            stds(pts_indx,1) = stds(pts_indx,1) + width(2);
            pts(pts_indx,2) = pts(pts_indx,2) - width(2);
          case 2
        end
        pts(pts_indx,:) = (pts(pts_indx,:) + center) / 2;

      case 'n'

        if (~isempty(pts))
          resolution = [xscale(1)/xscale(2) yscale(1)/yscale(2)];
          pts = bsxfun(@minus, pts, origin);
          pts = bsxfun(@times, pts, resolution);
          stds = stds(:,1) * resolution(2);

          lims = bsxfun(@minus, lims, origin);
          lims = bsxfun(@times, lims, resolution);

          if (data_type == 1)
            pts(end+1,:) = [lims(2,1) 0];
            stds(end+1) = stds(end);
          end

          all_data{end+1,1} = pts;
          all_data{end,2} = stds;
          all_data{end,3} = lims;
        end

        xscale = zeros(1,2);
        yscale = zeros(1,2);
        origin = zeros(1,2);
        lims = zeros(2,2);

        pts = zeros(0,2);
        stds = zeros(0,2);
        pts_indx = 0;

        data_type = values(2);

      otherwise
        switch values(1)
          case -1
            stds(end+1,:) = [values(5) values(1)];
            pts(end+1,:) = [values(2), values(3) + stds(end,1)];
          case 1
            stds(end+1,:) = [values(5) values(1)];
            pts(end+1,:) = [values(2), values(3)];
          case 2
            stds(end+1,:) = [values(5)/2 values(1)];
            pts(end+1,:) = [values(2), values(3) + stds(end,1)];
        end
    end
  end

  if (~isempty(pts))
    resolution = [xscale(1)/xscale(2) yscale(1)/yscale(2)];
    pts = bsxfun(@minus, pts, origin);
    pts = bsxfun(@times, pts, resolution);
    stds = stds(:,1) * resolution(2);

    lims = bsxfun(@minus, lims, origin);
    lims = bsxfun(@times, lims, resolution);

    if (data_type == 1)
      pts(end+1,:) = [lims(2,1) 0];
      stds(end+1) = stds(end);
    end

    all_data{end+1,1} = pts;
    all_data{end,2} = stds;
    all_data{end,3} = lims;
  end

  for i=1:size(all_data,1)
    pts = all_data{i,1};
    stds = all_data{i,2};
    lims = all_data{i,3};

    pos = [lims(1,1):range(lims(:,1))/100:lims(2,1)];

    smean = csaps(pts(:,1), pts(:,2));
    sstds = csaps(pts(:,1), stds(:,1));

    mval = ppval(smean, pos);
    sval = ppval(sstds, pos);

    figure;errorbar(pts(:,1), pts(:,2), stds);
    xlim(lims(:,1))
    ylim(lims(:,2))

    hold on
    plot(pos, mval, 'k');
    plot(pos, mval+sval, 'r');

    [m,s] = mymean(mval);

    if (abs(m)>1)
      mval = mval / 60;
      [m,s] = mymean(mval);
    end

    if (m > 0)
      %title([max(mval(:)) m s mean(mval(mval>m)) std(mval(mval>m)) mean(mval(mval>m+s)) std(mval(mval>m+s))])
      title([mean(mval(mval>m+s)) std(mval(mval>m+s))])
    else
      %title([min(mval(:)) m s mean(mval(mval<m)) std(mval(mval<m)) mean(mval(mval<m-s)) std(mval(mval<m-s))])
      title([mean(mval(mval<m-s)) std(mval(mval<m-s))])
    end

  end

  keyboard

  return;
end
