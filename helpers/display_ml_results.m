function p = display_ml_results(fname, file_pattern)

  %close all;
  if (isstruct(fname))
    p = fname;
  else

    if (nargin == 2)
      p = gather_ml_results(fname, file_pattern);
    elseif (strncmp(fname(end-2:end), 'txt', 3))
      p = gather_ml_results(fname);
    else
      p = parse_ml_results(fname, Inf, true);
    end
  end

  if (iscell(p))
    all_pts = [];
    avg_pts = [];
    for i=1:size(p,1)
      data = p{i, 2};
      pts = [];
      for j=1:size(data, 1)
        s = data{j,1}(end);
        if (isstruct(s))
          pts = [pts; [s.params s.score]];
        end
      end

      all_pts = [all_pts; [pts ones(size(pts, 1), 1)*i]];
      avg_pts = [avg_pts; mymean(pts, 1)];
    end

    %subplot(1,2,1)
    %hold off

    for q=1:size(all_pts, 2)-3
      for r=q+1:size(all_pts, 2)-2
        figure;
        scatter(all_pts(:,q), all_pts(:,r), 'b');
        hold on
        scatter(avg_pts(:,q), avg_pts(:,r), 'r');

        xlabel(num2str(q))
        ylabel(num2str(r))
      end
    end

    %figure;
    %scatter(all_pts(:,1), all_pts(:,2), 'b');
    %hold on
    %scatter(avg_pts(:,1), avg_pts(:,2), 'r');
    %subplot(1,2,2)
    %figure;
    %scatter(all_pts(:,1), all_pts(:,3), 'b');
    %hold on
    %scatter(avg_pts(:,1), avg_pts(:,3), 'r');

    %keyboard
  else
    for i=1:length(p)
      confs = NaN(0, 3);

      j = 0;
      for j=1:length(p(i).evolution)
        pts = p(i).evolution{j};

        if (~isempty(p(i).initial_condition))
          pts = [[NaN p(i).initial_condition(j, :)]; pts];
        end
        [npts, ngraphs] = size(pts);
        colors = jet(ngraphs);

        x = [-npts:-1].';

        if (j==1)
        %  figure;hold on;
          all_pts = [pts x ones(size(x))*j];
        else
          all_pts = [all_pts; [[NaN(1, ngraphs+1); [pts x]] ones(length(x)+1, 1)*j]];
        end

        %for k=1:ngraphs
        %  subplot(ceil(sqrt(ngraphs)), ceil(sqrt(ngraphs)),k);hold on;
        %  plot(x, pts(:,k), 'Color', colors(k,:));
        %end

        if (~isempty(p(i).config) & ~isempty(p(i).config{j}))
          confs = [confs; [p(i).config{j}{2,2} p(i).config{j}{6,2} j]];
        end
      end
      %if (j>0)
      %  mtit(['N=' num2str(j)]);
      %end
      
      if (~isempty(j))

        obj = mymean(p(i).goal, 1);
        obj = obj .* 10.^(-floor(log10(obj)));

        %for k=1:ngraphs
        %  [m, s, x] = mymean(all_pts(:, k), 1, all_pts(:,end-1));

        %  h = subplot(ceil(sqrt(ngraphs)), ceil(sqrt(ngraphs)),k, 'Parent', h1, 'NextPlot', 'add'););
        %  plot(h, x, m, 'Color', colors(k,:)*0.5);
        %  plot(h, x, m+s, 'Color', colors(k,:)*0.5);
        %  plot(h, x, m-s, 'Color', colors(k,:)*0.5);

        %  if (k>1)
        %    plot(h,  [min(x) max(x)], obj([k-1 k-1]), 'k')
        %  end
        %end

        [v, indx, jndx] = unique(confs(:, 1:2), 'rows');
        for j=1:length(indx)  
          h1 = figure('Visible', 'off');
          indxs = confs(jndx == j, end);
          indexes = ismember(all_pts(:,end), indxs);
          for k=1:ngraphs
            h = subplot(ceil(sqrt(ngraphs)), ceil(sqrt(ngraphs)),k, 'Parent', h1, 'NextPlot', 'add');

            [m, s, x] = mymean(all_pts(indexes, k), 1, all_pts(indexes,end-1));

            if (k == 1)
              max_thresh = max(m) + 2*max(s);
              thresh = double(all_pts(indexes, k) > max_thresh | isnan(all_pts(indexes, k)));
              thresh = mymean(thresh, 1, all_pts(indexes, end));

              bads = indxs(thresh == 1);
              if (~isempty(bads))
                indxs = indxs(thresh ~= 1);
                indexes = ismember(all_pts(:,end), indxs);
                [m, s, x] = mymean(all_pts(indexes, k), 1, all_pts(indexes,end-1));
                bads = ismember(all_pts(:,end), bads);
              end
            end

            if (any(bads))
              plot(h, all_pts(bads,end-1), all_pts(bads,k), 'Color', [0.5 0.5 0.5]);
            end
            plot(h, all_pts(indexes,end-1), all_pts(indexes,k), 'Color', colors(k,:));
            plot(h, x, m, 'Color', colors(k,:)*0.5);
            plot(h, x, m+s, 'Color', colors(k,:)*0.5);
            plot(h, x, m-s, 'Color', colors(k,:)*0.5);

            if (k>1)
              plot(h, [min(x) max(x)], obj([k-1 k-1]), 'k')
            end
          end
          mtit(h1, [num2str(v(j,:)) ', N=' num2str(length(indxs))]);

          print(h1, '-dpdf', '-r450', ['PNG/fit_synth_data-'  num2str(i) '-' num2str(j) '.pdf']);

          close(h1);
        end
      end
    end
  end

  return;
end
