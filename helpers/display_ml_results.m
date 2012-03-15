function p = display_ml_results(fname)

  p = parse_ml_results(fname, Inf, true);

  for i=1:length(p)
    confs = NaN(0, 3);

    for j=1:length(p(i).evolution)
      pts = p(i).evolution{j};
      pts = [[NaN p(i).initial_condition(j, :)]; pts];
      [npts, ngraphs] = size(pts);
      colors = jet(ngraphs);

      x = [-npts:-1].';

      if (j==1)
        figure;hold on;
        all_pts = [pts x ones(size(x))*j];
      else
        all_pts = [all_pts; [[NaN(1, ngraphs+1); [pts x]] ones(length(x)+1, 1)*j]];
      end

      for k=1:ngraphs
        subplot(ceil(sqrt(ngraphs)), ceil(sqrt(ngraphs)),k);hold on;
        plot(x, pts(:,k), 'Color', colors(k,:));
      end

      if (~isempty(p(i).config{j}))
        confs = [confs; [p(i).config{j}{2,2} p(i).config{j}{6,2} j]];
      end
    end
    
    if (~isempty(j))
      obj = mymean(p(i).goal, 1);
      obj = obj .* 10.^(-floor(log10(obj)));

      for k=1:ngraphs
        [m, s, x] = mymean(all_pts(:, k), 1, all_pts(:,end-1));

        subplot(ceil(sqrt(ngraphs)), ceil(sqrt(ngraphs)),k);hold on;
        plot(x, m, 'Color', colors(k,:)*0.5);
        plot(x, m+s, 'Color', colors(k,:)*0.5);
        plot(x, m-s, 'Color', colors(k,:)*0.5);

        if (k>1)
          plot([min(x) max(x)], obj([k-1 k-1]), 'k')
        end
      end

      [v, indx, jndx] = unique(confs(:, 1:2), 'rows')
      for j=1:length(indx)  
        figure;
        indexes = ismember(all_pts(:,end), confs(jndx == j, end));
        for k=1:ngraphs
          subplot(ceil(sqrt(ngraphs)), ceil(sqrt(ngraphs)),k);hold on;
          plot(all_pts(indexes,end-1), all_pts(indexes,k), 'Color', colors(k,:));

          [m, s, x] = mymean(all_pts(indexes, k), 1, all_pts(indexes,end-1));

          plot(x, m, 'Color', colors(k,:)*0.5);
          plot(x, m+s, 'Color', colors(k,:)*0.5);
          plot(x, m-s, 'Color', colors(k,:)*0.5);

          if (k>1)
            plot([min(x) max(x)], obj([k-1 k-1]), 'k')
          end
        end
        mtit(num2str(v(j,:)));
      end
    end
  end

  return;
end
