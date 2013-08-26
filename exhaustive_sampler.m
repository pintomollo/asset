function [best_pos, best_score] = exhaustive_sampler(fitfun, p0, opts)

  max_iter = opts.max_iter;
  nparams = length(p0);
  %val_range = log(abs(opts.range_size) + 1);
  val_range = opts.range_size;

  log_file = opts.log_file;
  disp_iter = opts.printint;

  switch (opts.sampling_type)
    case 'independent'
      nevals = ceil((max_iter/nparams)/2);
      pos = [0:(nevals-1)]*(1/(nevals-1));
      if (opts.is_log)
        pos = val_range.^([-pos(end:-1:2) pos(2:end)]);
      else
        pos = 1 + val_range.*([-pos(end:-1:2) pos(2:end)]);
      end
      nevals = length(pos);

      all_pos = ones(nevals*nparams, nparams);
      for i=1:nparams
        all_pos([1:nevals] + (i-1)*nevals,i) = pos(:);
      end
      all_pos(end+1,:) = 1;
      %all_pos = repmat({pos(:)}, 1, nparams);
      %all_pos = enumerate(all_pos{:});
    case 'dual'
      nevals = ceil(sqrt(2*max_iter/(nparams*(nparams-1)))/2);
      pos = [0:(nevals-1)]*(1/(nevals-1));
      if (opts.is_log)
        pos = val_range.^([-pos(end:-1:2) pos(2:end)]);
      else
        pos = 1 + val_range.*([-pos(end:-1:2) pos(2:end)]);
      end

      pos = enumerate(pos(:), pos(:));
      nevals = size(pos, 1);

      all_pos = ones(nevals*(nparams*(nparams-1))/2, nparams);
      count = 0;
      for i=1:nparams-1
        for j=i+1:nparams
          all_pos([1:nevals] + count*nevals,i) = pos(:, 1);
          all_pos([1:nevals] + count*nevals,j) = pos(:, 2);
          count = count + 1;
        end
      end
      all_pos(end+1,:) = 1;
    case 'together'
      nevals = max(ceil((max_iter^(1/nparams))/2),2);
      pos = [0:(nevals-1)]*(1/(nevals-1));
      if (opts.is_log)
        pos = val_range.^([-pos(end:-1:2) pos]);
      else
        pos = 1 + val_range.*([-pos(end:-1:2) pos]);
      end

      all_pos = repmat({pos(:)}, 1, nparams);
      all_pos = enumerate(all_pos{:});
      all_pos(end+1,:) = 1;
  end

  if (val_range == 1)
    all_pos = all_pos(1,:);
    warning('A range of 1 will evaluate only the current position');
  end

  nevals = size(all_pos, 1);
%  disp('Not randomized')

  all_pos = all_pos(randperm(nevals), :);
  p0 = bsxfun(@times, all_pos, p0(:).');

  nulls = all(p0==0,1);
  if (any(nulls))
    tmp_pos = all_pos;
    tmp_pos(tmp_pos < 1) = -1./tmp_pos(tmp_pos < 1);
    tmp_pos(tmp_pos==1) = 0;
    p0(:, nulls) = tmp_pos(:,nulls);
  end

  % Unique identifier for multiple writers in the same file
  uuid = ['EXHSMPLR' num2str(round(rand(1)*100)) ' '];
  do_log = ~(isempty(log_file));

  if (do_log)
    [fid, err] = fopen(['.' filesep log_file '.dat'], 'a');
  end

  best_score = Inf;
  best_pos = NaN(1, nparams);

  fprintf('Exhaustive sampling of [%f %f] in %d (%d x %d) iterations!\n', pos(1), pos(end), nevals, nparams, (nevals-1)/nparams);

  for i=1:nevals
    curr_p = p0(i,:);
    score = fitfun(curr_p);

    if (score < best_score)
      best_score = score;
      best_pos = curr_p;
    end

    if (do_log)
      fprintf(fid, [uuid '%ld : %e |'], i, score);
      fprintf(fid, ' %f', curr_p);
      fprintf(fid, '\n');
    end

    if (mod(i, disp_iter) == 0)
      fprintf('Iter=%d, %d%% done\n', i, round((i/nevals)*100));
    end
  end

  if (do_log)
    fclose(fid);
  end

  return;
end
