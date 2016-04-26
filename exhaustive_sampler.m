function [best_pos, best_score] = exhaustive_sampler(fitfun, p0, opts)

  max_iter = opts.max_iter;
  nparams = length(p0);
  %val_range = log(abs(opts.range_size) + 1);
  val_range = opts.range_size;

  log_file = opts.log_file;
  disp_iter = opts.printint;

  print_str = [' %.' num2str(opts.precision) 'f'];

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
    case {'dual', 'hessian'}
      if (opts.sampling_type(1)=='d')
        nevals = ceil(sqrt(2*max_iter/(nparams*(nparams-1)))/2);
      else
        nevals = 2;
      end

      pos = [0:(nevals-1)]*(1/(nevals-1));
      if (opts.is_log)
        pos = val_range.^([-pos(end:-1:2) pos(2:end)]);
      else
        pos = 1 + val_range.*([-pos(end:-1:2) pos(2:end)]);
      end

      %{
      poz = [0:(2*(nevals-1))]*(1/(2*(nevals-1)));
      if (opts.is_log)
        poz = (2*val_range).^(poz(2:end));
      else
        poz = 1 + (2*val_range).*(poz(2:end));
      end
      poz = enumerate(poz(:), poz(:));
      %}

      single_pos = pos(:);
      pos = enumerate(single_pos, single_pos);
      nevals = size(pos, 1);
      nsingle = length(single_pos);

      all_pos = ones(nevals*(nparams*(nparams-1))/2, nparams);
      single_vals = ones(nsingle*nparams, nparams);
      count = 0;
      for i=1:nparams-1
        for j=i+1:nparams
          all_pos([1:nevals] + count*nevals,i) = pos(:, 1);
          all_pos([1:nevals] + count*nevals,j) = pos(:, 2);
          count = count + 1;
        end
        single_vals([1:nsingle]+(i-1)*nsingle, i) = single_pos;
      end
      single_vals([1:nsingle]+(nparams-1)*nsingle, end) = single_pos;

      all_pos = [all_pos; single_vals];
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
    tmp_pos(tmp_pos < 1) = 1 - 1./tmp_pos(tmp_pos < 1);
    tmp_pos(tmp_pos >= 1) = tmp_pos(tmp_pos>=1) - 1;
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
      fprintf(fid, print_str, curr_p);
      fprintf(fid, '\n');
    end

    if (mod(i, disp_iter) == 0)
      fprintf('Iter=%d, %d%% done\n', i, round((i/nevals)*100));
    end
  end

  if (do_log)
    fclose(fid);
  end

  if (strncmp(opts.sampling_type,'hessian', 7))
    for nloop=1:max_iter
      fprintf('Looks like numerical Hessian estimation, checking that the sampling was adequat... (%d / %d)\n', nloop, max_iter);
      values = group_ml_results(['.' filesep log_file '.dat']);
      values = extract_model_parameters(values);
      [rC, C, rH] = correlation_matrix(values{1,2}{1,2}.evolution);

      precision_thresh = 10^(-opts.precision);

      if (any(abs(rH(:)) < precision_thresh) || any(isnan(rC(:))))
        [I, J] = find((abs(rH) + tril(Inf(nparams), -1) < precision_thresh) | ...
                      (isnan(rC) & triu(ones(nparams), 0)));
        problems = [I J];

        val_range = 2*val_range;
        new_pos = [-1 1].';
        if (opts.is_log)
          new_pos = val_range.^new_pos;
        else
          new_pos = 1 + val_range.*new_pos;
        end

        single_pos = new_pos(:);
        double_pos = enumerate(single_pos, single_pos);

        nsingle = size(single_pos, 1);
        ndouble = size(double_pos, 1);
        is_single = (problems(:,1) == problems(:,2));
        tmp_all = ones(length(single_pos)*sum(is_single) + size(double_pos, 1)*sum(~is_single), nparams);

        count = 0;
        for i=1:length(is_single)
          if (is_single(i))
            tmp_all(count+[1:nsingle], problems(i,1)) = single_pos;

            count = count + nsingle;
          else
            tmp_all(count+[1:ndouble], problems(i,:)) = double_pos;

            count = count + ndouble;
          end
        end

        p0 = median(p0, 1);
        new_p0 = bsxfun(@times, tmp_all, p0(:).');
        nulls = all(new_p0==0,1);

        if (any(nulls))
          tmp_pos = tmp_all;
          tmp_pos(tmp_pos < 1) = 1 - 1./tmp_pos(tmp_pos < 1);
          tmp_pos(tmp_pos >= 1) = tmp_pos(tmp_pos>=1) - 1;
          new_p0(:, nulls) = tmp_pos(:,nulls);
        end

        nevals = size(new_p0, 1);

        fprintf('Resampling %d parameters in %d iterations!\n', size(problems, 1), nevals);

        [fid, err] = fopen(['.' filesep log_file '.dat'], 'a');
        for i=1:nevals
          curr_p = new_p0(i,:);
          score = fitfun(curr_p);

          fprintf(fid, [uuid '%ld : %e |'], i, score);
          fprintf(fid, print_str, curr_p);
          fprintf(fid, '\n');

          if (mod(i, disp_iter) == 0)
            fprintf('Iter=%d, %d%% done\n', i, round((i/nevals)*100));
          end
        end
      else
        break;
      end

      fclose(fid);
    end
  end

  return;
end
