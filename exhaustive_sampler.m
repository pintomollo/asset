function [best_pos, best_score] = exhaustive_sampler(fitfun, p0, opts)

  max_iter = opts.max_iter;
  nparams = length(p0);
  val_range = log(abs(opts.range_size) + 1);

  log_file = opts.log_file;
  disp_iter = opts.printint;

  nevals = ceil((max_iter^(1/nparams))/2);
  pos = [0:(nevals-1)]*(val_range/(nevals-1));
  pos = exp([-pos(end:-1:2) pos]);

  all_pos = repmat({pos(:)}, 1, nparams);
  all_pos = enumerate(all_pos{:});

  nevals = size(all_pos, 1);
  all_pos = all_pos(randperm(nevals), :);
  all_pos = bsxfun(@times, all_pos, p0(:).');

  % Unique identifier for multiple writers in the same file
  uuid = ['EXHSMPLR' num2str(round(rand(1)*100)) ' '];
  do_log = ~(isempty(log_file));

  if (do_log)
    [fid, err] = fopen(['.' filesep log_file '.dat'], 'a');
  end

  best_score = Inf;
  best_pos = NaN(1, nparams);

  fprintf('Exhaustive sampling of [%f %f] in %d iterations!\n', pos(1), pos(end), nevals);

  for i=1:nevals
    curr_p = all_pos(i,:);
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
