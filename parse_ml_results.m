function params = parse_ml_results(fname)
  
  params = get_struct('ml_params',1);

  if (~exist(fname, 'file'))
    return;
  end

  fid = fopen(fname, 'rt');
  if (fid<0)
    fname = [pwd fname];
    fid = fopen(fname,'rt');
    if (fid<0)
      return;
    end
  end

  line = fgetl(fid);
  while ischar(line)
    if (length(line) > 0)
      tokens = regexp(line, '\s+', 'split');
      ntokens = length(tokens);

      if (ntokens > 5 && tokens{2}(1) ~= '%')
        if (tokens{1}(1) == 'C')
          iscmaes = true;
          nparams = ntokens - 7;
        else
          iscmaes = false;
          nparams = ntokens - 5;
        end

        score = str2double(tokens{4});

        if (length(params) < nparams)
          params(nparams) = get_struct('ml_params',1);
        end

        if (isempty(params(nparams).score) || score < params(nparams).score)
          params(nparams).score = score;
          params(nparams).ml_type = tokens{1};

          if (iscmaes)
            shift = 6;
          else
            shift = 5; 
          end

          myparams = zeros(1, nparams);
          for j = 1:nparams
            myparams(j) = str2double(tokens{j + shift});
          end
          params(nparams).params = myparams;
        end
      end
    end

    line = fgetl(fid);
  end

  fclose(fid);

  return;
end
