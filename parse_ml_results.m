function params = parse_ml_results(fname, nbests, keep_evolution)
  
  if (nargin == 1)
    nbests = 1;
    keep_evolution = false;
  else
    if (islogical(nbests))
      if (nargin == 3)
        [nbests, keep_evolution] = deal(keep_evolution, nbests);
      else
        keep_evolution = nbests;
        nbests = 1;
      end
    elseif (nargin == 2)
      keep_evolution = false;
    end
  end

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
  had_headers = false;
  if (ischar(line) & strncmp(line, '###', 3))
    had_headers = true;
  end

  first = true;

  curr_score = Inf;
  curr_id = NaN;
  curr_params = [];

  if (had_headers)
    curr_conf = {};
    curr_goal = [];
    curr_ic = [];
  end

  if (keep_evolution)
    evol = [];
    nbuffer = 500;
    evol_indx = 0;
  end

  nparams = NaN;
  lparams = NaN;
  while ischar(line)
    if (length(line) > 0)
      tokens = regexp(line, '\s+', 'split');
      ntokens = length(tokens);

      if (had_headers)
        if (ntokens > 0)
          switch tokens{1}
            case 'Config:'
              conf = cell(0, 2);
              for i=2:length(tokens)
                ts = regexp(tokens{i}, '(\w):([\w\.]+),?', 'tokens');
                if (numel(ts) > 0 & length(ts{1}) == 2)
                  ts = ts{1};
                  tmp = str2num(ts{2});
                  if (~isempty(tmp))
                    ts{2} = tmp;
                  end
                  conf(end+1, :) = ts;
                end
              end
            case 'Fitting'
              lparams = str2num(tokens{2});
            case 'IC'
              if (isnan(lparams))
                lparams = 20;
              end
              init_cond = NaN(1, lparams);
              goal = NaN(1, lparams);

              p_indx = 0;
              g_indx = 0;
              for t_indx=2:length(tokens)
                if (tokens{t_indx}(1) == '(')
                  if (p_indx == 0)
                    p_indx = 1;
                    init_cond(p_indx) = str2num(tokens{t_indx}(2:end));
                  else
                    g_indx = 1;
                    goal(g_indx) = str2num(tokens{t_indx}(2:end));
                  end
                elseif (tokens{t_indx}(end) == ')')
                  if (isfinite(p_indx))
                    init_cond(p_indx+1) = str2num(tokens{t_indx}(1:end-1));
                    p_indx = Inf;
                  else
                    goal(g_indx+1) = str2num(tokens{t_indx}(1:end-1));
                    break;
                  end
                else
                  if (isfinite(p_indx))
                    p_indx = p_indx + 1;
                    init_cond(p_indx) = str2num(tokens{t_indx});
                  elseif (g_indx > 0)
                    g_indx = g_indx + 1;
                    goal(g_indx) = str2num(tokens{t_indx});
                  end
                end
              end

              init_cond = init_cond(~isnan(init_cond));
              goal = goal(~isnan(goal));

            case '###'
              init_cond = [];
              goal = [];
              conf = {};

            otherwise
          end
        end
      end

      if (ntokens > 5 && numel(tokens{2}) > 0 && tokens{2}(1) ~= '%')
        ml_id = regexp(tokens{1}, '[a-zA-Z]+([0-9]+)', 'tokens');
        if (~isempty(ml_id))
          ml_id = str2num(ml_id{1}{1});
          
          if (ml_id ~= curr_id)
            if (~isnan(curr_id) & isfinite(curr_score))
              if (isempty(params(nparams).score)|| any(~isfinite(params(nparams).score)) || nbests==1)
                params(nparams).score = curr_score;
                params(nparams).ml_type = curr_type;
                params(nparams).params = curr_params;

                if (had_headers)
                  if (isempty(curr_goal))
                    params(nparams).goal = NaN(1, nparams);
                    params(nparams).initial_condition = NaN(1, nparams);
                  else
                    params(nparams).goal = curr_goal;
                    params(nparams).initial_condition = curr_ic;
                  end
                  params(nparams).config = {curr_conf};
                end

                if (keep_evolution)
                  params(nparams).evolution = {evol(any(~isnan(evol), 2),:)};
                end
              else
                tmp_score = [params(nparams).score; curr_score];
                [score, indexes] = unique(tmp_score);

                if (length(indexes) > nbests)
                  indexes = indexes(1:nbests);
                end
                tmp_ml = [params(nparams).ml_type; curr_type];
                tmp_params = [params(nparams).params; curr_params];

                params(nparams).score = tmp_score(indexes, 1);
                params(nparams).ml_type = tmp_ml(indexes, 1);
                params(nparams).params = tmp_params(indexes, :);

                if (had_headers)
                  if (isempty(curr_goal))
                    tmp_goal = [params(nparams).goal; NaN(1, nparams)];
                    tmp_ic = [params(nparams).initial_condition; NaN(1, nparams)];
                  else
                    tmp_goal = [params(nparams).goal; curr_goal];
                    tmp_ic = [params(nparams).initial_condition; curr_ic];
                  end
                  tmp_conf = [params(nparams).config; {curr_conf}];

                  params(nparams).goal = tmp_goal(indexes,:);
                  params(nparams).initial_condition = tmp_ic(indexes,:);
                  params(nparams).config = tmp_conf(indexes, 1);
                end

                if (keep_evolution)
                  tmp_evol = [params(nparams).evolution; {evol(any(~isnan(evol), 2),:)}];
                  params(nparams).evolution = tmp_evol(indexes, 1);
                end
              end
            end

            curr_id = ml_id;
            curr_score = Inf;
            curr_type = '';
            curr_params(:) = NaN;

            if (keep_evolution)
              evol = [];
              evol_indx = 0;
            end
          end

          if (tokens{1}(1) == 'C')
            iscmaes = true;
            lparams = ntokens - 7;
          else
            iscmaes = false;
            lparams = ntokens - 5;
          end

          score = str2double(tokens{4});

          if (score >= 0)
            if (score < curr_score | keep_evolution)
              if (length(params) < lparams)
                params(lparams) = get_struct('ml_params',1);
              end

              if (iscmaes)
                shift = 6;
              else
                shift = 5; 
              end

              myparams = zeros(1, lparams);
              for j = 1:lparams
                myparams(j) = str2double(tokens{j + shift});
              end

              if (score < curr_score)
                curr_score = score;
                curr_params = myparams;
                curr_type = tokens(1);
                nparams = lparams;

                if (had_headers)
                  curr_conf = conf;
                  curr_goal = goal;
                  curr_ic = init_cond;
                end
              end

              if (keep_evolution)
                evol_indx = evol_indx + 1;

                if (mod(evol_indx, nbuffer) == 1)
                  evol = [evol; NaN(nbuffer, lparams+1)];
                end
                evol(evol_indx, :) = [score myparams];
              end
            end
          elseif (first)
            warning('Found a negative score, make sure this is normal as they are ignored!');
            first = false;
          end
        end
      end
    end

    line = fgetl(fid);
  end

  if (~isnan(curr_id) & isfinite(curr_score))
    if (isempty(params(nparams).score)|| any(~isfinite(params(nparams).score)) || nbests==1)
      params(nparams).score = curr_score;
      params(nparams).ml_type = curr_type;
      params(nparams).params = curr_params;

      if (had_headers)
        if (isempty(curr_goal))
          params(nparams).goal = NaN(1, nparams);
          params(nparams).initial_condition = NaN(1, nparams);
        else
          params(nparams).goal = curr_goal;
          params(nparams).initial_condition = curr_ic;
        end
        params(nparams).config = {curr_conf};
      end

      if (keep_evolution)
        params(nparams).evolution = {evol(any(~isnan(evol), 2),:)};
      end
    else
      tmp_score = [params(nparams).score; curr_score];
      [score, indexes] = unique(tmp_score);

      if (length(indexes) > nbests)
        indexes = indexes(1:nbests);
      end
      tmp_ml = [params(nparams).ml_type; curr_type];
      tmp_params = [params(nparams).params; curr_params];

      params(nparams).score = tmp_score(indexes, 1);
      params(nparams).ml_type = tmp_ml(indexes, 1);
      params(nparams).params = tmp_params(indexes, :);

      if (had_headers)
        if (isempty(curr_goal))
          tmp_goal = [params(nparams).goal; NaN(1, nparams)];
          tmp_ic = [params(nparams).initial_condition; NaN(1, nparams)];
        else
          tmp_goal = [params(nparams).goal; curr_goal];
          tmp_ic = [params(nparams).initial_condition; curr_ic];
        end
        tmp_conf = [params(nparams).config; {curr_conf}];

        params(nparams).goal = tmp_goal(indexes,:);
        params(nparams).initial_condition = tmp_ic(indexes,:);
        params(nparams).config = tmp_conf(indexes, 1);
      end

      if (keep_evolution)
        tmp_evol = [params(nparams).evolution; {evol(any(~isnan(evol), 2),:)}];
        params(nparams).evolution = tmp_evol(indexes, 1);
      end
    end
  end

  fclose(fid);

  return;
end
