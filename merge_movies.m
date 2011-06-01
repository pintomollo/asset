function mymovie = merge_movies(mymovies)

  mymovie = get_struct('mymovie',1);
  trackings = get_struct('trackings',1);
  
  if (isempty(mymovies))
    return;
  end

  nframes = 0;
  for i=1:length(mymovies)
    current_movie = load(mymovies{i});
    [mymovie, nframes] = merge_recursive(mymovie, current_movie.mymovie, nframes);
    if (isfield(current_movie, 'trackings'))
      trackings = merge_recursive(trackings, current_movie.trackings, nframes);
    end
  end

  save(mymovie.experiment, 'mymovie', 'trackings');

  return;
end

function [mystruct, new_nframes] = merge_recursive(mystruct, twin_struct, nframes)

  new_nframes = nframes;

  fields = fieldnames(twin_struct);
  for i = 1:length(fields)
    field = fields{i};
    if (isfield(mystruct, field) & ~isempty(twin_struct.(field)))
      if (isempty(mystruct.(field)))
        mystruct.(field) = twin_struct.(field);
      else
        issingle = (prod(size(twin_struct.(field))) == 1);
          
        switch (get_type(twin_struct.(field)))
          case 'struct'
            switch field
              case {'child', 'files'}
                mystruct.(field) = match_struct(mystruct.(field), twin_struct.(field), nframes);
              otherwise
                if (issingle)
                  [mystruct.(field), tmp_nframes] = merge_recursive(mystruct.(field), twin_struct.(field), nframes);
                  if (tmp_nframes > new_nframes)
                    new_nframes = tmp_nframes;
                  end
                else
                  mystruct.(field) = mycat(mystruct.(field), twin_struct.(field));
                end
            end
          case 'char'
            switch field
              % mymovie
              case 'fname'
                [mystruct.(field), new_size] = append_file(mystruct.(field), twin_struct.(field), nframes);
                if (new_size > 0)
                  new_nframes = new_size;
                end
              case 'file'
                if (iscell(mystruct.(field)))
                  mystruct.(field) = [mystruct.(field), twin_struct.(field)];
                else
                  mystruct.(field) = {mystruct.(field), twin_struct.(field)};
                end

              % trackings
              case {'expr', 'name', 'experiment'}
                mystruct.(field) = common_substring(mystruct.(field), twin_struct.(field), 'MERGE');
            end
          case 'num'
            if (~issingle && ~strncmp(field, 'color', 5))
              mystruct.(field) = mycat(mystruct.(field), twin_struct.(field));
            end
          case 'java'
            mystruct.(field) = combine_metadata(mystruct.(field), twin_struct.(field));
        end
      end
    end
  end

  return;
end

function [fname, new_size] = append_file(new_file, old_file, frame_shift)
  
  [fpath, name, ext, vers] = fileparts(old_file);
  new_size = 0;

  switch ext
    case '.shapes'
      fname = common_substring(new_file, old_file, 'MERGE');

      if (~strncmp(fname, new_file, length(new_file)))
        copyfile(new_file, fname, 'f');
      end

      fold = fopen(old_file, 'rt');
      fid = fopen(fname, 'a+');
      fprintf(fid, '\nFrameShift=%d\n', frame_shift);
  
      line = fgets(fold);
      while ischar(line)
        fwrite(fid, line);
        line = fgets(fold);
      end

      fclose(fid);
      fclose(fold);
    case '.tmp'
      if (frame_shift == 0)
        fname = '';
        [nframes, imgsize] = size_data(new_file);
        for i=1:nframes
          fname = store_data(fname, load_data(new_file, i));
        end
      else
        fname = new_file;
      end

      [new_size, imgsize] = size_data(fname);

      [nframes, imgsize] = size_data(old_file);
      for i=1:nframes
        fname = store_data(fname, load_data(old_file, i));
      end
  end

  return;
end

function mystruct = match_struct(mystruct, twin_struct, nframes)

  for i=1:length(twin_struct)

    found = false;
    if (isfield(twin_struct(i), 'name'))
      name = regexp(twin_struct(i).name, '([^-]+).*', 'tokens');
      if (~isempty(name))
        name = name{1};
      end

      for j=1:length(mystruct)
        new_name = regexp(mystruct(j).name, '([^-]+).*', 'tokens');
        if (~isempty(new_name))
          new_name = new_name{1};
        end

        if (strcmp(new_name, name))
          [mystruct(j), new_size] = merge_recursive(mystruct(j), twin_struct(i), nframes);
          found = true;

          break;
        end
      end

    else
      indx = regexp(twin_struct(i).fname, '.*?-(\d+)\.shapes', 'tokens');
      indx  = str2double(indx{1}{1});

      for j=1:length(mystruct)
        new_indx = regexp(mystruct(j).fname, '.*?-(\d+)\.shapes', 'tokens');
        new_indx  = str2double(new_indx{1}{1});

        if (new_indx == indx)
          [mystruct(j), new_size] = merge_recursive(mystruct(j), twin_struct(i), nframes);

          found = true;
          break;
        end
      end
    end

    if (~found)
      mystruct(end+1) = duplicate(twin_struct(i), nframes);
    end
  end

  return;
end

function old_struct = duplicate(old_struct, frame_shift)

  fields = fieldnames(old_struct);

  for i = 1:length(fields)
    field = fields{i};

    if (~isempty(old_struct.(field)))
      [m,n,o,p] = size(old_struct.(field));
      switch (get_type(old_struct.(field)))
        case 'struct'
          switch field
            case {'mean', 'paths'}
              buffer = cell([m frame_shift o p]);
              old_struct.(field) = mycat(buffer, old_struct.(field));
            case 'shapes'
              buffer = repmat(struct('path',[]), [m frame_shift o p]);
              old_struct.(field) = mycat(buffer, old_struct.(field));
            case {'files', 'child'}
              for j=1:length(old_struct.(field))
                old_struct.(field)(j) = duplicate(old_struct.(field)(j), frame_shift);
              end
          end
        case 'num'
          buffer = NaN([m frame_shift o p]);
          old_struct.(field) = mycat(buffer, old_struct.(field));
        case 'char'
          switch field
            case 'fname'
              new_name = regexprep(old_struct.(field), '-\d+-','_MERGE_');
              if (exist(new_name, 'file'))
                delete(new_name);
              end
              old_struct.(field) = append_file(new_name, old_struct.(field), frame_shift);
            otherwise
              old_struct.(field) = regexprep(old_struct.(field), '-\d+-','_MERGE_');
          end
      end
    end
  end

  return;
end

function merge = mycat(array1, array2)

  [a, b, c, d] = size(array1);
  [m, n, o, p] = size(array2);

  if (a ~= m)
    merge = array1;

    return;
  end

  nbuff = c - o;

  if (nbuff ~= 0)
    switch get_type(array1)
      case 'num'
        buffer = NaN;
      case 'struct'
        buffer = array1(1,1,1);
      case 'cell'
        buffer = cell(1);
    end

    if (nbuff > 0)
      array2 = cat(3, array2, repmat(buffer, [m n abs(nbuff) p]));
    else
      array1 = cat(3, array1, repmat(buffer, [a b abs(nbuff) d]));
    end
  end

  try
    merge = [array1 array2];
  catch ME
    if (strncmp(ME.identifier, 'MATLAB:catenate:structFieldBad',30))
      merge = [array1 rmfield(array2, setdiff(fieldnames(array2), fieldnames(array1)))];
    else
      throw(ME);
    end
  end

  return;
end

function meta = combine_metadata(meta, new_meta)

  meta = char(meta);
  new_meta = char(new_meta);
  [starts, ends] = regexp(meta,'<Plane.*?>.*?</\Plane>');
  ntimes = length(starts);

  if (ntimes == 0)
    return;
  end

  insert = ends(end);

  [starts, ends] = regexp(new_meta,'<Plane.*?>.*?</\Plane>');
  substr = new_meta(starts(1):ends(end));
  ntimes = ntimes + length(starts);

  meta = [meta(1:insert) substr meta(insert+1:end)];
  meta = regexprep(meta, 'SizeT="\d+"', ['SizeT="' num2str(ntimes) '"']);
  meta = java.lang.String(meta);

  return;
end
