function migrate_data(orig_folder, new_folder, duplicate)

  if (nargin == 1)
    new_folder = '.';
    duplicate = [];
  elseif (nargin == 2)
    if (islogical(new_folder))
      duplicate = new_folder;
      new_folder = '.';
    else
      duplicate = [];
    end
  end

  if (isempty(duplicate))
    answer = questdlg('Do you want to keep the original files ?', 'Data Migration');
    switch answer
      case 'Yes'
        duplicate = true;
      case 'No'
        duplicate = false;
      case 'Cancel'
        return;
    end
  end

  orig_folder = absolutepath(orig_folder);
  new_folder = absolutepath(new_folder);

  if (exist(orig_folder, 'dir') ~= 7)
    error([orig_folder ' is not a folder']);
  end
  if (exist(new_folder, 'dir') ~= 7)
    error([new_folder ' is not a folder']);
  end

  files = dir(fullfile(orig_folder, '*.mat'));
  tmp_folder = fullfile(orig_folder, 'TmpData');
  new_tmp_folder = fullfile(new_folder, 'TmpData');

  migrate_data = (exist(tmp_folder, 'dir') == 7);

  for i=1:length(files)
    changed = false;
    fname = fullfile(orig_folder, files(i).name);
    new_fname = fullfile(new_folder, files(i).name);
    data = load(fname);

    if (isfield(data, 'mymovie'))
      fields = fieldnames(data.mymovie);

      if (exist(new_fname, 'file') == 2)
        new_data = load(new_fname);

        for f = 1:length(fields)
          field = fields{f};

          if (~isempty(data.mymovie.(field)) & isfield(data.mymovie.(field), 'fname'))
            data.mymovie.(field).fname = new_data.mymovie.(field).fname;
          end
          if (~isempty(data.mymovie.(field)) & isfield(data.mymovie.(field), 'file'))
            data.mymovie.(field).file = new_data.mymovie.(field).file;
          end
        end
      else
        for f = 1:length(fields)
          field = fields{f};

          if (~isempty(data.mymovie.(field)) & isfield(data.mymovie.(field), 'fname'))
            indx = regexp(data.mymovie.(field).fname, '[/\\]');
            if (~isempty(indx))
              name = data.mymovie.(field).fname(indx(end)+1:end);
            else
              name = fname;
            end
            orig_file = fullfile(tmp_folder, name);
            new_file = get_new_name('tmpmat(\d+)\.ome\.tiff?', new_tmp_folder);
            if (exist(orig_file, 'file') == 2)

              if (duplicate)
                copyfile(orig_file, new_file, 'f');
              else
                movefile(orig_file, new_file, 'f');
              end

              data.mymovie.(field).fname = new_file;

              data.mymovie.(field).file = strrep(data.mymovie.(field).file, '\', filesep);
              data.mymovie.(field).file = strrep(data.mymovie.(field).file, '/', filesep);

              changed = true;
            end
          end
        end
      end
    end

    mymovie = data.mymovie;
    opts = data.opts;
    if (isfield(data, 'trackings'))
      trackings = data.trackings;
    else
      trackings = [];
    end
  
    if (duplicate)
      save(new_fname, 'mymovie', 'trackings', 'opts');
    else
      save(fname, 'mymovie', 'trackings', 'opts');
      movefile(fname, new_folder, 'f');
    end
    
    disp(new_fname);
  end

  return;
end
