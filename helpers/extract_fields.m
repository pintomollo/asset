function [field_values] = extract_fields(mystruct, nframe, type, fields, crop)

  crop_size = NaN(1,2);
  crop_factor = 0;
  if (nargin < 5)
    crop = false;
  elseif (isnumeric(crop))
    if (numel(crop) > 1)
      crop_size = crop;
      crop = true;
    else
      crop_factor = crop;
      crop = true;
    end
  else
    if (crop)
      crop_factor = 2.2;
    end
  end
  if (nargin < 4 | isempty(fields))
    fields = {type, 'eggshell', 'cortex', 'eggshell', 'ruffles', 'centrosomes'; ...
              'img','carth',    'carth',  'axes',     'carth',   'carth'};
  end

  [nargs, nfields] = size(fields);
  field_values = repmat({NaN(1,2)}, nfields, 1);
  field_values(strncmp('centrosomes', fields(1,:), 11)) = {NaN(2,2)};

  if (~isfield(mystruct, type))
    return;
  end

  if (size(fields,1)==1)
    fields = [fields; repmat({'carth'},1,nfields)];
    nargs = nargs + 1;
  end
  
  if (crop | any(strncmp('axes', fields(2,:), 4)))
    if (isfield(mystruct.(type), 'centers'))
      tmp_type = type;
    elseif (isfield(mystruct, 'dic') & isfield(mystruct.dic, 'centers'))
      tmp_type = 'dic';
    else
      return;
    end

    center = mystruct.(tmp_type).centers(:,nframe);
    axes_length = mystruct.(tmp_type).axes_length(:,nframe);
    orient = mystruct.(tmp_type).orientations(:,nframe);
    if (crop_factor > 0)
      crop_size = crop_factor * axes_length.';
      crop_size = round(crop_size(1, [2 1]));
    end
  end

  for i=1:nfields
    try
      switch fields{2,i}
        case 'img'
          if (isfield(mystruct.(fields{1,i}), 'projection'));
            tmp_value = double(load_data(mystruct.(fields{1,i}).projection, nframe));
          else
            tmp_value = double(load_data(mystruct.(fields{1,i}), nframe));
          end
          tmp_value = imnorm(imhotpixels(tmp_value));
        case 'axes'
          tmp_value = bsxfun(@plus, [zeros(1,2); [cos(orient) -sin(orient)] * axes_length(1,1); [-sin(orient) -cos(orient)] * axes_length(2,1)], center');
        case 'carth'
          if (strncmp(fields{1,i}, 'centrosomes', 11))

          tmp_value = mystruct.(type).(fields{1,i})(nframe).(fields{2,i});
          if (nframe > 1)
            links = mystruct.(type).(fields{1,i})(nframe).cluster;
            linked = tmp_value(links(:,1), :);
            tmp_old = mystruct.(type).(fields{1,i})(nframe-1).(fields{2,i});
            linked_old = tmp_old(links(:,2),:);


          else
            linked = [];
            linked_old = [];

          end

          npts = size(linked, 1);
          ntotal = size(tmp_value, 1);
          full_mat = NaN(npts + 2*ntotal, 2);
          indexes = [0:npts-1]*3;
          full_mat(indexes+1,:) = linked;
          full_mat(indexes+2,:) = linked_old;

          remainings = setdiff([1:ntotal].', unique(links(:,1)));
          full_mat((3*npts)+1:2:end,:) = tmp_value(remainings, :);
          tmp_value = full_mat;
        else
          tmp_value = mystruct.(type).(fields{1,i})(nframe).(fields{2,i});
        end
        otherwise
          tmp_value = mystruct.(type).(fields{1,i})(nframe).(fields{2,i});
      end
      if (isstruct(tmp_value) & isfield(tmp_value, 'breaks'))
        tmp_value = fnval(tmp_value, tmp_value.breaks).';
      elseif (iscell(tmp_value))
        tmp_double = [];
        for c=1:numel(tmp_value)
          if (~isempty(tmp_value{c}) & isnumeric(tmp_value{c}))
            tmp_double = [tmp_double; tmp_value{c}];
            tmp_double(end+1,:) = NaN;
          end
        end
        if (isempty(tmp_double))
          continue;
        end
        tmp_value = tmp_double;
      elseif (isempty(tmp_value))
        tmp_value = NaN(1, 2);
      end
      if (crop)
        tmp_value = realign(tmp_value, crop_size, center, orient);
      end
      field_values{i} = tmp_value;
    catch ME
      % Do nothing
      if (~strncmp(ME.identifier,'MATLAB:nonExistentField',23))
        throw(ME);
      end
    end
  end

  %keyboard

  return;
end
