function mystruct = merge_structures(varargin)

  mystruct = [];
  complement = [];
  keywords = {};
  args = {};
  check_fields = true;

  for i = 1:length(varargin)
    if (isstruct(varargin{i}))
      if (isempty(mystruct))
        mystruct = varargin{i};
      elseif (isempty(complement))
        complement = varargin{i};
      else
        args = [args varargin(i)];
      end
    else
      keywords = varargin{i};
    end
  end

  if (isempty(mystruct) || isempty(complement))
    if (isempty(mystruct))
      mystruct = complement;
    end

    return;
  end

  if (isempty(keywords))
    check_fields = false;
  end

  destination_fields = fieldnames(mystruct);
  complement_fields = fieldnames(complement);

  existing_fields = ismember(complement_fields, destination_fields);
  for field = complement_fields(~existing_fields)'
    field_name = field{1};
    mystruct = setfield(mystruct, field_name, complement.(field_name));
  end

  for field = complement_fields(existing_fields)'
    field_name = field{1};
    if (isstruct(mystruct.(field_name)))
      if (numel(mystruct.(field_name)) == 1 & numel(complement.(field_name)) == 1)
        mystruct.(field_name) = merge_structures(mystruct.(field_name), complement.(field_name), keywords);
      else
        sort_keys = ismember(keywords, fieldnames(mystruct.(field_name)));

        %%% IDEA : Find keyword fields, and if they are all exactly the same in two structures, keep the mystruct one, otherwise, concatenate it.

        if (check_fields & any(sort_keys))
          for i =1:numel(mystruct.(field_name))
            
          end
        else
          [n, m, o] = size(mystruct.(field_name));
          [p, q, r] = size(complement.(field_name));

          if (n == p)
            mystruct.(field_name) = cat(1, mystruct.(field_name), complement.(field_name));
          elseif (m == q)
            mystruct.(field_name) = cat(2, mystruct.(field_name), complement.(field_name));
          elseif (o == r)
            mystruct.(field_name) = cat(3, mystruct.(field_name), complement.(field_name));
          else
            mystruct.(field_name) = cat(1, mystruct.(field_name)(:), complement.(field_name)(:));
          end
        end
      end
    elseif (isempty(mystruct.(field_name)))
      mystruct.(field_name) = complement.(field_name);
    end

    % Merge also all the other types as cells ??
  end

  if (numel(args) > 0)
    mystruct = merge_structures(mystruct, keywords, args{:});
  end

  return;
end
