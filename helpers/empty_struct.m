function is_empty = empty_struct(mystruct, field)
  
  is_empty = true;

  if (isfield(mystruct, field))
    nstruct = numel(mystruct);

    for i=1:nstruct
      if (~isempty(mystruct(i).(field)) & any(isfinite(mystruct(i).(field)(:))))
        is_empty = false;
        break;
      end
    end
  end

  return;
end
