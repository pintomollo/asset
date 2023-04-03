function mystruct = update_structure(mystruct, struct_type)
% UPDATE_STRUCTURE converts an older structure into an up-to-date one.
%
%   [NEW_STRUCT] = UPDATE_STRUCTURE(OLD_STRUCT, STRUCT_TYPE) converts OLD_STRUCT into
%   a NEW_STRUCT of type STRUCT_TYPE, as defined by get_struct.m
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 29.08.2014

  % We use a temporary file to save the current field values
  tmp_fname = absolutepath('./tmp_struct_update.tmp');
  save_parameters(mystruct, tmp_fname);

  % Then we get the new structure and load the values in it
  mystruct = get_struct(struct_type, size(mystruct));
  mystruct = load_parameters(mystruct, tmp_fname);

  % Finally, we remove this mention in the configuration files history, as well as
  % duplicates and empty calls
  if (isfield(mystruct, 'config_files'))
    files = mystruct.config_files(~cellfun('isempty', mystruct.config_files(1:end-1)));
    [values, indx, junk] = unique(files);
    indx = sort(indx);
    mystruct.config_files = files(indx);
  end

  % And delete the temporary file
  delete(tmp_fname);

  return;
end
