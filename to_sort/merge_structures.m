function [mystruct, is_same] = merge_structures(varargin)
% MERGE_STRUCTURES merges two structures, as well as their sub-structure, by adding
%   the additional fields of the second structure in the first one. It does not
%   combine the non-struct fields while it can concatenate different struct fields
%   if they have different values in the keyword fields.
%
%   MYSTRUCT = MERGE_STRUCTURES(MYSTRUCT, COMPLEMENT) puts all the missing fields
%   of COMPLEMENT into MYSTRUCT, doing the same recursively for the children structure
%   fields.
%
%   MYSTRUCT = MERGE_STRUCTURES(..., STRUCT1, STRUCT2, ...) incrementally merges MYSTRUCT
%   with the provided structures (STRUCT1, STRUCT2, ...).
%
%   MYSTRUCT = MERGE_STRUCTURES(..., KEYWORDS) identify similar structures by comparing
%   the value of their KEYWORD fields. KEYWORD can either be a cell array of field name
%   or individual string name(s), which value is compared to define different structures. 
%   In the case different structures are identified, they are concatenated, while 
%   similare structurs are merged recursively. Empty KEYWORD fields in MYSTRCUT are 
%   replaced by the corresponding value from COMPLEMENT.
%
%   [MYSTRUCT, IS_SAME] = MERGE_STRUCTURES(...) returns IS_SAME, which is true if 
%   MYSTRUCT was not modified during the merge.
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 18.03.2011

  % Prepare the input and output variables
  mystruct = [];
  complement = [];
  keywords = {};
  is_same = true;

  % Used to store the additional structures to be merged with mystruct
  args = {};

  % Extract the input variables from the varargin
  for i = 1:length(varargin)
    % Get the type of input parameter
    type = class(varargin{i});

    % Sort it accordingly
    switch type
      % We assume the first structure is mystruct, the second the complement and
      % all the following ones will be stored in args for later processing
      case 'struct'
        if (isempty(mystruct))
          mystruct = varargin{i};
        elseif (isempty(complement))
          complement = varargin{i};
        else
          args = [args varargin(i)];
        end

      % If it's a boolean, it's the recursive is_same parameter
      case 'logical'
        is_same = is_same & varargin{i};

      % If it's a cell array, it's the keywords
      case 'cell'
        keywords = [keywords varargin{i}];

      % If it's a string, we treat it as a single keyword and attach it to the list
      case 'char'
        keywords = [keywords varargin(i)];
    end
  end

  % If we have only one of the two structures, we return it directly
  if (isempty(mystruct) || isempty(complement))
    if (isempty(mystruct))
      mystruct = complement;

      is_same = false;
    end

    return;
  end

  % Extract the field names of both structures
  destination_fields = fieldnames(mystruct);
  complement_fields = fieldnames(complement);

  % Extract the fields that are in complement but not in mystruct
  missing_fields = ~ismember(destination_fields, complement_fields);

  % Assign an empty value to all of them so we get the full structure
  for field = destination_fields(missing_fields)'
    mystruct = setfield(mystruct, field{1}, []);
  end

  % Check whether there are some keywords fields
  have_keys = ismember(keywords, destination_fields);

  % If there are any, we will need to compare their values
  if (any(have_keys))
    % Retrieve the actual names
    sort_keys = keywords(have_keys);

    % Loop over the different complements (in case of an array of structures)
    for i=1:numel(complement)

      % Initially we have not found any similar structure in mystruct
      twins = -1;

      % Loop over the different mystruct
      for j=1:numel(mystruct)

        % Initialize the variables used to shorten the loops in case of success
        same_structs = true;
        empty_keys = true;

        % Loop over all the keywords
        for k=1:length(sort_keys)

          % Get the type of the current keyword in both structures
          destination_type = class(mystruct(j).(sort_keys{k}));
          complement_type = class(complement(i).(sort_keys{k}));

          % Need to check that at least one field is not empty
          if (empty_keys & ~isempty(mystruct(j).(sort_keys{k})))
            empty_keys = false;
          end

          % We have dissimilar keyword fields if:
          %  - both types are different
          %  - if both strings are not the same (necessary to avoid different string length errors)
          %  - if other types are not equal
          if (~strcmp(destination_type, complement_type) | ...
               (strcmp(destination_type, 'char') & ~strcmp(mystruct(j).(sort_keys{k}), complement(i).(sort_keys{k}))) | ...
               (mystruct(j).(sort_keys{k}) ~= complement(i).(sort_keys{k})))

            % With only one different field, both structures are different so we can stop
            same_structs = false;
            break;
          end
        end

        % If we have not found different fields, we assume they are the same structures
        if (same_structs)
          twins = j;

          break;
        end
      end

      % If we have found a similar structure, we recursively merge them
      if (twins > 0)
        [mystruct(j), is_same] = merge_single_structs(mystruct(j), complement(i), keywords, is_same, destination_fields, complement_fields);

      % Otherwise, if the keyword fields of mystruct are all empty, we replace it
      % by the current complement one (which cannot be empty as they would be equal otherwise)
      elseif (empty_keys)
        mystruct(j) = complement(i);
        is_same = false;

      % Otherwise, we concatenate them
      else
        mystruct = cat(2, mystruct, complement(i));
        is_same = false;
      end
    end

  % If we do not have keywords, we can simply merge them recursively
  elseif (numel(mystruct) == 1 & numel(complement) == 1)

    [mystruct, is_same] = merge_single_structs(mystruct, complement, keywords, is_same, destination_fields, complement_fields);
  end

  % If we have more structures to merge into mystruct, we do it recursively
  if (numel(args) > 0)
    mystruct = merge_structures(mystruct, keywords, is_same, args{:});
  end

  return;
end

% This function merges iteratively two structures, without taking into accound keywords
% as this was done in the main function.
function [mystruct, is_same] = merge_single_structs(mystruct, complement, keywords, is_same, destination_fields, complement_fields)

  % Get the fields of complement which already exist in mystruct
  existing_fields = ismember(complement_fields, destination_fields);

  % If there are some missing fields
  if (any(~existing_fields))
    is_same = false;

    % Loop over the missing fields and copy it from complement
    for field = complement_fields(~existing_fields)'
      field_name = field{1};
      mystruct = setfield(mystruct, field_name, complement.(field_name));
    end
  end

  % Now loop over the existing common fields
  for field = complement_fields(existing_fields)'
    field_name = field{1};

    % If they are structures, we recursively merge them
    if (isstruct(mystruct.(field_name)))
      [mystruct.(field_name), is_same] = merge_structures(mystruct.(field_name), complement.(field_name), keywords, is_same);

    % Otherwise, if the field in complement is a structure and if the field in mystruct
    % is empty, we copy it
    elseif (isempty(mystruct.(field_name)) & isstruct(complement.(field_name)))
      mystruct.(field_name) = complement.(field_name);
      is_same = false;
    end
  end

  return
end
