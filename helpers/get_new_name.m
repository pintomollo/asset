function new_name = get_new_name(pattern, folder)
% GET_NEW_NAME returns the next available name for a file in an incrementally increasing
% sequence, as defined by the provided regular expression. Returns the relative path.
%
%   [NEW_NAME] = GET_NEW_NAME(PATTERN, FOLDER) uses the provided PATTERN to identify the
%   existing files in FOLDER and returns the next available NEW_NAME by incrementing its
%   "counter". Note that a "counter" is defined as "(\d+)" in a regular expression.
%
%   [...] = GET_NEW_NAME(PATTERN) searches for files in the current folder.
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 19.05.2011

  % Check whether the pattern is valid
  if (findstr(pattern, '*') > 0|isempty(findstr(pattern, '(\d+)')))
    error(['Regular expression ''' pattern ''' cannot be incremented']);
  end

  % Default value for the folder
  if (nargin == 1)
    folder = '.';
  end

  % Get the full path
  cd_dir = pwd;

  % Create the folder if need be
  if (~exist(folder, 'dir'))
    mkdir(folder);
  end

  % Create the full path to the folder
  folder_dir = absolutepath(folder, cd_dir);

  % Initialization
  new_tokens = [];

  % Get the list of files contained in the folder, and parse each of them
  ls_dir = dir(folder_dir);
  for d = 1:length(ls_dir)

    % Check whether the current file matches the provided pattern
    [tmp, tokens] = regexp(ls_dir(d).name, pattern, 'match', 'tokens');

    % Loop over the detected counters and store them
    ntokens = length(tokens);
    for i = 1:ntokens
      new_tokens = [new_tokens str2double(char(tokens{i}))];
    end
  end

  % If we have not found any counter, start at 1
  if (length(new_tokens) == 0)
    indx = '1';

  % Otherwise, we need to figure out which is the first available value
  else
    % Sort the values to identify the "missing" ones
    tokens = sort(new_tokens);

    % Build the list of possible values (including the next one at the end)
    indxs = [1:tokens(end)+1];

    % Extract the ones missing in the folder
    indxs = setxor(tokens, indxs);

    % Get the first missing index
    indx = num2str(indxs(1));
  end

  % Replace the counter pattern with the proper index value
  new_name = regexprep(pattern, '\(.+\)', indx);

  % Remove the regular expression "artifacts"
  new_name = regexprep(new_name, '[\\\?]', '');

  % Build the full filename
  new_name = fullfile(folder_dir, new_name);

  % Get only the relative path
  new_name = relativepath(new_name);

  return;
end
