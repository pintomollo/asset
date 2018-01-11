function [new_str, end_str] = common_substring(str1, str2, replacement_str)
% COMMON_SUBSTRING returns the longest initial and final substrings.
%
%   [INIT_STR, END_STR] = COMMON_SUBSTRING(STR1, STR2) returns the longest common
%   substrings that are located at the begining of both STR1 and STR2 (INIT_STR) and
%   at the end of both (END_STR). In other words, it will return the first N and the last
%   M common characters between both strings.
%
%   [NEW_STRING] = COMMON_SUBSTRING(STR1, STR2, REPL_STR) replaces the non-common central
%   substring by REPL_STR. Separators ('_') are added betweent the common and the 
%   replacement strigns. In other words, it returns the following concatenation:
%   [INIT_STR '_' REPL_STR '_' END_STR].
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 19.05.2011

  % Initialize the output variables
  new_str = '';
  end_str = '';

  % If one string is empty, we have nothing to do
  if (isempty(str1) | isempty(str2))
    return;
  end

  % Get the minimal string size to avoid indexing problems
  nstr = min(length(str1), length(str2));

  % Do a character by character comparison, aligning the strings either with their
  % begining or end, and look for the first/last different one.
  sindx = find(str1(1:nstr)~=str2(1:nstr), 1, 'first');
  eindx = find(str1(end-nstr+1:end)~=str2(end-nstr+1:end), 1, 'last');

  % If we have not found any different character, everything is the same, otherwise go
  % back by one to point on the last identical character, rather than the first different one.
  if (isempty(sindx))
    sindx = nstr;
  else
    sindx = sindx - 1;
  end
  if (isempty(eindx))
    eindx = 1;
  else
    eindx = eindx + 1;
  end
  
  % Extract the common substrings
  new_str = str1(1:sindx);
  end_str = str1(eindx:end);

  % This is a replacement problem
  if (nargin == 3)

    % We add some separators between the common parts and the replacement string,
    % but check if there are some already.
    sep1 = '';
    sep2 = '';
    if (sindx > 0 & str1(sindx) ~= '_')
      sep1 = '_';
    end
    if (eindx < nstr & str1(eindx) ~= '_')
      sep2 = '_';
    end

    % Concatenate the final result
    new_str = [new_str sep1 replacement_str sep2 end_str];
  end

  return;
end
