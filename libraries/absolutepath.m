function  abs_path = absolutepath( rel_path, act_path, throwErrorIfFileNotExist )
%ABSOLUTEPATH  returns the absolute path relative to a given startpath.
%
%   Syntax:
%      abs_path = ABSOLUTEPATH( rel_path )
%
%   Parameters:
%      rel_path           - Relative path
%
%
%   See also:  RELATIVEPATH PATH

%   Jochen Lenz

%   Jonathan karr 12/17/2010
%   - making compatible with linux
%   - commented out lower cases
%   - switching findstr to strfind
%   - fixing mlint warnings
%   Jonathan karr 1/11/2011
%   - Per Abel Brown's comments adding optional error checking for absolute path of directories that don't exist
%   Jonathan karr 1/12/2011
%   - fixing bugs and writing test

%   Simon Blanchoud 11/20/2014
%   - Modified to accept already absolute paths
%   - Allowed non-existing absolute paths

%   Simon Blanchoud 11/11/2020
%   - Switched to using Octave build-in functions

%build absolute path
if (is_absolute_filename(rel_path))
    abs_path = rel_path;
else
    abs_path = make_absolute_filename(tilde_expand(rel_path));
end
