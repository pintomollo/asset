function is_closed = getClosed(obj)
% GETCLOSED returns whether the polygon is closed or not.
% 
%   BOOL = getClosed(POLY)      returns BOOL defining whether the
%                               IMPOLY polygon POLY is closed or not.
%   BOOL = POLY.getClosed()     subsref for of the same call.
%
% Simon Blanchoud
% University of Fribourg
% 01.10.18

  % Make sure that the graphic object still exists
  if (~ishandle(obj.hPolygon))
    error('impoly:invalidObject', 'impoly graphic object deleted.');
  end

  % Retrieve the data and the status of the polygon
  data = get(obj.hPolygon, 'userdata');
  is_closed = data.is_closed;

  return;
end
