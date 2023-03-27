function pos = getPosition(obj)
% GETPOSITION returns the list of vertices that composes the polygon.
% 
%   POS = getPosition(POLY)     returns the coordinates POS of the
%                               IMPOLY polygon POLY. Note that the first
%                               vertex is NOT duplicated in closed polygons.
%   POS = POLY.getPosition()    subsref for of the same call.
%
% Simon Blanchoud
% University of Fribourg
% 01.10.18

  % Make sure that the graphic object still exists
  if (~ishandle(obj.hPolygon))
    error('impoly:invalidObject', 'impoly graphic object deleted.');
  end

  % Retrieve the data and the positions
  data = get(obj.hPolygon, 'userdata');
  pos = data.position;

  return;
end
