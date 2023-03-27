function delete(obj)
% Deletes the impoly object
%
% Simon Blanchoud
% University of Fribourg
% 01.10.18
    
  % Check if the object handle has not yet been deleted, if
  % not delete it (the parameter structure is store inside it).
  if (ishandle(obj.hPolygon))
    delete(obj.hPolygon);
  end

  return;
end
