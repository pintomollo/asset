function done = save_data(fname, img)

  if(isstruct(fname) & isfield(fname, 'fname'))
    fname = fname.fname;
  end

  done = false;
  waiting_time = 0;
  while (~done)
    try
      imwrite(img, fname, 'TIFF', 'WriteMode', 'append');
      done = true;
    catch ME
      if (waiting_time < 20)
        nsecs = rand(1);
        waiting_time = waiting_time + nsecs;
        pause(nsecs);
      else
        rethrow(ME)
      end
    end
  end

  return
end
