function [ssize, frames, frame_size, h] = size_data(fname)
  
  if(isstruct(fname) & isfield(fname, 'fname'))
    fname = fname.fname;
  end

  h = fopen(fname,'r+');

  frewind(h);
  ssize = fread(h,3,'double')';
  datasize = prod(ssize);
  doublesize = ftell(h) / 3;

  fseek(h,0,'eof');
  frames = ftell(h);
  frames = ((frames / doublesize) - 3) / datasize;

  frame_size = datasize*doublesize;

  if(nargout>3)
    fseek(h,doublesize*3,'bof');
  else
    fclose(h);
  end

  return;
end
