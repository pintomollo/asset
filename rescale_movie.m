function [mymovie] = rescale_movie(mymovie, overwrite)
  
  was_struct = true;
  if (isstruct(mymovie))
    if (isfield(mymovie, 'fname'))
      mymovie = struct('data', mymovie);
      was_struct = false;
    end
  else
    was_struct = false;
    mymovie = struct('data', struct('fname', mymovie, 'file', mymovie, 'detrend', false)); 
  end

  fields = fieldnames(mymovie);

  hwait = waitbar(0,'','Name','CellCoord Info');
  for f = 1:length(fields)
    field = fields{f};

    for k=1:length(mymovie.(field))
      indx = findstr(mymovie.(field)(k).file, '/');
      if (isempty(indx))
        indx = 0;
      else
        indx = indx(end) + 1;
      end
      waitbar(0, hwait, ['Rescaling Movie ' strrep(mymovie.(field)(k).file(indx:end),'_','\_')]);

      new_field = mymovie.(field)(k);
      new_field.fname = '';

      [tmpsize, nframes] = size_data(mymovie.(field)(k));

      if (mymovie.(field)(k).detrend)
        for i=1:nframes
          img = load_data(mymovie.(field)(k),i);
          if (mymovie.(field)(k).hot_pixels)
            img = imhotpixels(img);
          end
          img = imdetrend(img);

          if(overwrite)
            new_field = modify_data(mymovie.(field)(k), img, i);
          else
            new_field = store_data(new_field, img);
          end
          waitbar(i/nframes,hwait);
        end
      else

        for i=1:nframes

          img = load_data(mymovie.(field)(k),i);
          if (mymovie.(field)(k).hot_pixels)
            img = imhotpixels(img);
          end
          minimg = min(min(img));
          maximg = max(max(img));

          if(minimg < mymovie.(field)(k).min)
            mymovie.(field)(k).min = minimg;
          end
          if(maximg > mymovie.(field)(k).max)
            mymovie.(field)(k).max = maximg;
          end

          waitbar(i/(2*nframes),hwait);
        end

        for i=1:nframes
          img = load_data(mymovie.(field)(k),i);
          img = imnorm(img, mymovie.(field)(k).min, mymovie.(field)(k).max);

          if(overwrite)
            new_field = modify_data(mymovie.(field)(k), img, i);
          else
            new_field = store_data(new_field, img);
          end
          waitbar(0.5 + i/(2*nframes),hwait);
        end
      end

      mymovie.(field)(k) = new_field;
    end
  end

  if (~was_struct)
    mymovie = mymovie.data;
  end

  close(hwait);

  clear img;

  return;
end


