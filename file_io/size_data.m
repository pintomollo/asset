function [nframes, ssize, infos] = size_data(fname)

    %Substitute imfinfos to handle natively missing files
    if (isempty(fname))
      error('MATLAB:imfinfo:fileOpen', 'No file name was provided.');
    elseif(isstruct(fname) && isfield(fname, 'fname'))
      if (isempty(fname.fname))
        error('MATLAB:imfinfo:fileOpen', 'No file name was provided.');
      else
        fid = fopen(fname.fname, 'r');
      end
    else
      [fid, m] = fopen(fname, 'r');
    end

    if (fid == -1)
      if (isstruct(fname))
        if (isfield(fname, 'eggshell'))
          nframes = length(fname.eggshell);
        elseif (isfield(fname, 'cortex'))
          nframes = length(fname.cortex);
        elseif (isfield(fname, 'centers'))
          nframes = size(fname.centers, 2);
        else
          error('MATLAB:imfinfo:fileOpen', ...
              'Unable to find suitable fields to get the size of the recording.');
        end
      else
        error('MATLAB:imfinfo:fileOpen', ...
            'Unable to open file "%s" for reading.', fname);
      end

      ssize = NaN(1,2);
%      pixelType = 'unknown';
    else
      filename = fopen(fid);  % Get the full pathname if not in pwd.
      fclose(fid);
      
      idx = find(filename == '.');
      format = '';
      if (~isempty(idx))
        extension = lower(filename(idx(end)+1:end));
        % Look up the extension in the file format registry.
        fmt_s = imformats(extension);
        tf = feval(fmt_s.isa, filename);
            
        if (tf)
          format = fmt_s.ext{1};
        end
      end

      if (isempty(format))
        infos = imfinfo(fname);
      else
        infos = feval(fmt_s.info, filename);
      end

      nframes = length(infos);
      ssize = [infos(1).Height infos(1).Width];
      %if (nargout > 2)
      %  type = ['uint' num2str(infos(1).BitDepth)];
      %  switch lower(infos(1).SampleFormat(1))
      %    case 'u'
      %      pixelType = ['uint' num2str(infos(1).BitDepth)];
      %    case 'i'
      %      pixelType = ['int' num2str(infos(1).BitDepth)];
      %    otherwise
%            pixelType = 'double';
      %  end
      %end
    end
%  end

  return;
end
