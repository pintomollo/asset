function [result, r] = load_data(fid, indexes)
% A script for opening microscopy images in MATLAB using Bio-Formats.
%
% The function returns a list of image series; i.e., a cell array of cell
% arrays of (matrix, label) pairs, with each matrix representing a single
% image plane, and each inner list of matrices representing an image
% series. See below for examples of usage.
%
% Portions of this code were adapted from:
% http://www.mathworks.com/support/solutions/en/data/1-2WPAYR/
%
% This method is ~1.5x-2.5x slower than Bio-Formats's command line
% showinf tool (MATLAB 7.0.4.365 R14 SP2 vs. java 1.6.0_20),
% due to overhead from copying arrays.
%
% Thanks to all who offered suggestions and improvements:
%     * Ville Rantanen
%     * Brett Shoelson
%     * Martin Offterdinger
%     * Tony Collins
%     * Cris Luengo
%     * Arnon Lieber
%     * Jimmy Fong
%
% NB: Internet Explorer sometimes erroneously renames the Bio-Formats library
%     to loci_tools.zip. If this happens, rename it back to loci_tools.jar.
%
% Here are some examples of accessing data using the bfopen function:
%
%     % read the data using Bio-Formats
%     data = bfopen('C:/data/experiment.lif');
%
%     % unwrap some specific image planes from the result
%     numSeries = size(data, 1);
%     series1 = data{1, 1};
%     series2 = data{2, 1};
%     series3 = data{3, 1};
%     metadataList = data{1, 2};
%     % ...etc.
%     series1_numPlanes = size(series1, 1);
%     series1_plane1 = series1{1, 1};
%     series1_label1 = series1{1, 2};
%     series1_plane2 = series1{2, 1};
%     series1_label2 = series1{2, 2};
%     series1_plane3 = series1{3, 1};
%     series1_label3 = series1{3, 2};
%     % ...etc.
%
%     % plot the 1st series's 1st image plane in a new figure
%     series1_colorMaps = data{1, 3};
%     figure('Name', series1_label1);
%     if isempty(series1_colorMaps{1})
%         colormap(gray);
%     else
%         colormap(series1_colorMaps{1});
%     end
%     imagesc(series1_plane1);
%
%     % Or if you have the image processing toolbox, you could use:
%     % imshow(series1_plane1, []);
%
%     % Or animate as a movie (assumes 8-bit unsigned data)
%     v = linspace(0, 1, 256)';
%     cmap = [v v v];
%     for p = 1:series1_numPlanes
%         M(p) = im2frame(uint8(series1{p, 1}), cmap);
%     end
%     movie(M);
%
%     % Query some metadata fields (keys are format-dependent)
%     subject = metadataList.get('Subject');
%     title = metadataList.get('Title');

  if (isjava(fid) | nargout == 2)
    [nframes, ssize, pixelType, r] = size_data(fid);

    indexes = indexes(indexes > 0 & indexes <= nframes);

    bpp = loci.formats.FormatTools.getBytesPerPixel(pixelType);
    fp = loci.formats.FormatTools.isFloatingPoint(pixelType);
    sgn = loci.formats.FormatTools.isSigned(pixelType);
    type = char(loci.formats.FormatTools.getPixelTypeString(pixelType));

    bppMax = power(2, bpp * 8);
    little = r.isLittleEndian();
    result = zeros([ssize length(indexes)], type);

    [junk, junk, realEndian] = computer;
    realEndian = strncmp(realEndian, 'L', 1);
    
    for i = 1:length(indexes)
      arr = r.openBytes(indexes(i) - 1);

      if (realEndian ~= little)
        arr = swapbytes(arr);
      end
      arr = loci.common.DataTools.makeDataArray(arr, bpp, fp, little);

      % Java does not have explicitly unsigned data types;
      % hence, we must inform MATLAB when the data is unsigned
      if ~sgn
        switch class(arr)
            case 'int8'
                arr = typecast(arr, 'uint8');
            case 'int16'
                arr = typecast(arr, 'uint16');
            case 'int32'
                arr = typecast(arr, 'uint32');
            case 'int64'
                arr = typecast(arr, 'uint64');
            case 'double'
                arr = typecast(arr, 'double');
            case {'single', 'float'}
                arr = typecast(arr, 'single');
        end
      end

      arr = reshape(arr, ssize([2 1]))';
      result(:,:,i) = arr;
    end

    if (~isjava(fid) & nargout == 1)
      r.close();
    end
  else
    if(isstruct(fid) & isfield(fid, 'fname'))
      fid = fid.fname;
    end

    [nframes, ssize, type] = size_data(fid);

    indexes = indexes(indexes > 0 & indexes <= nframes);

    result = zeros([ssize length(indexes)], type);

    for i = 1:length(indexes)
      result(:,:,i) = imread(fid, indexes(i));
    end
  end
end
