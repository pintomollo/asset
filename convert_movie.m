function [newfile, policy] = convert_movie(fname, policy, opts)
% CONVERT_MOVIE converts (potentially) any type of movie into OME-TIFF using the
% Bio-Formats library (http://www.loci.wisc.edu/software/bio-formats).
%
%   [NEW_MOVIE] = CONVERT_MOVIE(FNAME, OPTS) converts the movie FNAME into the OME-TIFF
%   type and returns its new file name (NEW_MOVIE).
%
%   [...] = CONVERT_MOVIE(FNAME, POLICY, OPTS) provides in addition the action to be
%   performed in case the OME-TIFF version of FNAME already exists.
%   The possible actions are:
%     0   Ask the user
%     1   Overwrite
%     2   Overwrite for all the files
%     3   Use the existing one
%     4   Use all the existing ones
%
%   [NEW_MOVIE, POLICY] = CONVERT_MOVIE(...) returns in addition the policy. Useful
%   when looping over various channels of the same movie to remember the user's choice.
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 19.05.2011

  % Initialization
  newfile = '';

  % Create a null policy if it is not provided
  if (nargin == 2)
    opts = policy;
    policy = 0;
  end

  % If there is no movie file, stop here
  if(isempty(fname) | ~ischar(fname))
    return;
  end

  % Import the Java classes
  import loci.formats.ImageReader;
  import loci.formats.out.OMETiffWriter;
  import loci.formats.gui.LegacyQTTools;
  import loci.formats.MetadataTools;
  import loci.formats.ChannelMerger;

  % We need the absolute path for Java to wrk properly
  fname = absolutepath(fname);

  % Create the Metadata structure necessary to copy them from the movie
  omexmlMeta = MetadataTools.createOMEXMLMetadata();

  % Create the file reader pointing to our file
  r = ImageReader();

  % If required, try to merge separated files
  if (opts.merge_input_files)
    r = ChannelMerger(r);
  end

  % Set the metadata structure, this will copy into it the metadata
  r.setMetadataStore(omexmlMeta);

  % Split the filename using the provided pattern to change its extension to OME-TIFF
  [tokens,junk]=regexp(fname, opts.file_regexpr, 'tokens');
  name = tokens{1}{1};
  suffix = tokens{1}{2};
  ext = tokens{1}{3};

  % Remove the extension
  name = [name suffix];

  % Identify the filename VS the path
  [slash] = findstr(name, filesep);
  if(length(slash)>0)
    slash = slash(end) + 1;
  else
    slash = 1;
  end

  % Creat the fancy name for display (otherwise it thinks they are LaTeX commands)
  printname = strrep(name(slash:end),'_','\_');

  % Initialize the movie
  r.setId(fname);

  % Get the current file format
  format = char(r.getReader().getFormat());

  % We have to particular actions: if it's QuickTime or if it's already OME-TIFF
  if (strncmp(format,'QuickTime',9))

    % In case it's a QT movie, check if the legacy is usable
    qt = LegacyQTTools;
    qtjava = qt.canDoQT();

    % If it is, use it
    if (qtjava)
      r.getReader().setLegacy(true);

    % Otherwise, display why it is not
    else
      qt.checkQTLibrary();
    end
  elseif (strncmp(format,'OME-TIFF',8))

    % If it's already an OME-TIFF file, we do not need to convert it, so stop here
    newfile = fname;
    newfile = relativepath(newfile);

    % Close the file reader
    r.close();

    return;
  end

  % Create the file writer
  w = OMETiffWriter();

  % Copy the metadata into it
  w.setMetadataRetrieve(omexmlMeta);

  % Should speed-up the writing process
  w.setWriteSequentially(true);

  % Check whether the conversion is feasible
  if(~w.isSupportedType(r.getPixelType()))

    r.close();
    error(['Conversion from ' format ' to OME-TIFF is not feasible !']);

    return
  end

  % We create an OME-TIFF file
  newname = [name '.ome.tiff'];

  % If the file already exists, we have several possibilities (c.f. help section)
  if(exist(newname,'file'))

    % We initially do not know what to do
    answer = 0;

    % If there is no provided policy, we need to ask the user
    if (policy == 0)

      % We do not accept "empty" answers
      while (answer == 0)
        answer = menu(['The OME-TIFF version of ' printname ' already exists, overwrite it ?'],'Yes', 'Yes to All','No','No to All');
      end

    % Otherwise, we use the provided policy
    else
      answer = policy;
    end

    % Store the answer in case we need to output it
    policy = answer;

    % If it's a "to All" answer, act as it's a standard yes/no
    if (mod(answer, 2) == 0)
      answer = answer - 1;

    % Otherwise, reset the policy for the potential next file
    else
      policy = 0;
    end

    % Act accorindly to the policy
    switch answer

      % Delete the current files (did not dare overwriting it directly)
      case 1
        delete(newname);

      % Otherwise we can stop here
      case 3
        % Store the new name
        newfile = newname;
        newfile = relativepath(newfile);

        % Close the file handler
        r.close();

        return;
    end
  end

  % Initialize the file writer
  w.setId(newname);

  % Retrieve the compression mode from the options
  compression = opts.compression;

  % If none is provided, ask the user
  if(isempty(compression))

    % Get the available compressions
    compress = w.getCompressionTypes();
    [sel, ok] = listdlg('PromptString', 'Select a Compression',...
                  'SelectionMode', 'single',...
                  'ListString', char(compress));

    % If none is selected, use the first one
    if(ok == 0)
      sel = 1;
    end

    % Get the string
    compression = compress(sel);
  end

  % Apply the compression
  w.setCompression(java.lang.String(compression));

  % A nice status bar
  str = [' Converting ' printname ' from [' format '] to [' char(compression) ' OME-TIFF]'];
  hwait = waitbar(0, str,'Name','Bio-Formats Library');

  % Simply copy everything
  numSeries = r.getSeriesCount();
  for s=0:(numSeries-1)

    % Bio-formats requirement
    r.setSeries(s);
    
    % Loop over the frames
    numImages = r.getImageCount();
    for i=0:numImages-1

      % Simply read and write everything
      img = r.openBytes(i);
      w.saveBytes(i, img);
      
      % Show our progress
      waitbar(double(s*numSeries + (i+1))/(numSeries*numImages),hwait);
    end
  end

  % Close both handlers
  r.close();
  w.close();

  % Close the status bar
  close(hwait);
  
  % Store the new name
  newfile = newname;
  newfile = relativepath(newfile);

  return;
end
