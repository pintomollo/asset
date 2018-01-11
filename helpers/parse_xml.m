function xml_tree = parse_xml(data)
% PARSEXML converts an XML file to a MATLAB structure.
%
%   XML = PARSE_XML(DATA) parses the content of the string DATA,
%   converting in into an XML structure.
%
%   XML = PARSE_XML(FNAME) converts the content of the XML file FNAME.
%
% Based on the code provided by Matlab, with a number
% of simplifactions in the resulting structure.
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 09.05.2015

  % Initialize some variables
  xml_tree = [];

  node = struct('Name', '', 'Attributes', [],    ...
                'Data', '', 'Children', []);

  leaf = struct('Name', '', 'Value', '');

  % Work only on string
  if (ischar(data) && ~isempty(data))

    % Prepare some more data
    xml_input = [];

    % Clean the input string from non-ASCII characters and excessive spaces
    data = data(double(data)<128);
    data = strtrim(data);
    data = regexprep(data, '>\s+<', '><');

    % Either we have a string or a filename
    if (data(1) == '<')
      xml_input = org.xml.sax.InputSource();
      xml_input.setCharacterStream(java.io.StringReader(data));
    elseif (exist(data, 'file') == 2)
      xml_input = data;
    end

    if (~isempty(xml_input))
      % Load the data
      xml_tree = xmlread(xml_input);

      % Recurse over child nodes. This could run into problems 
      % with very deeply nested trees.
      xml_tree = parseChildNodes(xml_tree);
    end
  end

  return;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Code adapted from MATLAB:                             %
  % http://www.mathworks.com/help/matlab/ref/xmlread.html %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % ----- Local function PARSECHILDNODES -----
  function [children,attrs] = parseChildNodes(theNode)
  % Recurse over node children.
    children = [];
    attrs = [];
    if theNode.hasChildNodes
       childNodes = theNode.getChildNodes;
       numChildNodes = childNodes.getLength;

       if (numChildNodes == 1)
         children = makeStructFromNode(childNodes.item(0));
       else
         children = repmat(node, 1, numChildNodes);
         attrs = repmat(leaf, 1, numChildNodes);

          cc = 1;
          ca = 1;
          for count = 1:numChildNodes
              theChild = childNodes.item(count-1);
              item = makeStructFromNode(theChild);

              if (ischar(item))
                attrs(ca).Name = item;
                ca = ca + 1;
              elseif (isfield(item, 'Value'))
                attrs(ca) = item;
                ca = ca + 1;
              elseif (numel(item) > 1)
                children = [item, children];
                cc = cc + length(item);
              elseif (numel(item) == 1)
                children(cc) = item;
                cc = cc + 1;
              end
          end

          children = children(1:cc-1);
          attrs = attrs(1:ca-1);
        end
    end

    return;
  end

  % ----- Local function MAKESTRUCTFROMNODE -----
  function nodeStruct = makeStructFromNode(theNode)
  % Create structure of node info.

    nodeStruct = node;
    nodeStruct.Name = char(theNode.getNodeName);
    nodeStruct.Attributes = parseAttributes(theNode);
    [nodeStruct.Children, attrs] = parseChildNodes(theNode);

    if (~isempty(attrs))
      if (isempty(nodeStruct))
        nodeStruct.Attributes = attrs;
      else
        nodeStruct.Attributes = [nodeStruct.Attributes, attrs];
      end
    end

    if any(strcmp(methods(theNode), 'getData'))
       nodeStruct.Data = char(theNode.getData); 
    else
       nodeStruct.Data = '';
    end

    switch lower(nodeStruct.Name)
      case '#text'
        nodeStruct = nodeStruct.Data;
      case {'children', 'value'}
        nodeStruct = nodeStruct.Children;
      case 'imagevalue'
        for indx=length(nodeStruct.Attributes):-1:1
          if (isempty(nodeStruct.Attributes(indx).Value))
            nodeStruct.Data = nodeStruct.Attributes(indx).Name;
            nodeStruct.Attributes(indx) = [];
          elseif (strncmpi(nodeStruct.Attributes(indx).Name, 'Name', 4))
            nodeStruct.Name = nodeStruct.Attributes(indx).Value;
            nodeStruct.Attributes(indx) = [];
          end
        end
      case 'originalmetadata'
        for indx=length(nodeStruct.Attributes):-1:1
          if (isempty(nodeStruct.Attributes(indx).Value))
            nodeStruct.Data = nodeStruct.Attributes(indx).Name;
            nodeStruct.Attributes(indx) = [];
          elseif (strncmpi(nodeStruct.Attributes(indx).Name, 'Key', 4))
            nodeStruct.Name = nodeStruct.Attributes(indx).Value;
            nodeStruct.Attributes(indx) = [];
          end
        end
      otherwise
        if (isempty(nodeStruct.Attributes) && isempty(nodeStruct.Data))
          if (ischar(nodeStruct.Children))
            tmp = leaf;
            tmp.Name = nodeStruct.Name;
            tmp.Value = nodeStruct.Children;
            nodeStruct = tmp;
          elseif (isempty(nodeStruct.Children))
            nodeStruct = nodeStruct.Name;
          elseif (~isfield(nodeStruct.Children, 'Children'))
            nodeStruct.Attributes = nodeStruct.Children;
            nodeStruct.Children = [];
          elseif (numel(nodeStruct.Children)==1 && strcmp(nodeStruct.Name, nodeStruct.Children.Name))
            nodeStruct = nodeStruct.Children;
          end
        elseif (ischar(nodeStruct.Children) && isempty(nodeStruct.Data))
            nodeStruct.Data = nodeStruct.Children;
            nodeStruct.Children = [];
        end
    end

    return;
  end

  % ----- Local function PARSEATTRIBUTES -----
  function attributes = parseAttributes(theNode)
  % Create attributes structure.

    attributes = [];
    if theNode.hasAttributes
       theAttributes = theNode.getAttributes;
       numAttributes = theAttributes.getLength;

        if (numAttributes == 1)
          name = char(theAttributes.item(0).getName);

          if (~strncmp(name, 'i:nil', 5))
            attributes = leaf;
            attributes.Name = name;
            attributes.Value = char(theAttributes.item(0).getValue);
          end
        else
         attributes = repmat(leaf, 1, numAttributes);

         for count = 1:numAttributes
            attrib = theAttributes.item(count-1);
            attributes(count).Name = char(attrib.getName);
            attributes(count).Value = char(attrib.getValue);
         end
       end
    end

    return;
  end
end
