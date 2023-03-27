function display(obj)

  if (ishandle(obj.hPolygon))
    data = get(obj.hPolygon, 'userdata');

    printf ("handle =\n");
    display(obj.hPolygon)
    printf ("position =\n");
    display(data.position)
    printf ("is_closed =\n");
    display(data.is_closed)
    printf ("color =\n");
    display(data.color)

  else
    display([]);
  end

  return;
end
