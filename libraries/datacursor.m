## Copyright (C) 2017 Pantxo Diribarne
## 
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*- 
## @deftypefn {} {@var{h} =} datacursor ()
## @deftypefnx {} {@var{h} =} datacursor (@var{hax})
##
## @seealso{}
## @end deftypefn

function [htmp] = datacursor (hf = [])
  if (isempty (hf))
    hf = get (0, "currentfigure")
    if (isempty (hf))
      error ("datacursor: no figure")
    endif
  elseif (! ishandle (hf))
    error ("datacursor: expect HF to be a valid figure handle")
  endif
  
  htmp = install_all_listeners (hf);
endfunction

function htmp = install_all_listeners (hf)
  htmp = annotation (hf, "textarrow", [0 0], [eps eps], ...
                     "headlength", 6, "headwidth", 6, ...
                     "visible", "off");
  install_listeners (hf, htmp);
  cursors = getappdata (hf, "__datacursor__");
  setappdata (hf, "__datacursor__", [cursors htmp]);
endfunction

function install_listeners (hf, htmp)
  haxes = get (hf, "children");
  haxes(! isaxes (haxes)) = [];
  tag = get (haxes, "tag");
  idx = strcmp (tag, "legend") | strcmp (tag, "colorbar") | ...
        strcmp (tag, "scribeoverlay");
  haxes(idx) = [];
  for hax = haxes
    dellistener (hax, "currentpoint")
    
    lsn = {@update_nearest, htmp};
    addlistener (hax, "currentpoint", lsn)
    addlistener (htmp, "beingdeleted", {@cleanlistener, hax, lsn})
  endfor
endfunction

function cleanlistener (h, e, h2, lsn)
  if (ishandle (h2))
    dellistener (h2, "currentpoint", lsn);
  endif
endfunction

function update_nearest (hax, e, ht)
  hf = ancestor (hax, "figure");
  if (getappdata (hf, "__datacursor__")(end) != ht)
    return;
  endif
  hax = get (hf, "currentaxes");
  hch = get (hax, "children");
  hch(! strcmp (get (hch, "type"), "line")) = [];
  if (isempty (hch))
    error ("datacursor: no line objects in the axes object")
  endif
  xd = get (hch, "xdata");
  yd = get (hch, "ydata");
  if (iscell (xd))
    xd = [xd{:}];
    yd = [yd{:}];
  endif
  cp = get (hax, "currentpoint")(1,1:2);
  if (! isempty (xd))
    [~, idx] = min ((cp(1)-xd).^2 + (cp(2)-yd).^2);
    [x, y] = axes2normalized (hax, xd(idx), yd(idx));
    set (ht, "string", sprintf ("(%g, %g)", xd(idx), yd(idx)), ...
         "x", [x+0.02 x], "y", [y+0.02 y],
         "visible", "on", "userdata", [xd(idx) yd(idx)])
  endif
endfunction

function [x0, y0] = axes2normalized (hax, x, y)
  xl = xlim (hax);
  yl = ylim (hax);
  
  if (strcmp (get (hax, "xscale"), "log"))
    xl = log10 (xl);
    x0 = (log10(x) - xl(1)) / diff(xl);
  else
    x0 = (x - xl(1)) / diff(xl);
  endif
    
  if (strcmp (get (hax, "yscale"), "log"))
    yl = log10 (yl);
    y0 = (log10(y) - yl(1)) / diff(yl);
  else
    y0 = (y - yl(1)) / diff(yl);
  endif
  
  pos = get (hax, "position");
  x0 = x0*pos(3) + pos(1);
  y0 = y0*pos(4) + pos(2);
endfunction
