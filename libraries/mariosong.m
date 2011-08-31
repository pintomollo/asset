function mariosong(working)
% mariosong.m  : mario brother's theme song..   
% programming  : James Humes
% transcription: Stewart Bozarth
% http://www.mathworks.com/matlabcentral/fileexchange/8442

  if (nargin == 0)
    working = true;
  end

  t = 0.17;

  if (working)
    %Mario's intro
    %treble
    keyst = [ 56 56 0 56 0 52 56 0 59 0 0 47 0 0 ];
    tdur = [t t t t t t t t t t 2*t t t 2*t ];
    %bass
    keysb = [ 30 30 0 30 0 30 30 0 47 0 0 35 0 0];
    bdur = [ t t t t t t t t t t 2*t t t 2*t];
    %alto
    keysa= [ 46 46 0 46 0 46 46 0 51 0 0 47 0 0 ];
    adur = [ t t t t t t t t t t 2*t t t 2*t];
  else
    %Mario's lose life
    %treble
    keyst = [52 53 54 0 51 57 0 57 57 56 54 52 0 35 0 28 0];
    tdur = [(1/4)*t (1/4)*t (1/2)*t 3*t t t t t (2/3)*t (2/3)*t (2/3)*t t t t t t t 3*t];
    %bass
    keysb = [0 0 0 0 35 0 0 35 35 37 39 40 0 35 0 28 0];
    bdur = [(1/4)*t (1/4)*t (1/2)*t 3*t t t t t (2/3)*t (2/3)*t (2/3)*t t t t t t t 3*t];
    %alto
    keysa = [0 0 0 0 47 54 0 54 52 51 47 44 0 44 40 0];
    adur = [(1/4)*t (1/4)*t (1/2)*t 3*t t t t t (2/3)*t (2/3)*t (2/3)*t t t t t t t 3*t];
  end

  fs = 11025;
  xt = zeros(1, ceil(sum(tdur)*fs)+1);
  xb = zeros(1, ceil(sum(bdur)*fs)+1);
  xa = zeros(1, ceil(sum(adur)*fs)+1);

  n1=1;
  for kk = 1:length(keyst)
    keynum=keyst(kk);
    tone=note(keyst(kk), tdur(kk));
    n2=n1 + length(tone)-1;
    xt(n1:n2) = xt(n1:n2) + tone;
    n1 = n2;
  end
  
  n1 = 1;
  for kk = 1:length(keysb)
    keynum=keysb(kk);
    tone=note(keysb(kk), bdur(kk));
    n2=n1 +length(tone)-1;
    xb(n1:n2) = xb(n1:n2) + tone;
    n1=n2;
  end

  n1=1;
  for kk = 1:length(keysa)
    keynum=keysa(kk);
    tone=note(keysa(kk), adur(kk));
    n2=n1 +length(tone)-1;
    xa(n1:n2) = xa(n1:n2) + tone;
    n1=n2;
  end

  %the "mixing board"

  xx=xa+xb+xt;
  soundsc(xx, fs)

  return;
end

function tone=note(keynum, dur)

  fs=11025;
  tt = 0:(1/fs):dur;

  %generates rests
  if keynum == 0 
    tone([1:length(tt)]) = 0;
    return;
  end

  %adding these octaves rounded out the sound.  
  freq=440*2^((keynum-49)/12);
  freq3=freq*3;
  freq5=freq*5;
  freq9=freq*9;
  freq7=freq*7;
  tone1 = .75*sin(2*pi*freq*tt);
  tone3 = .65*sin(2*pi*freq3*tt);
  tone5=.5*sin(2*pi*freq5*tt);
  tone9 = .222*sin(2*pi*freq9*tt);
  tone7 = .12*sin(2*pi*freq7*tt);

  tone=tone1+tone3+tone5+tone7+tone9;

  return;
end
