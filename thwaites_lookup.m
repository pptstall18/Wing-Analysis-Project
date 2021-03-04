function H = thwaites_lookup(m)
%
%  function H = thwaites_lookup(m)
%
%  Returns interpolated value of shape factor 
%  for the specified value of Thwaites' 'm' parameter.  Based on tabulated
%  values in Young, p90.
%

if m <= -0.25     % outside range of table, return values at -0.25

%  l = 0.5;
  H = 2;

elseif m <= 0.09    % interpolate from tabulated values

  mtab = [-.25 -.2 -.14 -.12 -.1 -.08 -.064 -.048 -.032 -.016 .0 .016 .032 ...
             .04 .048 .056 .06 .064 .068 .072 .076 .080 .084 .086 .088 .09];

%  ltab = [.5 .463 .404 .382 .359 .333 .313 .291 .268 .244 .22 .195 .168 ... 
%            .153 .138 .122 .113 .104 .095 .085 .072 .056 .038 .027 .015 .0];

  Htab = [2 2.07 2.18 2.23 2.28 2.34 2.39 2.44 2.49 2.55 2.61 2.67 2.75 ...
           2.81 2.87 2.94 2.99 3.04 3.09 3.15 3.22 3.3 3.39 3.44 3.49 3.55];

%  l = interp1(mtab,ltab,m);
  H = interp1(mtab,Htab,m);

else      %  m > 0.09; b-l separated, return m = 0.09 values

%  l = 0;
  H = 3.55;

end

