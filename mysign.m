function out = mysign(a,b)

if b >= 0
  out = abs(a);
else
  out = -abs(a);
end
end