function et = en_kin(etastart,myY,kx,kz,u,v,w)

Np = 64;
ystart = 10;

if kx == 0
  lx = 2*pi;
else
  lx = abs(2*pi/kx);
end

if etastart>0
  iy0 = find(myY>ystart,1,'last');
  iy1 = find(myY>etastart,1,'last');
  idy = iy1:1:iy0;
  ly = myY(iy1)-myY(iy0);
  y = myY(idy);
else
  iy0 = find(myY>ystart,1,'last');
  idy = 1:iy0;
  ly = myY(1)-myY(iy0);
  y = myY(idy);
end
ny = length(idy);

if kz == 0 
  lz = 2*pi;
else
  lz = abs(2*pi/kz);
end

x = zeros(Np,1); z = x;
x(:,1) = (0:Np-1)/Np*lx;
z(:,1) = (0:Np-1)/Np*lz;

phi = kx*x + kz*z.';
cp = cos(phi);
sp = sin(phi);

U = zeros(length(x),length(z),ny);
V = U;
W = U;

for iy = 1:ny
 U(:,:,iy) = real(u(idy(iy)))*cp - imag(u(idy(iy)))*sp;
 V(:,:,iy) = real(v(idy(iy)))*cp - imag(v(idy(iy)))*sp;
 W(:,:,iy) = real(w(idy(iy)))*cp - imag(w(idy(iy)))*sp;
end

inK = -trapz(y,U.^2,3) -trapz(y,V.^2,3) -trapz(y,W.^2,3);

et = sum(sum(inK));
dv = 1/(ly*Np^2);
et = 0.5*et*dv;

end

