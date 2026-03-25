clear; syms L1 L2 L3 L4 L5 L6 L7 L8 L9 x y z u
eqn1= L1*x+L2*y==0;
eqn2= L3*x+L4*y+L5*z==0;
eqn3= L6*x+L7*y+L8*z==L9*u;
sol=solve(eqn1,eqn2,eqn3,x,y,z); G=sol.y/u
pause;


% G = -(L1*L5*L9)/(L1*L4*L8 - L1*L5*L7 - L2*L3*L8 + L2*L5*L6)


syms sig b xbar ybar zbar s
G=subs(sol.y/u,{L1,L2,L3,L4,L5,L6,L7,L8,L9},{s+sig,-sig,zbar,s+1,xbar,-ybar,-xbar,s+b,-b})
[numG,denG] = numden(G);      % this extracts out the num and den of G
numG=coeffs(numG,s);          % this extracts the powers of s in the num and den
denG=coeffs(denG,s);
numG=simplify(numG/denG(end)); % this makes the den monic
denG=simplify(denG/denG(end));

numG=numG(end:-1:1)   % this reverses the order of the vector of coefficients.
denG=denG(end:-1:1)

% Poles
% (c, c, -1)
polycoeff1 = [1 6 52 376];
poles1 = RR_roots(polycoeff1)
% (-c,-c,-1)
polycoeff2 = [1 6 52 0];
poles2 = RR_roots(polycoeff2)
% (0, 0, -u)
polycoeff3 = [1 6 -183 -188];
poels3 = RR_roots(polycoeff3)

% roots
% (c, c, -1)
numcoeff1 = [6.8556 27.4226];
zeros1 = RR_roots(numcoeff1)
% (-c,-c,-1)
numcoeff2 = [-6.8556 -27.4226];
zeros2 = RR_roots(numcoeff2)
% (0, 0, -u)
numcoeff3 = [0 0];
zeros3 = RR_roots(numcoeff3)
