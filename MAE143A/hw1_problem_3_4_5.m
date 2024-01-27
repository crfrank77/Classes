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

% 4
% Poles
% (c, c, -1)
polycoeff1 = [1 6 52 376];
poles1 = RR_roots(polycoeff1)
% (-c,-c,-1)
polycoeff2 = [1 6 52 376];
poles2 = RR_roots(polycoeff2)
% (0, 0, -u)
polycoeff3 = [1 6 -183 -188];
poles3 = RR_roots(polycoeff3)

% Roots
% (c, c, -1)
numcoeff1 = [6.8556 27.4226];
zeros1 = RR_roots(numcoeff1)
% (-c,-c,-1)
numcoeff2 = [-6.8556 -27.4226];
zeros2 = RR_roots(numcoeff2)
% (0, 0, -u)
numcoeff3 = [0 0];
zeros3 = RR_roots(numcoeff3)


% 5
% 1-> (c,c,-1) 2-> (-c,-c,-1) 3-> (0,0,-u)
A(1) = (numcoeff1(1)*(poles1(1))+numcoeff1(2)) / ( (poles1(1)) * (poles1(1)-poles1(2)) * (poles1(1)-poles1(3)) );
A(2) = (numcoeff2(1)*(poles2(1))+numcoeff2(2)) / ( (poles2(1)) * (poles2(1)-poles2(2)) * (poles2(1)-poles2(3)) );
A(3) = (numcoeff3(1)*(poles3(1))+numcoeff3(2)) / ( (poles3(1)) * (poles3(1)-poles3(2)) * (poles3(1)-poles3(3)) )

B(1) = (numcoeff1(1)*(poles1(2))+numcoeff1(2)) / ( (poles1(2)) * (poles1(2)-poles1(1)) * (poles1(2)-poles1(3)) );
B(2) = (numcoeff2(1)*(poles2(2))+numcoeff2(2)) / ( (poles2(2)) * (poles2(2)-poles2(1)) * (poles2(2)-poles2(3)) );
B(3) = (numcoeff3(1)*(poles3(2))+numcoeff3(2)) / ( (poles3(2)) * (poles3(2)-poles3(1)) * (poles3(2)-poles3(3)) )

C(1) = (numcoeff1(1)*(poles1(3))+numcoeff1(2)) / ( (poles1(3)) * (poles1(3)-poles1(1)) * (poles1(3)-poles1(2)) );
C(2) = (numcoeff2(1)*(poles2(3))+numcoeff2(2)) / ( (poles2(3)) * (poles2(3)-poles2(1)) * (poles2(3)-poles2(2)) );
C(3) = (numcoeff3(1)*(poles3(3))+numcoeff3(2)) / ( (poles3(3)) * (poles3(3)-poles3(1)) * (poles3(3)-poles3(2)) )
 
D(1) = (numcoeff1(1)*(0)+numcoeff1(2)) / ( (-poles1(1)) * (-poles1(2)) * (-poles1(3)) );
D(2) = (numcoeff2(1)*(0)+numcoeff2(2)) / ( (-poles2(1)) * (-poles2(2)) * (-poles2(3)) );
D(3) = (numcoeff3(1)*(0)+numcoeff3(2)) / ( (-poles3(1)) * (-poles3(2)) * (-poles3(3)) )

ydot1t = @(t) A(1)*exp(poles1(1)*t) + D(1) + exp(real(poles1(2))*t)*((B(1)+C(1))*cos(imag(poles1(2))*t)+(i*(B(1)-C(1)))*sin(imag(poles1(3))*t)); 
ydot2t = @(t) A(2)*exp(poles2(1)*t) + D(2) + exp(real(poles2(2))*t)*((B(2)+C(2))*cos(imag(poles2(2))*t)+(i*(B(2)-C(2)))*sin(imag(poles2(3))*t)); 
ydot3t = 0;
% RR_tff(double(NumY(1)), double(denY))


counter = 1;
for k = 1:.1:10
    plotter(counter) = ydot1t(k);
    plotter2(counter) = ydot2t(k);
    counter = counter +1;
end
xval = 1:0.1:10;
figure(1);
subplot(2,1,1)
plot(xval,plotter);
subplot(2,1,2);
plot(xval,plotter2);


%plot([1,10],plotter);


