function xnext = Rocketdyn_LinearVx(TS,m,L,J,g,x,u)
% x = [theta, w, h, v, x, vx]
% u = [Fe Fth Ftl]
% Small Angle Aproximation has been made
%Linearized model matrices
A=[1  TS   0   0   0   0;
   0   1   0   0   0   0;
   0   0   1   TS  0   0;
   0   0   0   1   0   0;
   0   0   0   0   1   TS;
   0   0   0   0   0   1 ];

B=[   0        0               0;
      0   TS*L/(4*J)    -TS*L/(2*J);
      0        0               0;
    TS/m       0               0;
      0        0               0;
      0      -TS/m           -TS/m];
  
C=[   0;
      0;
      0;
      -g*TS;
      0;
      0];
  
  xnext=A*x+B*u+C;
end