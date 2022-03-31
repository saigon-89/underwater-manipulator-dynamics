n = 2; % link number
sigma = [0; 0]; % 0 if revolute, 1 if prismatic
l = [1; 1]; % link length [m]
r = [0.1; 0.1]; % radius [m]
m = [1; 1]; % mass of each link [kg]
B = [0; 0]; % buoyancy [N] 
Rot0 = rotx(90); % rotation matrix

q = sym('q', [n 1], 'real'); % generalized coordinates vector

%%%%%% a  alpha d   q %%%%%%
DH = [l(1)  0   0  q(1); ...
      l(2)  0   0  q(2)];   
      
c{1} = [-l(1)/2; 0; 0];
c{2} = [-l(2)/2; 0; 0];
