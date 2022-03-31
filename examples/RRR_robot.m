n = 3; % link number
sigma = [0; 0; 0]; % 0 if revolute, 1 if prismatic
l = [0.75; 1; 1]; % link length [m]
r = [0.1; 0.1; 0.1]; % radius [m]
m = [0.75; 1; 1]; % mass of each link [kg]
B = [0; 0; 0]; % buoyancy [N] 
orient = ['z'; 'x'; 'x'];
Rot0 = eye(3); % rotation matrix

q = sym('q', [n 1], 'real'); % generalized coordinates vector

%%%%%% a    alpha d    q %%%%%%
DH = [ 0    pi/2  l(1) q(1); ...
       l(2) 0     0    q(2); ...
       l(3) 0     0    q(3)];   
       
c{1} = [0; -l(1)/2; 0];
c{2} = [-l(2)/2; 0; 0];
c{3} = [-l(3)/2; 0; 0];
