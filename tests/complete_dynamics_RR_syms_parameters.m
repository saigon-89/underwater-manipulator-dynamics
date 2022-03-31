clear; clc;

n = 2; % число звеньев

l = sym('l', [n 1], 'real'); % длины звеньев
r = sym('r', [n 1], 'real'); % радиусы звеньев
m = sym('m', [n 1], 'real'); % mass of each link [kg]
R0 = rotx(0); % начальная матрица поворота

q = sym('q', [n 1], 'real'); % обобщенные координаты (углы соединений)
g = sym('g', 'real');

% Координаты центров масс каждого звена в его собственной системе отсчета
c{1} = [-l(1)/2; 0; 0];
c{2} = [-l(2)/2; 0; 0];

% cell array of your homogeneous transformations; each Ti{i} is a 4x4 symbolic transform matrix
% NOTE: for symbolic arrays: q(1) = q1, q(2) = q2, etc.
Ti = cell(n,1);

A = @(a, alpha, d, theta)...
        ([cos(theta) -sin(theta)*cos(alpha) sin(theta)*sin(alpha) a*cos(theta);...
        sin(theta) cos(theta)*cos(alpha) -cos(theta)*sin(alpha) a*sin(theta);...
        0 sin(alpha) cos(alpha)  d;...
        0 0 0 1]);

%      a  alpha d   q
DH = [l(1)  0   0  q(1); ...
      l(2)  0   0  q(2)];     % DH Parameter Matrix

% Homogeneous transformations solution
% https://automaticaddison.com/how-to-find-denavit-hartenberg-parameter-tables/
% DH = [0    pi/2  l(1)  q(1); ...
%      l(2)   0     0    q(2)];     % DH Parameter Matrix

T0 = eye(4); T0(1:3,1:3) = R0;
T = T0;
for i = 1:n
    temp = A(DH(i,1),DH(i,2),DH(i,3),DH(i,4));
    T = T*temp;
    Ti{i} = T;
end

%% Центры масс
com = cell(n,1);
for i=1:n    
P = Ti{i}*[[1;0;0;0] [0;1;0;0] [0;0;1;0] [c{i};1]];
    x = P(1,4);
    y = P(2,4);
    z = P(3,4);
    com{i} = [x; y; z];  
end

%% Potential energy
% Potential energy solution
P = eye(4);
PE = 0;

%% Inertia matrix and kinetic energy
qd = sym('qd', [n 1], 'real'); % "q dot" - the first derivative of the q's in time (joint velocities)

% inertia tensor for each link relative to the inertial frame stored in an nx1 cell array
I = cell(1,n);
syms I1_xx I1_yy I1_zz I2_xx I2_yy I2_zz
I{1} = diag([I1_xx I1_yy I1_zz]);
I{2} = diag([I2_xx I2_yy I2_zz]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Velocity Jacobians
Jv = cell(1,n);
Jw = cell(1,n);
for i = 1:n
    z = T0(1:3,3);
    o = T0(1:3,4);
    Jvt = sym(zeros(3,n));
    Jwt = sym(zeros(3,n));
    for j = 1:i
        Jvt(:,j) = cross(z, com{i}-o);
        Jwt(:,j) = z;
        z = Ti{j}(1:3,3);
        o = Ti{j}(1:3,4);
    end
    Jv{i} = Jvt;
    Jw{i} = Jwt;
end

%% ВТОРОЙ СПОСОБ (РЕЗУЛЬТАТ ТОТ ЖЕ!) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Jv = cell(1,n);
% Jw = cell(1,n);
% for i = 1:n
%     z = T0(1:3,3);
%     Jvt = sym(zeros(3,n));
%     Jwt = sym(zeros(3,n));
%     for j = 1:i
%         Jvt(:,j) = diff(com{i}, q(j));
%         Jwt(:,j) = z;
%         z = Ti{j}(1:3,3);
%     end
%     Jv{i} = Jvt;
%     Jw{i} = Jwt;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% M_sym = Inertia matrix solution & PE = Poterntial Energy
M_sym = 0;
for i = 1:n
    R = Ti{i}(1:3,1:3);
    M_sym = M_sym + (m(i)*Jv{i}'*Jv{i} + Jw{i}'*R*I{i}*R'*Jw{i});
    PE = PE + m(i)*g*com{i}(2);
end

%% Equations of motion
qdd = sym('qdd', [n 1], 'real'); % "q double dot" - the second derivative of the q's in time (joint accelerations)

% The Christoffel symbols
c = zeros(n,n,n,'sym');
for k = 1:n
    for i = 1:n
        for j = 1:n
            c(i,j,k) = 0.5 * (diff(M_sym(k,j),q(i)) + diff(M_sym(k,i),q(j)) - diff(M_sym(i,j),q(k)));
        end
    end
end

% The coriolis matrix
C_sym = zeros(n,n,'sym');
for k = 1:n
    for j = 1:n
        temp = 0;
        for i = 1:n
            temp = temp + c(i,j,k)*qd(i);
        end
        C_sym(k,j) = temp; %% ТУТ поменял индексы
    end
end

% The gravitation terms
g_sym = zeros(n,1,'sym');
for k = 1:n
    g_sym(k) = diff(PE,q(k));
end

%% ВТОРОЙ СПОСОБ (РЕЗУЛЬТАТ ТОТ ЖЕ!) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% g_sym = zeros(n,1,'sym');
% g_accel = [0; -9.81; 0];
% for i = 1:n
%     for k = 1:n
%         g_sym(i) = g_sym(i) - Jv{k}(:,i)' * m(k) * g_accel;
%     end
% end

M_sym = simplify(M_sym)
C_sym = simplify(C_sym)
g_sym = simplify(g_sym)


%% TO CHECK RESULTS [презентация - 57 слайд]: http://oramosp.epizy.com/teaching/18/robotics/lectures/Topic11_Dynamics_I.pdf?i=3
m11 = m(1)*(l(1)/2)^2 + m(2)*( l(1)^2 + 2*l(1)*l(2)*0.5*cos(q(2)) + (l(2)/2)^2 ) + I1_zz + I2_zz;
m12 = m(2)*( l(1)*(l(2)/2)*cos(q(2)) + (l(2)/2)^2 ) + I2_zz;
m21 = m12;
m22 = m(2)*(l(2)/2)^2 + I2_zz;

M_ch = [ m11, m12; 
         m21, m22 ];

h = -m(2)*l(1)*l(2)*0.5*sin(q(2));

C_ch = [  h*qd(2),  h*(qd(1)+qd(2));  
         -h*qd(1),  0  ];

g_ch = [ (m(1)*l(1)*0.5 + m(2)*l(1))*cos(q(1)) + m(2)*l(2)*0.5*cos(q(1)+q(2)); m(2)*l(2)*0.5*cos(q(1)+q(2)) ].*g;

simplify(M_ch - M_sym)
simplify(C_ch - C_sym)
simplify(g_ch - g_sym)