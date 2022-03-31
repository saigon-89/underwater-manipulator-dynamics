clear; clc;

n = 2; % число звеньев
l = [1; 1]; % длины звеньев
r = [0.1; 0.1]; % радиусы звеньев
m = [1; 1]; % mass of each link [kg]
R0 = rotx(90); % начальная матрица поворота

q = sym('q', [n 1], 'real'); % обобщенные координаты (углы соединений)

% Координаты центров масс каждого звена в его собственной системе отсчета
c{1} = [-l(1)/2; 0; 0];
c{2} = [-l(2)/2; 0; 0];

%      a  alpha d   q
DH = [l(1)  0   0  q(1); ...
      l(2)  0   0  q(2)];     % DH Parameter Matrix

% Homogeneous transformations solution
% https://automaticaddison.com/how-to-find-denavit-hartenberg-parameter-tables/
% DH = [0    pi/2  l(1)  q(1); ...
%      l(2)   0     0    q(2)];     % DH Parameter Matrix

% inertia tensor for each link relative to the inertial frame stored in an nx1 cell array
I = cell(1,n);
I{1} = m(1).*diag([0.5*r(1)^2, (3*r(1)^2 + l(1)^2)/12, (3*r(1)^2 + l(1)^2)/12]);
I{2} = m(2).*diag([0.5*r(2)^2, (3*r(2)^2 + l(2)^2)/12, (3*r(2)^2 + l(2)^2)/12]);

g_accel = 9.81;


% cell array of your homogeneous transformations; each Ti{i} is a 4x4 symbolic transform matrix
% NOTE: for symbolic arrays: q(1) = q1, q(2) = q2, etc.
Ti = cell(n,1);

A = @(a, alpha, d, theta)...
        ([cos(theta) -sin(theta)*cos(alpha) sin(theta)*sin(alpha) a*cos(theta);...
        sin(theta) cos(theta)*cos(alpha) -cos(theta)*sin(alpha) a*sin(theta);...
        0 sin(alpha) cos(alpha)  d;...
        0 0 0 1]);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
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

%% ВТОРОЙ СПОСОБ (РЕЗУЛЬТАТ ТОТ ЖЕ!) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    PE = PE + m(i)*g_accel*com{i}(3);
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

M_sym = simplify(M_sym);
C_sym = simplify(C_sym);
g_sym = simplify(g_sym);

matlabFunction(M_sym,'File','get_M','Vars',{q});
matlabFunction(C_sym,'File','get_C','Vars',{[q;qd]});
matlabFunction(g_sym,'File','get_g','Vars',{q});

%% TO DO: added mass + added Coriolis mass

%% ode45 calculation: свободное падение при отклонении!
q0 = deg2rad([-95; 0]); 
dq0 = zeros(n,1);
t_end = 60; dt = 0.01;
[t,Y] = ode45(@(t,y)odefcn(t,y), 0:dt:t_end, [q0; dq0]);

%% Plot graphs
figure
q_vect = Y(:,1:n);
title('Значения q_i(t)'), hold on, grid on
plot(t, rad2deg(q_vect)), xlabel('t, сек'), ylabel('Положение, deg'), xlim([0 t_end])
labels = [];
for i = 1:n
    labels = [labels; strcat('q_', num2str(i), '(t)')];
end
legend(cellstr(labels))

%% Построение манипулятора 
figure
for i = [1, length(q_vect)]
    x = 0; y = 0; z = 0;
    for j = 1:n
        tmp = subs(Ti{j}(1:3,4), q, q_vect(i,:)');
        x = [x; tmp(1)]; y = [y; tmp(2)]; z = [z; tmp(3)];
    end
    plot3(x, y, z, 'r*'), hold on
    plot3(x, y, z, 'r-'), hold on
end
xlabel('X'), ylabel('Y'), zlabel('Z')
grid on
axis equal

function dy = odefcn(t,y)
    n = numel(y)/2;
    dy = zeros(2*n,1);
    dy(1:n) = y(n+1:end);
    dy(n+1:end) = get_M(y(1:n)) \ ...
        ( -get_C(y)*y(n+1:end) - get_g(y(1:n)) );
end