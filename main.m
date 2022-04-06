clear; clc;

%% GUI Input Form
answer = inputdlg('link number (n)','Enter link number', [1 30]);
n = str2num(answer{1}); % link number

prompt = {'\sigma [0 if revolute, 1 if prismatic]:',
    'link length [m]:',
    'link radius [m]:',
    'mass of each link [kg]:',  
    'buoyancy of each link [N]:', 
    'friction of each joint [const]:', 
    'R_0 rotation matrix [ex.: rotx(deg)]:'};
dlgtitle = 'Input parameters';
dims = [1 45];
definput = {'zeros(n,1)','[]','[]','[]','[]','zeros(n,1)','eye(3)'};
options.Interpreter = 'tex';
answer = inputdlg(prompt,dlgtitle,dims,definput,options);

% sigma = sym('sigma', [n 1], 'integer');
% l = sym('l', [n 1], 'real'); % link length [m]
% r = sym('r', [n 1], 'real'); % radius [m]
% m = sym('m', [n 1], 'real'); % mass of each link [kg]
% B = sym('B', [n 1], 'real'); % buoyancy [N] 

sigma = str2num(answer{1});
l = str2num(answer{2});
r = str2num(answer{3});
m = str2num(answer{4});
B = str2num(answer{5});
mu = str2num(answer{6});
Rot0 = str2num(answer{7});

q = sym('q', [n 1], 'real'); % generalized coordinates vector

%% D-H Parameter Matrix
answer = inputdlg('you can use q(i) and l(i) values  [ a | alpha | d | q ]', ...
    'Enter D-H table', [n 35], {'[]'});
eval(strcat('DH = ',answer{1}));   

%% Mass centers coordinates in link-based frame
c = cell(n,1);
prompt = cell(n,1);
for i = 1:n
    prompt{i} = strcat(num2str(i),'-st center');
end
answer = inputdlg(prompt,'Link-based frame mass centers', [1 55]);
for i = 1:n
    c{i} = str2num(answer{i}); 
end

dq = sym('dq', [n 1], 'real'); % joint velocities

%% Constants
g = 9.81; % gravitational acceleration constant
rho = 1000; % fluid density

%% Function handlers
% Rotation matrix from D-H table
A = @(a, alpha, d, q) ...
        ([cos(q) -sin(q)*cos(alpha)  sin(q)*sin(alpha) a*cos(q); ...
          sin(q)  cos(q)*cos(alpha) -cos(q)*sin(alpha) a*sin(q); ...
          0       sin(alpha)         cos(alpha)        d; ...
          0       0                  0                 1]);
% Added mass of cylindrical body
Add = @(m, r, l) ...
        (diag([0.1*m; ...
          pi*rho*r^2*l; ...
          pi*rho*r^2*l; ...
          0; ...
          (pi*rho*r^2*l^3)/12; ...
          (pi*rho*r^2*l^3)/12]));
% Rotation matrix
Rot = @(eta) ...
        [cos(eta(6))*cos(eta(5)) sin(eta(6))*cos(eta(5)) -sin(eta(5)); ...
        -sin(eta(6))*cos(eta(4)) + cos(eta(6))*sin(eta(5))*sin(eta(4)) ...
        cos(eta(6))*cos(eta(4)) + sin(eta(6))*sin(eta(5))*sin(eta(4)) ...
        sin(eta(4))*cos(eta(5));
        sin(eta(6))*sin(eta(4)) + cos(eta(6))*sin(eta(5))*cos(eta(4)) ...
        -cos(eta(6))*sin(eta(4)) + sin(eta(6))*sin(eta(5))*cos(eta(4)) ...
        cos(eta(4))*cos(eta(5))];
% Inertia tensors of cylindrical body
Iox = @(m, r, l) (m.*diag([0.5*r^2, (3*r^2 + l^2)/12, (3*r^2 + l^2)/12]));
Ioy = @(m, r, l) (m.*diag([(3*r^2 + l^2)/12, 0.5*r^2, (3*r^2 + l^2)/12]));
Ioz = @(m, r, l) (m.*diag([(3*r^2 + l^2)/12, (3*r^2 + l^2)/12, 0.5*r^2]));
    
%% Base generalized coordinates (6-DOF)
syms x y z phi theta psi
eta = [ x; y; z; phi; theta; psi ];

%% Homogeneous transformations solution
Tr = cell(n,1);
Tr0 = eye(4,'sym'); Tr0(1:3,1:3) = Rot0 * Rot(eta); Tr0(1:3,4) = eta(1:3);
T = Tr0;
for i = 1:n
    T = T * A(DH(i,1),DH(i,2),DH(i,3),DH(i,4));
    Tr{i} = simplify(T);
end

%% Mass centers for each link
r_c_m = cell(n,1);
for i = 1:n    
    temp = Tr{i}*[[1;0;0;0] [0;1;0;0] [0;0;1;0] [c{i};1]];
	r_c_m{i} = temp(1:3,4);  
end

%% Inertia tensor for each link relative to the inertial frame 
% computation depends on COM offset
I = cell(1,n);
for i = 1:n
    [~, index] = max(abs(c{i}));
    if (index == 1) 
        I{i} = Iox(m(i), r(i), l(i));
    elseif (index == 2)
        I{i} = Ioy(m(i), r(i), l(i));
    else
        I{i} = Ioz(m(i), r(i), l(i));
    end
end

%% Velocity Jacobians
Jv = cell(1,n);
Jw = cell(1,n);
for i = 1:n
    z = Tr0(1:3,3);
    o = Tr0(1:3,4);
    Jvt = sym(zeros(3,n));
    Jwt = sym(zeros(3,n));
    for j = 1:i
        if (sigma(i) == 0)
            Jvt(:,j) = cross(z, r_c_m{i}-o);
            Jwt(:,j) = z;
        else
            Jvt(:,j) = z;
            Jwt(:,j) = 0;
        end
        z = Tr{j}(1:3,3);
        o = Tr{j}(1:3,4);
    end
    Jv{i} = simplify(Jvt);
    Jw{i} = simplify(Jwt);
end

%% Same result calculation
% Jv = cell(1,n);
% Jw = cell(1,n);
% for i = 1:n
%     z = Tr0(1:3,3);
%     Jvt = sym(zeros(3,n));
%     Jwt = sym(zeros(3,n));
%     for j = 1:i
%         Jvt(:,j) = diff(r_c_m{i}, q(j));
%         if (sigma(i) == 0)  
%             Jwt(:,j) = z;
%         else
%             Jwt(:,j) = 0;
%         end
%         z = Ti{j}(1:3,3);
%     end
%     Jv{i} = Jvt;
%     Jw{i} = Jwt;
% end

%% Potential energy solution / Inertia matrix solution
PE = 0;
M_sym = 0;
for i = 1:n
    R = Tr{i}(1:3,1:3);
    A = Add(m(i), r(i), l(i));
    M_sym = M_sym + (Jv{i}'*(m(i)*eye(3)+A(1:3,1:3))*Jv{i} + ...
        Jw{i}'*R*(I{i}+A(4:6,4:6))*R'*Jw{i});
    %M_sym = M_sym + (m(i)*Jv{i}'*Jv{i} + Jw{i}'*R*I{i}*R'*Jw{i});
    %M_sym = M_sym + (Jv{i}'*A(1:3,1:3)*Jv{i} + Jw{i}'*R*A(4:6,4:6)*R'*Jw{i});
    PE = PE + (m(i)*g - B(i))*r_c_m{i}(3);
end

%% The Christoffel symbols
c = zeros(n,n,n,'sym');
for k = 1:n
    for i = 1:n
        for j = 1:n
            c(i,j,k) = 0.5 * (diff(M_sym(k,j),q(i)) + ...
                diff(M_sym(k,i),q(j)) - diff(M_sym(i,j),q(k)));
        end
    end
end

%% The Coriolis matrix
C_sym = zeros(n,n,'sym');
for k = 1:n
    for j = 1:n
        temp = 0;
        for i = 1:n
            temp = temp + c(i,j,k)*dq(i);
        end
        C_sym(k,j) = temp;
    end
end

%% Friction terms
D_sym = diag(mu);

%% The gravitation terms
g_sym = zeros(n,1,'sym');
for k = 1:n
    g_sym(k) = diff(PE,q(k));
end

%% Same result calculation
% g_sym = zeros(n,1,'sym');
% g_vect = [0; 0; -g];
% for i = 1:n
%     for k = 1:n
%         g_sym(i) = g_sym(i) - Jv{k}(:,i)' * m(k) * g_vect;
%     end
% end

%% Dynamics matrices
% M_sym = simplify(M_sym);
% C_sym = simplify(C_sym);
% g_sym = simplify(g_sym);

%% Generate numeric functions
matlabFunction(M_sym,'File','get_M','Vars',{[eta;q]});
matlabFunction(C_sym,'File','get_C','Vars',{[eta;q;dq]});
matlabFunction(g_sym,'File','get_g','Vars',{[eta;q]});

%% Motion equations
% ddq = sym('ddq', [n 1], 'real'); 
% tau = M_sym*ddq + C_sym*dq + D_sym*dq + g_sym
