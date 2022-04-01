# Underwater mobile-base manipulator dynamics calculator
1. [Assumptions](https://github.com/saigon-89/underwater-mobile-manipulator-dynamics/blob/main/README.md#assumptions)
2. [Description](https://github.com/saigon-89/underwater-mobile-manipulator-dynamics/blob/main/README.md#description)
    1. [Input GUI](https://github.com/saigon-89/underwater-mobile-manipulator-dynamics/blob/main/README.md#input-gui)
    2. [Constants](https://github.com/saigon-89/underwater-mobile-manipulator-dynamics/blob/main/README.md#constants)
    3. [Function handlers](https://github.com/saigon-89/underwater-mobile-manipulator-dynamics/blob/main/README.md#function-handlers)
    4. [Base generalized coordinates (6-DOF)](https://github.com/saigon-89/underwater-mobile-manipulator-dynamics#base-generalized-coordinates-6-dof)
    5. [Homogeneous transformations solution](https://github.com/saigon-89/underwater-mobile-manipulator-dynamics#homogeneous-transformations-solution)
    6. [Mass centers for each link](https://github.com/saigon-89/underwater-mobile-manipulator-dynamics#mass-centers-for-each-link)
    7. [Inertia tensors](https://github.com/saigon-89/underwater-mobile-manipulator-dynamics#inertia-tensors)
    8. [Jacobians](https://github.com/saigon-89/underwater-mobile-manipulator-dynamics#jacobians)
    9. [Potential energy and inertia matrix](https://github.com/saigon-89/underwater-mobile-manipulator-dynamics#potential-energy-and-inertia-matrix)
    10. [The Christoffel symbols](https://github.com/saigon-89/underwater-mobile-manipulator-dynamics#the-christoffel-symbols)
    11. [The Coriolis matrix](https://github.com/saigon-89/underwater-mobile-manipulator-dynamics#the-coriolis-matrix)
    12. [The gravitation terms](https://github.com/saigon-89/underwater-mobile-manipulator-dynamics#the-gravitation-terms)
    13. [Numeric functions generation](https://github.com/saigon-89/underwater-mobile-manipulator-dynamics#numeric-functions-generation)
3. [D-H parameters](https://github.com/saigon-89/underwater-mobile-manipulator-dynamics#d-h-parameters)
4. [Sources](https://github.com/saigon-89/underwater-mobile-manipulator-dynamics#sources)

## Assumptions
- each link approximated as cylindrical body
- mass centers of the links are equal to the centers of buoyancy
- friction forces aren't included in the calculations
- end-effector dynamics aren't considered
- damping forces aren't considered
- numeric functions for dynamics matrices generated with MATLAB `symbolic -> numeric` conversion

## Description
### Input GUI
...
### Constants
User can vary constants depending on environment
```matlab
g = 9.81; % gravitational acceleration constant
rho = 1000; % fluid density
```
### Function handlers
Function handlers are used to simplify readability of code

Rotation matrix from D-H table row [1]
```matlab
A = @(a, alpha, d, q) ...
        ([cos(q) -sin(q)*cos(alpha)  sin(q)*sin(alpha) a*cos(q); ...
          sin(q)  cos(q)*cos(alpha) -cos(q)*sin(alpha) a*sin(q); ...
          0       sin(alpha)         cos(alpha)        d; ...
          0       0                  0                 1]);
```

Added mass of cylindrical body [3]
```matlab
Add = @(m, r, l) ...
        (diag([0.1*m; ...
          pi*rho*r^2*l; ...
          pi*rho*r^2*l; ...
          0; ...
          (pi*rho*r^2*l^3)/12; ...
          (pi*rho*r^2*l^3)/12]));
```

[Rotation matrix](https://en.wikipedia.org/wiki/Rotation_matrix) from three Euler angles of base [1]
```matlab
Rot = @(eta) ...
        [cos(eta(6))*cos(eta(5)) sin(eta(6))*cos(eta(5)) -sin(eta(5)); ...
        -sin(eta(6))*cos(eta(4)) + cos(eta(6))*sin(eta(5))*sin(eta(4)) ...
        cos(eta(6))*cos(eta(4)) + sin(eta(6))*sin(eta(5))*sin(eta(4)) ...
        sin(eta(4))*cos(eta(5));
        sin(eta(6))*sin(eta(4)) + cos(eta(6))*sin(eta(5))*cos(eta(4)) ...
        -cos(eta(6))*sin(eta(4)) + sin(eta(6))*sin(eta(5))*cos(eta(4)) ...
        cos(eta(4))*cos(eta(5))];
```

[Inertia tensors](https://vk.com/away.php?to=https%3A%2F%2Fen.wikipedia.org%2Fwiki%2FList_of_moments_of_inertia) of cylindrical body for each position along axes 
```matlab
Iox = @(m, r, l) (m.*diag([0.5*r^2, (3*r^2 + l^2)/12, (3*r^2 + l^2)/12]));
Ioy = @(m, r, l) (m.*diag([(3*r^2 + l^2)/12, 0.5*r^2, (3*r^2 + l^2)/12]));
Ioz = @(m, r, l) (m.*diag([(3*r^2 + l^2)/12, (3*r^2 + l^2)/12, 0.5*r^2]));
```

### Base generalized coordinates (6-DOF)
Base position described with 6-DOF. Underwater vehicle is possible manipulator base
```matlab
syms x y z phi theta psi
eta = [ x; y; z; phi; theta; psi ];
```

### Homogeneous transformations solution
Matrix `Tr0` contains rotation and translation terms of base 
```matlab
Tr = cell(n,1);
Tr0 = eye(4,'sym'); Tr0(1:3,1:3) = Rot0 * Rot(eta); Tr0(1:3,4) = eta(1:3);
T = Tr0;
for i = 1:n
    T = T * A(DH(i,1),DH(i,2),DH(i,3),DH(i,4));
    Tr{i} = simplify(T);
end
```

### Mass centers for each link
Mass centers are variables of `q(i)` values, so we need to find general solution with rotation terms for each link COM (`r_c_m{i}`)
```matlab
r_c_m = cell(n,1);
for i = 1:n    
    temp = Tr{i}*[[1;0;0;0] [0;1;0;0] [0;0;1;0] [c{i};1]];
	r_c_m{i} = temp(1:3,4);  
end
```

### Inertia tensors
Inertia tensor for each link relative to the inertial frame 

computation depends on COM offset
```matlab
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
```

### Jacobians
Velocity Jacobians calculations can be found in [1]
```matlab
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
```

### Potential energy and inertia matrix
Potential energy calculated with hight along **OZ** axis
```matlab
PE = 0;
M_sym = 0;
for i = 1:n
    R = Tr{i}(1:3,1:3);
    A = Add(m(i), r(i), l(i));
    M_sym = M_sym + (Jv{i}'*(m(i)*eye(3)+A(1:3,1:3))*Jv{i} + ...
        Jw{i}'*R*(I{i}+A(4:6,4:6))*R'*Jw{i});
    PE = PE + (m(i)*g - B(i))*r_c_m{i}(3);
end
```

### The Christoffel symbols
[The Christoffel symbols](https://en.wikipedia.org/wiki/Christoffel_symbols) are used to calculate Coriolis terms

*Note that calculation requires only inertia matrix at input*
```matlab
c = zeros(n,n,n,'sym');
for k = 1:n
    for i = 1:n
        for j = 1:n
            c(i,j,k) = 0.5 * (diff(M_sym(k,j),q(i)) + ...
                diff(M_sym(k,i),q(j)) - diff(M_sym(i,j),q(k)));
        end
    end
end
```

### The Coriolis matrix
```matlab
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
```

### The gravitation terms
Potential energy gradient [1]
```matlab
g_sym = zeros(n,1,'sym');
for k = 1:n
    g_sym(k) = diff(PE,q(k));
end
```

### Numeric functions generation
Symbollic expressions replacing by numeric with [matlabFunction](https://www.mathworks.com/help/symbolic/matlabfunction.html)
```matlab
matlabFunction(M_sym,'File','get_M','Vars',{[eta;q]});
matlabFunction(C_sym,'File','get_C','Vars',{[eta;q;dq]});
matlabFunction(g_sym,'File','get_g','Vars',{[eta;q]});
```

## D-H parameters
To describe the robot kinematics, the Denavit-Hartenberg representation is used [2]

**D-H table of RRR-robot example:** (can be found in `examples/RRR_robot.m`)
| Link | a | ⍺ | d | q |
|:-:|:-:|:-:|:-:|:-:|
| 1 | 0 | pi/2 | l(1) | q(1) |
| 2 | l(2) | 0 | 0 | q(2) |
| 3 | l(3) | 0 | 0 | q(3) |

*Kinematic model can be found [here](https://www.wolframcloud.com/objects/demonstrations/DenavitHartenbergParametersForAThreeLinkRobot-source.nb)*

# Sources
1. Mark W. Spong, Seth Hutchinson, and Mathukumalli Vidyasagar. Robot modeling and control. John Wiley and Sons, Hoboken (N.J.), 2006.
2. Bruno Siciliano and Oussama Khatib. Springer Handbook of Rob. Springer Berlin Heidelberg, 2008.
3. Antonelli, G. Underwater Robots. 3rd. Springer, 2014.

# TO DO List
- [X] Symbollic check of method
- [ ] Hydrodynamic part calculation
    - [X] Added mass 
    - [ ] Damping effects
