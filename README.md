# Underwater mobile-base manipulator dynamics calculator
## Assumptions
- mass centers of the links are equal to the centers of buoyancy
- friction forces aren't included in the calculations
- end-effector dynamics aren't considered
- damping forces aren't considered
- numeric functions for dynamics matrices generated with MATLAB symbolic -> numeric conversion

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

Rotation matrix from D-H table [1]
```matlab
A = @(a, alpha, d, q) ...
        ([cos(q) -sin(q)*cos(alpha)  sin(q)*sin(alpha) a*cos(q); ...
          sin(q)  cos(q)*cos(alpha) -cos(q)*sin(alpha) a*sin(q); ...
          0       sin(alpha)         cos(alpha)        d; ...
          0       0                  0                 1]);
```

Added mass of cylindrical body [Antonelli/Fossen]
```matlab
Add = @(m, r, l) ...
        (diag([0.1*m; ...
          pi*rho*r^2*l; ...
          pi*rho*r^2*l; ...
          0; ...
          (pi*rho*r^2*l^3)/12; ...
          (pi*rho*r^2*l^3)/12]));
```

Rotation matrix [1]
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

Inertia tensors of cylindrical body [Wiki]
```matlab
Iox = @(m, r, l) (m.*diag([0.5*r^2, (3*r^2 + l^2)/12, (3*r^2 + l^2)/12]));
Ioy = @(m, r, l) (m.*diag([(3*r^2 + l^2)/12, 0.5*r^2, (3*r^2 + l^2)/12]));
Ioz = @(m, r, l) (m.*diag([(3*r^2 + l^2)/12, (3*r^2 + l^2)/12, 0.5*r^2]));
```

## D-H parameters
To describe the robot kinematics, the Denavit-Hartenberg representation is used

**D-H table of RRR-robot example:** *(can be found in `examples/RRR_robot.m`)*
| Link | a | ⍺ | d | q |
|:-:|:-:|:-:|:-:|:-:|
| 1 | 0 | pi/2 | l(1) | q(1) |
| 2 | l(2) | 0 | 0 | q(2) |
| 3 | l(3) | 0 | 0 | q(3) |

*Kinematic model can be found [here](https://www.wolframcloud.com/objects/demonstrations/DenavitHartenbergParametersForAThreeLinkRobot-source.nb)*

## TO DO List
- [X] Symbollic check of method
- [ ] Hydrodynamic part calculation
    - [X] Added mass 
    - [ ] Damping effects
