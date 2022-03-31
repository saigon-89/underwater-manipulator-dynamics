# Underwater mobile-base manipulator dynamics calculator
## Assumptions
- mass centers of the links are equal to the centers of buoyancy
- friction forces aren't included in the calculations
- end-effector dynamics aren't considered
- damping forces aren't considered
- numeric functions for dynamics matrices generated with MATLAB symbolic -> numeric conversion
## D-H parameters
To describe the robot kinematics, the Denavit-Hartenberg representation is used

**D-H table of RRR-robot example:**
| Link | a | ⍺ | d | q |
|:-:|:-:|:-:|:-:|:-:|
| 1 | 0 | pi/2 | l(1) | q(1) |
| 2 | l(2) | 0 | 0 | q(2) |
| 3 | l(3) | 0 | 0 | q(3) |

## TO DO List
- [X] Symbollic check of method
- [ ] Hydrodynamic part calculation
    - [X] Added mass 
    - [ ] Damping effects
