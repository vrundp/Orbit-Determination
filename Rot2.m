function R = Rot2(ang)
%+========================================================================+
% Rot2: Rotation matrix for a rotation about the second axis by some angle.
%
% INPUTS:
%
% ang [rad] rotation angle
%
% OUTPUTS:
%
% R [-] 3x3 rotation matrix 
% 
%+------------------------------------------------------------------------+
% References: Schaub, H & Junkins, J. Analytical Mechanics of Space
%             Systems. Fourth Edition. Pg 93. 
%
% Author: Vrund Patel
%+========================================================================+

    R = [cos(ang), 0.0, -sin(ang); ...
              0.0, 1.0,       0.0; ...
         sin(ang), 0.0,  cos(ang)];

end