function R = Rot3(ang)
%+========================================================================+
% Rot3: Rotation matrix for a rotation about the third axis by some angle.
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

    R = [ cos(ang), sin(ang), 0.0; ...
         -sin(ang), cos(ang), 0.0; ...
               0.0,      0.0, 1.0];

end