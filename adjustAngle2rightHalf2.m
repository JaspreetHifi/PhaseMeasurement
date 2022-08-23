function deltaRH = adjustAngle2rightHalf2( delta )

deltaRH = delta;

index = find(delta < -pi/2);
deltaRH(index) = deltaRH(index)+pi;

index = find(delta > pi/2);
deltaRH(index) = deltaRH(index)-pi;

% if delta < -pi/2
%     deltaRH = delta + pi;
% elseif delta > pi/2
%     deltaRH = delta - pi;
% else
%     deltaRH = delta;
% end
