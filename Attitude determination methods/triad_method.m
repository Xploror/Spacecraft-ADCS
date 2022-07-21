function [phi, theta, psi, Rot_mat] = triad_method(s_v, m_v)
% Format : s_v ---> [sun_vector inertial, sun_vector body]
%          m_v ---> [magnetic north inertial, magnetic north body]

sn = s_v(:,1);
sb = s_v(:,2);
mn = m_v(:,1);
mb = m_v(:,2);

% Initializing triad T frame
t1_b = sb;
t1_n = sn;

t2_b = cross(sb, mb)/norm(cross(sb, mb),2);
t2_n = cross(sn, mn)/norm(cross(sn, mn),2);

t3_b = cross(t1_b, t2_b);
t3_n = cross(t1_n, t2_n);

BT = [t1_b, t2_b, t3_b];
NT = [t1_n, t2_n, t3_n];

Rot_mat = BT*NT';

theta = -asin(Rot_mat(3,1));
phi = asin(Rot_mat(3,2)/cos(theta));
psi = asin(Rot_mat(2,1)/cos(theta));

end