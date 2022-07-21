function [phi, theta, psi, Rot_matrix] = davenport(s_v, m_v)

sn = s_v(:,1);
sb = s_v(:,2);
mn = m_v(:,1);
mb = m_v(:,2);

w1 = 1;
w2 = 2;

X = w1*sb*sn' + w2*mb*mn';
sig = trace(X);
Z = [X(2,3)-X(3,2); X(3,1)-X(1,3); X(1,2)-X(2,1)];
S = X + X';

K = [sig Z'; Z (S-sig*eye(3))];
[e_vec, e_val] = eig(K);
[max_val, ind] = max(max(e_val));
beta_opt = e_vec(:,ind);
b1 = beta_opt(1);
b2 = beta_opt(2);
b3 = beta_opt(3);
b4 = beta_opt(4);
Rot_matrix = [(b1^2+b2^2-b3^2-b4^2) 2*((b2*b3)+(b1*b4)) 2*((b2*b4)-(b1*b3));
              2*((b2*b3)-(b1*b4)) (b1^2-b2^2+b3^2-b4^2) 2*((b3*b4)+(b1*b2));
              2*((b2*b4)+(b1*b3)) 2*((b3*b4)-(b1*b2)) (b1^2-b2^2-b3^2+b4^2)];
          
theta = -asin(Rot_matrix(3,1));
phi = asin(Rot_matrix(3,2)/cos(theta));
psi = asin(Rot_matrix(2,1)/cos(theta));
end