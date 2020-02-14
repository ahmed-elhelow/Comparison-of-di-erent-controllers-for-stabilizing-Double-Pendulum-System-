clear all;
syms u x_dot phi1 phi2;

m0 = 0.8; m1 = 0.5; m2 = 0.3;
l1 = 0.3; l2 = 0.2;
b = 0.1; g = 9.8; J = 0.006;

a0 = m0 +m1 + m2;
a1 = m1*l1 + m2*l1;
a2 = m1*l1^2 + m2*l1^2 + J;
a3 = m2*l2;
a4 = m2*l1*l2;
a5 = m2*l2^2 + J;

Z = [a0 a1 a3; 
     a1 a2 a4; 
     a3 a4 a5];
N = [u - b*x_dot; 
     a1*g*phi1; 
     a3*g*phi2];

% numerical_answer is a vector of x_d_dot, phi1_d_dot, phi2_d_dot
numerical_answer = Z\N;

A = zeros(6,6);
B = zeros(6,1);
C = eye(6);
D = zeros(6,1);

A(1,2) = 1; A(3,4) = 1 ; A(5,6) = 1;

v = [x_dot phi1 phi2 u];
r = [2 4 6];
c = [2 3 5];

for i = 1:4
    v_subs = subs(v,v(:),zeros(4,1));
    v_subs(i) = 1;
    col = double(subs(numerical_answer,v,v_subs));
    for j = 1:3
        if i<4
            A(r(j),c(i)) = col(j);
        else
            B(r(j)) = col(j);
        end
    end
end

P1 = [-5 -5.5 -6 -6.5 -7 -7.5];
K_pp_1 = place(A,B,P1);

P2 = [-3 -3.5 -5 -5.5 -1.5 -2];
K_pp_2 = place(A,B,P2);

Q = diag([10 10 500 100 500 100]);
R = 1;
K_lqr = lqr(A,B,Q,R);
