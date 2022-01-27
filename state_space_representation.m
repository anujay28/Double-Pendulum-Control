syms theta_1 theta_2 theta_1dot theta_2dot theta_1ddot theta_2ddot
syms l1 l2 m1 m2 g r1 r2 i1 i2

%Co-ordinates and Velocity of m1
xm1 = r1*sin(theta_1);
ym1 = r1*cos(theta_1);
xm1_dot = diff(xm1,theta_1)*theta_1dot;
ym1_dot = diff(ym1 ,theta_1)*theta_1dot;
v1 = sqrt(xm1_dot^2 + ym1_dot^2);

%Co-ordinates and Velocity of m2
xm2 = l1*sin(theta_1) + r2*sin(theta_1 + theta_2) ;
ym2 = r2*cos(theta_1 + theta_2) + l1*cos(theta_1);
J1 = jacobian([xm2,ym2],[theta_1,theta_2]);
A = [theta_1dot;theta_2dot ];
B = J1*A;
xm2_dot = B(1);
ym2_dot = B(2);
v2 = sqrt(xm2_dot^2 + ym2_dot^2);

%Defining Kinetic Energy and Potential Energy
KE = 1/2*m1*(v1^2) + 1/2*m2*(v2^2) + 1/2*i1*(theta_1dot^2) + 1/2*i2*((theta_1dot + theta_2dot)^2);
PE = m1*g*r1*cos(theta_1) + m2*g*( l1*cos(theta_1) + r2*cos(theta_1 + theta_2) );

%Lagrangian Function 
L  = KE - PE ;

%Use of Jacobian to find partial derivatives 
J2 = jacobian(L , [theta_1 ,theta_1dot,theta_2,theta_2dot]);
J3 = jacobian( [J2(1,2), J2(1,4)], [theta_1dot , theta_2dot, theta_1 ,theta_2 ] );
C =[theta_1ddot;theta_2ddot ; theta_1dot ; theta_2dot];
D = J3*C;


%part a - Equations of Motion (Torque)
eq1 = D(1,1)-J2(1,1);
eq2 = D(2,1) - J2(1,3);

%part b - Finding the state space representation
sol = solve([eq1==0,eq2==0],[theta_1ddot,theta_2ddot] );
sol1 = sol.theta_1ddot;
sol2 = sol.theta_2ddot;

% Theta_1ddot and Theta_2ddot are used in ODE solver to define the state
% space representation

%part c
tspan = [0,10];
init_state = [30*2*pi/360 ; 45*2*pi/360 ; 0 ; 0 ];

[tsol,ysol] = ode45(@(t,y) ode_motion(t,y) , tspan , init_state );

%figure;
%plot(tsol,ysol(:,1),tsol,ysol(:,2),tsol,ysol(:,3),tsol,ysol(:,4));
plot(tsol,ysol,'-x');
ylabel('Theta and Omega  (rad & rad/s)');
xlabel('Time  (s)');
legend('theta_1','theta_2','theta_1dot','theta_2dot')