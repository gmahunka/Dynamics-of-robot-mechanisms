close all;
clear all;
warning off;
syms mu p d
%% #2
% Define the matrix
A = [ mu, -1,  0,  0;
      0, mu-1, -1, -1/2;
      0, 0, mu-1, -1;
     -d, (p+d), 0, mu];

% Compute the determinant
D = det(A);

% Display the result
disp('The determinant of the matrix is:');
simplify(D)

A2 = [ 1, 1, 1, 1, 1;
      4, 2, 0, -2, -4;
      6, 0, -2, 0, 6;
      4, -2, 0, 2, -4;
      1, -1, 1, -1, 1];

% Define the column vector
B = [1; -2; 1 + (d + p)/2; p/2; -d/2];

% Perform matrix multiplication
result = A2 * B;

% Display the result
disp('The result of the multiplication is:');
disp(simplify(result));

% Define the coefficients
a4 = p;
a3 = 2*d - p;
a2 = 4 - 4*d - p;
a1 = 8 + 2*d + p;
a0 = 4;

% Construct the Routh-Hurwitz matrix
H = [a3, a1, 0, 0;
     a4, a2, a0, 0;
     0, a3, a1, 0;
     0, a4, a2, a0];

% Compute all principal minors
minor_1x1 = det(H(1:1, 1:1));
minor_2x2 = det(H(1:2, 1:2));
minor_3x3 = det(H(1:3, 1:3));
minor_4x4 = det(H);

% Display the result
disp('The Routh-Hurwitz matrix is:');
disp(simplify(H));
h2 =@(p,d) - 8*d^2 + 8*d - 12*p;
h3 =@(p,d) - 16*d^3 - 8*d^2*p - 64*d^2 + 64*d - 16*p^2 - 96*p;
a2_plot=@(p,d) 4 - 4*d - p;
a3_plot=@(p,d) 2*d - p;

%Stability map
% figure(1);
% hold on;
% fimplicit(h2,[0, 1.5, 0, 1.1]);
% fimplicit(h3,[0, 1.5, 0, 1.1],'LineWidth',2);
% fimplicit(a2_plot,[0, 1.5,0,1.1]);
% fimplicit(a3_plot,[0, 1.5,0,1.1]);
% plot([0,0],[0,1]);
% title('Stability map');
% xlabel('p [-]'); ylabel('d [-]');
% grid on;
% text(0.05, 0.2, 'Stable Region', 'FontSize', 12, ...
%      'Color', 'red', 'FontWeight', 'bold', 'Rotation', 90);
% legend('det(H2) = 0', 'det(H3) = 0', 'a2 = 0', 'a3 = 0','a4=0');

%% #3
m=1.5;
tau=3e-3;
C=12;


% p_values = linspace(0, 0.3, 1000);  
% d_values = zeros(size(p_values));  
% options = optimset('Display', 'off'); 
% for i = 1:length(p_values)
%     p = p_values(i);
%     d_values(i) = fsolve(@(d) h3(p, d), 1, options);
% end
% p_d_pairs = [p_values' d_values'];
% figure(2);
% plot(p_values,d_values)

% We get these values by examining the plot
p1=0.1445;
d1=0.4247;

P1=p1*m/tau^2;
D1=d1*m/tau;

Q11=-(P1+D1/tau);
Q12=D1/tau;

Delta_1=C/P1;
Delta_v=4*C*tau^2/m;

%% #4
syms p d
alpha1=subs(sqrt(a1/a3),[p,d],[p1,d1]);

p1d=(-12 - d^2 + sqrt(d^4 - 16*d^3 - 40*d^2 + 64*d + 144)) / 4;
dd=linspace(0, 0.8284, 1000);
pp = zeros(size(dd));
f = zeros(size(dd));

for i = 1:length(dd)
     pp(i) = double(subs(p1d, d, dd(i))) ;
     alpha=(sqrt((8+2*dd(i)+pp(i))/(2*dd(i)-pp(i))));
     omega=atan2(2*alpha,alpha^2-1);
     if alpha==Inf
         omega=0;
     end
     f(i)=omega/(2*pi*tau);
end
figure;
subplot(1,2,2);
hold on;
xlabel('f [Hz]'); ylabel('d [-]');
plot(f,dd);
plot([0,45],[d1,d1],'r');
plot([0,45],[0.8284,0.8284], 'r--');
plot([0,45],[0,0], 'r--');
subplot(1,2,1);
hold on;
xlabel('p [-]'); ylabel('d [-]');
plot(pp,dd)
plot([0,0.2],[d1,d1]);
plot([0,0.2],[0.8284,0.8284], 'r--');
plot([0,0.2],[0,0], 'r--');

%% #5
syms rho p d

eq1 = rho^4 - 2*rho^3 + (1 + p/2 + d/2) * rho^2 + (p/2) * rho - d/2 == 0;
eq2 = 4*rho^4 - 4*rho^3 - p*rho + 2*d == 0;
eq3 = 6*rho^4 - (d + p + 2) * rho^2 - 3*d == 0;
sol = solve([eq1, eq2, eq3], [rho, p, d]);
sol_rho = vpa(sol.rho,4);
sol_p = vpa(sol.p,4);
sol_d = vpa(sol.d,4);
rho_min=sol_rho(2);
p2=sol_p(2);
d2=sol_d(2);

P2=p2*m/tau^2;
D2=d2*m/tau;

Q21=-(P2+D2/tau);
Q22=D2/tau;

Delta_2=C/P2;
%% #6
x0 = 10 * Delta_2;   % Initial position in meters
Tmax = 0.2;         % Simulation time in seconds
dt = 0.0001;         % Time step

t = 0:dt:Tmax;
N = length(t);
x = zeros(1, N);
v = zeros(1, N);
a = zeros(1, N);
Q2 = zeros(1, N);

% Initial condition
x(1) = x0;

% Simulation loop using numerical integration
for i = 2:N
    idx_tau = max(1, i - round(tau/dt));
    idx_2tau = max(1, i - round(2*tau/dt));
    Q2(i) = Q21 * x(idx_tau) + Q22 * x(idx_2tau);
    a(i) = Q2(i) / m;
    v(i) = v(i-1) + a(i) * dt;
    x(i) = x(i-1) + v(i) * dt;
end

% Plot results
% figure;
% plot(t, x, 'b', 'LineWidth', 1.5);
% xlabel('Time [s]'); ylabel('Position x(t) [m]');
% title('Position'); grid on;
% 
% figure;
% plot(t, v, 'r', 'LineWidth', 1.5);
% xlabel('Time [s]'); ylabel('Velocity v(t) [m/s]');
% title('Velocity'); grid on;
% 
% figure;
% plot(t, a, 'g', 'LineWidth', 1.5);
% xlabel('Time [s]'); ylabel('Acceleration a(t) [m/s^2]');
% title('Acceleration'); grid on;
% 
% figure;
% plot(t, Q2, 'k', 'LineWidth', 1.5);
% xlabel('Time [s]'); ylabel('Control Force Q2(t) [N]');
% title('Control Force'); grid on;
%% Additional tasks
x0 = 10*Delta_2;   % Initial position in meters
Tmax = 0.2;         % Simulation time in seconds
dt = 0.001;         % Time step

t = 0:dt:Tmax;
N = length(t);
x = zeros(1, N);
v = zeros(1, N);
a = zeros(1, N);
Q1 = zeros(1, N);
x(1) = x0;
% with friction
for i = 2:N
    idx_tau = max(1, i - round(tau/dt));
    idx_2tau = max(1, i - round(2*tau/dt));

    Q1(i) = Q11 * x(idx_tau) + Q12 * x(idx_2tau);
    friction = C * sign((x(i)-x(i-1))/tau);

    a(i) = (Q1(i) / m) - friction/m;
    v(i) = v(i-1) + a(i) * dt;
    x(i) = x(i-1) + v(i) * dt;
end

% figure;
% plot(t, x, 'b', 'LineWidth', 1.5);
% xlabel('Time [s]'); ylabel('Position x(t) [m]');
% title('With friction'); grid on;

%without friction
x = zeros(1, N);
v = zeros(1, N);
a = zeros(1, N);
Q1 = zeros(1, N);
x(1) = x0;

for i = 2:N
    idx_tau = max(1, i - round(tau/dt));
    idx_2tau = max(1, i - round(2*tau/dt));

    Q1(i) = Q11 * x(idx_tau) + Q12 * x(idx_2tau);
    
    a(i) = (Q1(i) / m);
    v(i) = v(i-1) + a(i) * dt;
    x(i) = x(i-1) + v(i) * dt;
end
% figure;
% plot(t, x, 'b', 'LineWidth', 1.5);
% xlabel('Time [s]'); ylabel('Position x(t) [m]');
% title('Without friction'); grid on;
