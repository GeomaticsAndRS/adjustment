%% Free Network Adjustment using Least Squares %%

clc;
clear all;
close all;
format long g


%% Observations and initial values for unknowns %%

directions = load('data\Directions.txt');
distances = load('data\Distances.txt');
points = load('data\Points.txt');

% Vector of observations
L = [distances; directions*pi/200];

% Vector of unknowns, extended by orientation unknowns
X_0 = [reshape([points(:,2) points(:,1)]',8,1); zeros(4,1)];

no_n = length(L);
no_u = length(X_0);
no_b = 3; %Number of constraints

r = no_n-no_u+no_b; %Redundancy


%% Stochastic Model %%

%VC Matrix of the observations
S_LL = diag([0.02^2*ones(5,1); (0.0015*pi/200)^2*ones(18,1)]);

%Theoretical standard deviation
sigma_0 = 1;

%Cofactor matrix of the observations
Q_LL = 1/sigma_0^2*S_LL;

%Weight matrix
P = inv(Q_LL);


%% Adjustment %%

syms x1 y1 x6 y6 x9 y9 x15 y15 o1 o6 o9 o15

%distances
L_sym(1) = euclideandist(x1,x15,y1,y15);
L_sym(2) = euclideandist(x6,x15,y6,y15);
L_sym(3) = euclideandist(x6,x1,y6,y1);
L_sym(4) = euclideandist(x9,x15,y9,y15);
L_sym(5) = euclideandist(x9,x6,y9,y6);

%directions
L_sym(6) = direction(x1,x15,y1,y15,o15);
L_sym(7) = direction(x9,x15,y9,y15,o15);
L_sym(8) = direction(x15,x1,y15,y1,o1);
L_sym(9) = direction(x6,x1,y6,y1,o1);
L_sym(10) = direction(x9,x6,y9,y6,o6);
L_sym(11) = direction(x15,x6,y15,y6,o6);
L_sym(12) = direction(x1,x6,y1,y6,o6);
L_sym(13) = direction(x15,x9,y15,y9,o9);
L_sym(14) = direction(x6,x9,y6,y9,o9);

L_sym(15:23) = L_sym(6:14);

%design matrix
J = jacobian(L_sym, [x6 y6 x9 y9 x1 y1 x15 y15 o6 o9 o1 o15]);

%Datum [6 9 1 15]
datum = [1 1 1 1]';

%Centroid
xc = sum(points(:,2).*datum)/sum(datum);
yc = sum(points(:,1).*datum)/sum(datum);

% Reduced Coordinates to the centroid
x_dash = points(:,2)-xc;
y_dash = points(:,1)-yc;

% Constraint matrix G
G = repmat(eye(2),1,4);

L_p = length(G);
Scale = sqrt(L_p);
G = G./Scale;

G_1 = sqrt(sum((y_dash.^2)+(x_dash.^2)));
G = [G; reshape([(y_dash./G_1) -(x_dash./G_1)]',1,8)];

datum = ones(3,1)*reshape((datum*ones(1,2))',1,8);
G = G.*datum;

%Final G-Matrix extended by 4x4 matrix of zeros for the omega's
G = [G zeros(3,4)];

%break-off condition
epsilon = 10^-10;
max_x_hat = Inf;

delta = 10^-12;
check2 = Inf;

iteration = 0;

while max_x_hat>epsilon ||  check2>delta

x6 = X_0(1);
y6 = X_0(2);
x9 = X_0(3);
y9 = X_0(4);
x1 = X_0(5);
y1 = X_0(6);
x15 = X_0(7);
y15 = X_0(8);
o6 = X_0(9);
o9 = X_0(10);
o1 = X_0(11);
o15 = X_0(12);
  
%Vector of reduced observations
L_0 = double(subs(L_sym));
l = L-L_0';

%Design Matrix
A = double(subs(J));

%Normal Matrix
N = A'*P*A;

%Vector of absolute values
n = A'*P*l;

%Extension of normal matrix
Next = [N G'; G zeros(3)];
next = [n;zeros(3,1)];

%Inversion of normal matrix / Cofactor matrix of the unknowns
Q_xx_ext = inv(Next);

%Solution of normal equation
x_hat = Q_xx_ext*next;

%Adjusted unknowns
X_hat = X_0+x_hat(1:end-3);

X_0 = X_hat;

%Check 1
max_x_hat = max(abs(x_hat(1:end-3)));

%Vector of residuals
v = A*x_hat(1:no_u)-l;
%Vector of adjusted observations
L_hat = L+v;

Psi(1) = euclideandist(x1,x15,y1,y15);
Psi(2) = euclideandist(x6,x15,y6,y15);
Psi(3) = euclideandist(x6,x1,y6,y1);
Psi(4) = euclideandist(x9,x15,y9,y15);
Psi(5) = euclideandist(x9,x6,y9,y6);

Psi(6) = direction(x1,x15,y1,y15,o15);
Psi(7) = direction(x9,x15,y9,y15,o15);
Psi(8) = direction(x15,x1,y15,y1,o1);
Psi(9) = direction(x6,x1,y6,y1,o1);
Psi(10) = direction(x9,x6,y9,y6,o6);
Psi(11) = direction(x15,x6,y15,y6,o6);
Psi(12) = direction(x1,x6,y1,y6,o6);
Psi(13) = direction(x15,x9,y15,y9,o9);
Psi(14) = direction(x6,x9,y6,y9,o9);

Psi(15:23) = Psi(6:14);

%Check 2
check2 = max(abs(L_hat-Psi'));

iteration = iteration+1;
end

%vector of residuals
v = A*x_hat(1:no_u)-l;

%vector of adjusted observations
L_hat = L+v;

%empirical reference standard deviation
s_0 = sqrt((v'*P*v)/r);

Q_xx = Q_xx_ext(1: end-3, 1:end-3);
S_XX_hat=s_0^2*Q_xx;
s_X = sqrt(diag(S_XX_hat)); %standard deviation of the adjusted unknows

Q_LL_hat = A*Q_xx*A';
S_LL_hat = s_0^2*Q_LL_hat;
s_L_hat = sqrt(diag(S_LL_hat)); %standard deviation of the adjusted observations

Q_vv = Q_LL-Q_LL_hat;
S_vv = s_0^2*Q_vv;
s_v = sqrt(diag(S_vv)); %standard deviation of the residuals


%% Plotting Adjusted Network %%

X_Plot=[X_hat(7); X_hat(5); X_hat(1); X_hat(3); X_hat(7)];
Y_Plot=[X_hat(8); X_hat(6); X_hat(2); X_hat(4); X_hat(8)];

figure
set(gcf,'name','Adjusted Network','numbertitle','off')

plot(X_Plot, Y_Plot, 'b');
hold on
scatter(X_Plot, Y_Plot, 'fill', 'b');
hold on
strim=num2str(15);
 text(X_Plot(1),Y_Plot(1),strim,'VerticalAlignment','bottom','FontSize',15)
 strim=num2str(1);
 text(X_Plot(2),Y_Plot(2),strim,'VerticalAlignment','bottom','FontSize',15)
  strim=num2str(6);
 text(X_Plot(3),Y_Plot(3),strim,'VerticalAlignment','bottom','FontSize',15)
  strim=num2str(9);
 text(X_Plot(4),Y_Plot(4),strim,'VerticalAlignment','bottom','FontSize',15)
title('Adjusted Network (Free Network)')
xlabel('Easting [m]')
ylabel('Northing [m]')
%legend('6 & 9 datum','Location','SouthOutside' );
hold all
     
%% Test Network %%

%empirical reference standard deviation - +/- 20-30% of sigma_0

if s_0 <= sigma_0+(0.3*sigma_0) & s_0 >= sigma_0-(0.3*sigma_0)
    display('Everything Ok')
else
    display('There was a problem')
end

%disp(s_0)