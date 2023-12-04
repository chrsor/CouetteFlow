tclear all
close all
clc
% @author: Christopher Sorg
% @date: July, 2022
% Under supervision of Dr. Daniel Veldman, FAU Erlangen-Nürnberg

% Source for viscosity: 
% https://en.wikipedia.org/wiki/List_of_viscosities
% If you choose your own parameters, recall:
% For our Stokes-setting you should choose a Newtonian fluid 
% with high viscosity

mu = 48.8;      % Default: Sunflower oil                                [mPa*s]
R1 = 5;         % Radius of the inner cylinder                          [m]
R2 = 6;         % Radius of the outer (R2>R1) cylinder                  [m]
Omega1 = 20;    % Angular velocity of the inner cylinder                [m/s]
Omega2 = 10;    % Angular velocity of the outer cylinder                [m/s]
meshr   = 25;   % Number of elements of the mesh in r-direction
meshphi = 50;   % Number of elements of the mesh in phi-direction
Phi = 2*pi;     % Angle in the mesh (please in [0,2*pi])
f = {@(r,phi) r^2+phi^2; @(r,phi) r};   % External forces (f_r,f_phi)   [N]

%%%% Algorithm %%%%
% Please do not change anything here
nelem = meshr*meshphi;
R = R2 - R1;

% Exact analytical homogenous solution for uphi with u_r = 0 (Th. 1.4)
uphiEx = R*(R2^2 * Omega2 - R1^2 * Omega1)/(R2^2 - R1^2) + ...
    (Omega1 - Omega2)*(R1^2 * R2^2)/(R*(R2^2 - R1^2));
fprintf("Analytical homogenous solution: %f\n", uphiEx);

% Step 1: Assign number to all nodes
fprintf("Starting Algorithm...\n");
Nodesr = 2*meshr+1;                    % Number of nodes in the r direction
Nodesphi = 2*meshphi+1;                % Number of nodes in the phi direction
nodesn = Nodesr*Nodesphi;              % Total number of nodes in velocity space
nodesM = (meshr+1)*(meshphi+1);        % Total number of nodes in pressure space

Rgrid  = linspace(R1,R2,Nodesr);       % Grid for r 
dr = Rgrid(2) - Rgrid(1);
Phigrid  = linspace(0,Phi,Nodesphi);   % Grid for phi
dphi = Phigrid(2) - Phigrid(1);

node_nmbrs = 1:nodesn;
node_nmbrs = reshape(node_nmbrs, Nodesr, Nodesphi);

node_nmbrsM = 1:nodesM;
node_nmbrsM = reshape(node_nmbrsM, [], meshphi+1);

% Step 2: Element list
elem_list = zeros(nelem, 9);        % Quadractic element with one node in the middle
elem_nmbrs = zeros(meshr, meshphi);
e = 0;
for ii = 1:meshr
    for jj = 1:meshphi
        e = e+1;
        elem_list(e, :) = [node_nmbrs(2*ii-1,2*jj-1), node_nmbrs(2*ii+1,2*jj-1), node_nmbrs(2*ii+1,2*jj+1), node_nmbrs(2*ii-1,2*jj+1), ...
            node_nmbrs(2*ii,2*jj-1), node_nmbrs(2*ii+1,2*jj), node_nmbrs(2*ii,2*jj+1), node_nmbrs(2*ii-1,2*jj), node_nmbrs(2*ii,2*jj)];
        elem_nmbrs(ii,jj) = e;                   
    end
end

elem_listM = zeros(nelem, 4);       % Rectangles
e = 0;
for ii = 1:meshr
    for jj = 1:meshphi
        e = e+1;
        elem_listM(e, :) = [node_nmbrsM(ii,jj), node_nmbrsM(ii+1,jj), node_nmbrsM(ii+1,jj+1), node_nmbrsM(ii,jj+1)];
        elem_nmbrs(ii,jj) = e;        
    end
end

% Step 3: Element matrices
syms r phi
p0  =@(x) (1-x)*(1-2*x);
p1  =@(x) (2*x-1)*x;
p12 =@(x) 4*(1-x)*x;
Ne = [p0(r)*p0(phi), p1(r)*p0(phi), p1(r)*p1(phi), p0(r)*p1(phi), ...
           p12(r)*p0(phi), p1(r)*p12(phi), p12(r)*p1(phi), p0(r)*p12(phi), p12(r)*p12(phi)]; % Velocity space
Me = [(1-r)*(1-phi), r*(1-phi), r*phi, (1-r)*phi];                                           % Pressure space

drNe = diff(Ne,r); 
dphiNe = diff(Ne,phi);
drMe = diff(Me,r);
dphiMe = diff(Me,phi);

% Transformation to the reference element
A11e1  = double(int(int((drNe.')*drNe     , r, [0, 1]), phi, [0, 1]));
A11e2  = double(int(int((dphiNe.')*dphiNe , r, [0, 1]), phi, [0, 1]));
A13e   = double(int(int((drNe.')*Me       , r, [0, 1]), phi, [0, 1]));
A23e   = double(int(int((dphiNe.')*Me     , r, [0, 1]), phi, [0, 1]));
A13eT  = double(int(int((Me.')*drNe       , r, [0, 1]), phi, [0, 1]));
A23eT  = double(int(int((Me.')*dphiNe     , r, [0, 1]), phi, [0, 1]));
Fe     = double(int(int(Ne.'              , r, [0, 1]), phi, [0, 1]));

% Step 4: Stiffness matrix A and load vector F
nDOFs = 2*nodesn+nodesM;
A = sparse(nDOFs,nDOFs);
for ii = 1:meshr
    for jj = 1:meshphi
        e = elem_nmbrs(ii,jj);
        nod = elem_list(e, :);
        rloc = (Rgrid(ii) + Rgrid(ii+1))/2; % Approximation
        A(nod,nod)               = A(nod,nod)               + ...
            mu*rloc*A11e1/(dr*dr)*(dr*dphi) + mu/rloc*A11e2/(dphi*dphi)*(dr*dphi);
        A(nodesn+nod,nodesn+nod) = A(nodesn+nod,nodesn+nod) + ...
            mu*rloc*A11e1/(dr*dr)*(dr*dphi) + mu/rloc*A11e2/(dphi*dphi)*(dr*dphi);
        nodM = elem_listM(e, :);
        A(nod, 2*nodesn+nodM)    = A(nod, 2*nodesn+nodM)    + (-rloc*A13e)/dr*(dr*dphi);
        A(nodesn + nod, 2*nodesn+nodM)    = A(nodesn + nod, 2*nodesn+nodM)    + (-A23e)/dphi*(dr*dphi);
        A(2*nodesn+nodM, nod)    = A(2*nodesn+nodM, nod)    + (-A13eT)/dr*(dr*dphi);
        A(2*nodesn+nodM, nodesn + nod)    = A(2*nodesn+nodM, nodesn + nod)    + (-A23eT)/dphi*(dr*dphi);
    end
end

F = sparse(nDOFs,1);
for ii = 1:meshr
    for jj = 1:meshphi
        e = elem_nmbrs(ii,jj);
        nod = elem_list(e, :);
        rloc = (Rgrid(ii) + Rgrid(ii+1))/2;       % Approximation
        philoc = (Phigrid(ii) + Phigrid(ii+1))/2; % Approximation
        flocR = f{1}(rloc, philoc);
        flocPhi = f{2}(rloc, philoc);
        F(nod,1)                     = F(nod,1)        + rloc*flocR*Fe;
        F(nodesn+nod,1)              = F(nodesn+nod,1) + rloc*flocPhi*Fe;
    end
end

% Step 5: Dirichlet BC
nodes_inner = node_nmbrs(1,:);
nodes_outer = node_nmbrs(end,:);
cdofs_inner = [nodes_inner, nodesn + nodes_inner];
cdofs_outer = [nodes_outer, nodesn + nodes_outer];
fdofs = setdiff(1:nDOFs, [cdofs_inner, cdofs_outer]);

F = F - A(:,nodesn+nodes_inner)*ones(length(nodes_inner),1)*Omega1 - ...
    A(:,nodesn+nodes_outer)*ones(length(nodes_outer),1)*Omega2;

Aff = A(fdofs, fdofs);
Ff  = F(fdofs);

% Step 6: Solve the system
UP = zeros(nDOFs,1);
UP(nodesn+nodes_inner) = Omega1;
UP(nodesn+nodes_outer) = Omega2;
UP(fdofs) = Aff \ Ff;

ur   = UP(1:nodesn);
uphi = UP(nodesn + (1:nodesn));
p    = UP(2*nodesn+1:end);

fprintf("Calculations done\n");

% Step 7: Ploting
fprintf("Now ploting..\n");

figure
surf(Rgrid, Phigrid, uphi(node_nmbrs).')
xlabel 'r'
ylabel 'phi'
zlabel 'u_\phi'

figure
surf(Rgrid, Phigrid, ur(node_nmbrs).')
xlabel 'r'
ylabel 'phi'
zlabel 'u_r'

figure
[RRR,PHI] = meshgrid(Rgrid, Phigrid);
X = RRR.*cos(PHI);
Y = RRR.*sin(PHI);
surf(X, Y, uphi(node_nmbrs).')
xlabel 'x'
ylabel 'y'
zlabel 'u_\phi'

figure
[RRR,PHI] = meshgrid(Rgrid, Phigrid);
X = RRR.*cos(PHI);
Y = RRR.*sin(PHI);
surf(X, Y, ur(node_nmbrs).')
xlabel 'x'
ylabel 'y'
zlabel 'u_r'