clear; close all;

%parameters 
rho = 0.107;
D = 0.255; 
C_m = 10e7

%domain and step
Lx = 100;
Ly = 100;
Lz = 100; 

Nx = 101;
Ny = 101;
Nz = 101;
Nt = 365;

dx = Lx/(Nx-1);
dy = Ly/(Ny-1);
dz = Lz/(Nz-1);
dt = 1;

%field variables
Tn = zeros(Nx, Ny, Nz); %temperature
x = linspace(0, Lx, Nx); %x-distance
y = linspace(0, Ly, Ny); %y-distance
z = linspace(0, Lz, Nz); %z-distance

[X,Y,Z] = meshgrid(x,y,z);

%insulate the boundary in x,y,z direction
K= ones(Nx, Ny, Nz);
K([1 end], :, :) = 0.0001;
K(:, [1 end], :) = 0.0001;
K(:, :, [1 end]) = 0.0001;

% %initial conditions
Tn(:,:,:) = 0;
t=0;
 
% %insulated boundary
Tn(1,:,:) = Tn(2,:,:);
Tn(end,:,:) = Tn(end-1,:,:);
 
Tn(:,1,:) = Tn(:,2,:);
Tn(:,end,:) = Tn(:,end-1,:);
 
Tn(:,:,1) = Tn(:,:,2);
Tn(:,:, end) = Tn(:,:,end-1);

figh = figure;

total_tb = linspace(0,Nt, Nt);

Tn(30,30,30) = 1;
for n=1:Nt
    Tc = Tn; 
    t=t+dt;
    for i=2:Nx-1
        for j=2:Ny-1
             for k=2:Nz-1
                Tn(i,j,k) = Tc(i,j,k) + D*(dt/dx^2)*...
                    ((Tc(i+1,j,k)-2*Tc(i,j,k) + Tc(i-1,j,k))/dx/dx +...
                     (Tc(i,j+1,k)-2*Tc(i,j,k) + Tc(i,j-1,k))/dy/dy +...
                     (Tc(i,j,k+1)-2*Tc(i,j,k) + Tc(i-1,j,k-1))/dz/dz)+ rho*Tc(i,j,k);
             end
        end
    end
   
   if t == 305
       uc_growth = sum(Tn,'all');
       figure (1)
       slice(X,Y,Z, Tn, [30 ],[],[])
       figure (2)
       slice(X,Y,Z, Tn, [],[30],[])
       figure (3)
       slice(X,Y,Z, Tn, [],[],[30])
   end
   
end

