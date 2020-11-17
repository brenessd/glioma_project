clear; close all;

%parameters 
rho = 0.107;
D = 0.255; 
C_m= 10e7;

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

%set initial conditions
Tn(:,:,:) = 0;
t=0;
 
% %insulated boundary
Tn(1,:,:) = Tn(2,:,:);
Tn(end,:,:) = Tn(end-1,:,:);
 
Tn(:,1,:) = Tn(:,2,:);
Tn(:,end,:) = Tn(:,end-1,:);
 
Tn(:,:,1) = Tn(:,:,2);
Tn(:,:, end) = Tn(:,:,end-1);

%define radiation therapy
num =30;
dose = 1.8; 
alpha = 0.027;
therapy = exp(-alpha*(num*dose*(1- (dose/10))));

figh = figure;

%initialize tumor population
Tn(30,30,30) = 1;

%method of finite differences
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
                     
               
                if t > 257 && t <= 262
                    Tn(i,j,k) = Tn(i,j,k) - therapy*Tc(i,j,k);
                end
                if t > 264 && t <= 269
                    Tn(i,j,k) = Tn(i,j,k) - therapy*Tc(i,j,k);
                end
                if t > 271 && t <= 276
                    Tn(i,j,k) = Tn(i,j,k) - therapy*Tc(i,j,k);
                end
                if t > 278 && t <= 283
                    Tn(i,j,k) = Tn(i,j,k) - therapy*Tc(i,j,k);
                end
                if t > 285 && t <= 290
                    Tn(i,j,k) = Tn(i,j,k) - therapy*Tc(i,j,k);
                end
                if t > 292 && t <= 297
                    Tn(i,j,k) = Tn(i,j,k) - therapy*Tc(i,j,k);
                end
                if t > 299 && t <= 304
                    Tn(i,j,k) = Tn(i,j,k) - therapy*Tc(i,j,k);
                end
             end 
        end
    end
    
   %measure tumor burden after treatment is done 
   if t == 305
       uc_growth = sum(Tn,'all');
       figure (2)
       slice(X,Y,Z, Tn, [30 ],[],[])
       figure (3)
       slice(X,Y,Z, Tn, [],[30],[])
       figure (4)
       slice(X,Y,Z, Tn, [],[],[30])
   end
   
   %Get frames for movie
   slice(X,Y,Z, Tn, [30 ],[],[30])
   title(sprintf('Uncontrolled Tumor growth at t = %f Days', t));
   movieVector(n) = getframe(figh);
end


%Step 5: Make & Save Movie
%Create a VideoWriter object and set properties
myWriter = VideoWriter('tumor-burden','MPEG-4');   %create an .mp4 file
myWriter.FrameRate = 20;

%Open the VideoWriter object, write the movie, and close the file
open(myWriter);
writeVideo(myWriter, movieVector);
close(myWriter);

disp('DONE!')
