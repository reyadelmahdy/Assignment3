close all
clear
clc
%% Part 2
% a)
% Setting Constants
mElec = 9.11e-31; % Electron rest mass (kg)
mEff = 0.26*mElec; % Effective mass (kg)
kb = 1.381e-23; % Boltzmann's Constant (J/K)
T = 300; % Room Temperature (K)
L = 200e-9; % Length (m)
W = 100e-9; % Width (m)
q = 1.60218e-19; % Elemntary Charge (C)
Lb = 40e-9; % Box length
Wb = 40e-9; % Box width
s = 1e-9; % mesh spacing
nx = ceil(L/s + 1);
ny = ceil(W/s + 1);

cond = 1; % Conductive
res = 1e-2; % Resistive

% Generating the conductivity map
condMap = zeros(nx,ny);
for i = 1:nx
   for j = 1:ny
       if (i-1>0.5*(L-Lb)/s) && ((i-1)<0.5*(L+Lb)/s) && (((j-1)<Wb/s||(j-1)>(W-Wb)/s))
           condMap(i,j) = res;
       else
           condMap(i,j) = cond;
       end
   end
end

G = sparse(nx*ny);
Bc = zeros(1,nx*ny);

for i = 1:nx
    for j = 1:ny
        n = j +(i-1)*ny; %Mapping Eq
        
        %Setting Boundary conditions
        if i == 1 || i == nx
            G(n,n) = 1;
            if i == 1
                Bc(n) = 1;
            end
        elseif j == 1 || j == ny
            nxm = j+(i-2)*ny;
            nxp = j+i*ny;
            nyp = j+1+(i-1)*ny;
            nym = j-1+(i-1)*ny;
            
            %Resistances from the conductivity map
            rxm = (condMap(i,j)+condMap(i-1,j))/2;
            rxp = (condMap(i,j)+condMap(i+1,j))/2;
            
            if j == 1
                ryp = (condMap(i,j)+condMap(i,j+1))/2;
                G(n,n) = -(rxm+rxp+ryp);
                G(n,nyp) = ryp; 
            else
                rym = (condMap(i,j)+condMap(i,j-1))/2;
                G(n,n) = -(rxm+rxp+rym);
                G(n,nym) = rym;
            end
            G(n,nxm) = rxm;
            G(n,nxp) = rxp; 
            
        % Rest of the nodes
        else
            nxm = j+(i-2)*ny;
            nxp = j+i*ny;
            nym = j-1+(i-1)*ny;
            nyp = j+1+(i-1)*ny;
            
            % Getting the resistances from the conduction map
            rxm = (condMap(i,j)+condMap(i-1,j))/2;
            rxp = (condMap(i,j)+condMap(i+1,j))/2;
            ryp = (condMap(i,j)+condMap(i,j+1))/2;
            rym = (condMap(i,j)+condMap(i,j-1))/2;
            
            % Assigning the node equations
            G(n,n) = -(rxm+rxp+rym+ryp);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nym) = rym;
            G(n,nyp) = ryp; 
        end
    end
end

%solving and plotting
V = G\Bc';
mappedV = zeros(nx,ny);
for i = 1:nx
    for j = 1:ny
        n = j+(i-1)*ny;
        mappedV(i,j) = V(n);
    end
end

% a)

[xaxis, yaxis] = meshgrid(0:s:L,0:s:W);
figure(1)
surf(xaxis',yaxis',mappedV)
hold on
imagesc([0 L],[0 W],mappedV')
title('Potential Over the Region')
hold off

% b)
 
[Ex, Ey] = gradient(-mappedV',1e-9);

figure(2)
quiver(xaxis',yaxis',Ex',Ey')
xlim([0 L])
ylim([0 W])
title('E Field Quiver')

% c)
% Reusing the code from Q1 with Ex and Ey as inputs

% Finding the thermal velocity
Vth = sqrt((2*kb*T)/mEff); %Thermal Velocity in m/s
mt = 0.2e-12; % mean time (s)

% Modeling the electron paths

numPar = 10000; % Number of particles

% Assigning particle positions
posX = L.*rand(numPar,2);
posY = W.*rand(numPar,2);
% Reassigning electron position so none of them spawn inside the contacts
for i = 1:numPar  
    while ((posX(i,2) > 80e-9 && posX(i,2) < 120e-9) && (posY(i,2) > 60e-9 ||  posY(i,2) < 40e-9))
        posX(i,2) = L.*rand(1,1);
        posY(i,2) = W.*rand(1,1);
    end
end
posX(:,1) = posX(:,2);
posY(:,1) = posY(:,2);

spacialStep = sqrt(L^2+W^2)/1000;
stepTime = 0.5*spacialStep/Vth; %using a reduced steptime for this part

% Assigning each particle a random velocity (acceleration will be taken into account during the iteration)
Vx = randn(numPar,2)*sqrt((kb*T)/mEff);
Vy = randn(numPar,2)*sqrt((kb*T)/mEff);

% Displacement per step of each electron
dispX = stepTime*Vx(:,1);
dispY = stepTime*Vy(:,1);

colors = rand(numPar,3);

pScatter = 1 - exp(-stepTime/mt); % Probability of an electron scattering

% Positions before collision
collX = posX(:,1);
collY = posY(:,1);

ic = 0; % iteration counter for animated plot

for i = 1:1000
    for j = 1:numPar
        if((round(posX(j,1)*1e9) > 200) || (round(posY(j,1)*1e9)) > 100)
            Fx = Ex(floor(posY(j,1)*1e9),floor(posX(j,1)*1e9)).*q;
            Fy = Ey(floor(posY(j,1)*1e9),floor(posX(j,1)*1e9)).*q;
        elseif((round(posY(j,1)*1e9)) <=0 || (round(posX(j,1)*1e9) <= 0))
            Fx = Ex(ceil(posY(j,1)*1e9),ceil(posX(j,1)*1e9)).*q;
            Fy = Ey(ceil(posY(j,1)*1e9),ceil(posX(j,1)*1e9)).*q;
        else
            Fx = Ex(round(posY(j,1)*1e9),round(posX(j,1)*1e9)).*q;
            Fy = Ey(round(posY(j,1)*1e9),round(posX(j,1)*1e9)).*q;
        end
        ax = Fx/mEff; % Horizontal Acceleration
        ay = Fy/mEff; % Vertical Acceleration
        % Adding the extra velocity from the acceleration (x only)
        Vx(j,1) = Vx(j,1) + ax*stepTime;
        Vy(j,1) = Vy(j,1) + ay*stepTime;
        
        if(rand < pScatter)
            % generating new velocity/displacment values
            Vx(j,:) = randn(1)*sqrt((kb*T)/mEff) + ax*stepTime;
            Vy(j,:) = randn(1)*sqrt((kb*T)/mEff) + ay*stepTime;
            dispX(j,:) = stepTime*Vx(j,1);
            dispY(j,:) = stepTime*Vy(j,1);
            collX = posX(:,1);
            collY = posY(:,1);
        end
        
        if (posX(j,1)+dispX(j) > L)
            dispX(j,:) = stepTime*Vx(j,1);
            posX(j,2) = posX(j,1)+dispX(j)-L;
        elseif (posX(j,1)+dispX(j) < 0)
            dispX(j,:) = stepTime*Vx(j,1);
            posX(j,2) = posX(j,1)+dispX(j)+L;
        else
            dispX(j,:) = stepTime*Vx(j,1);
            posX(j,2) = posX(j,1)+dispX(j);
        end
        
        % Adding extra boundary conditions for particle collision with the
        % contacts
        
        % Top or bottom of the region
        if ((posY(j,1)+dispY(j) > W)|| (posY(j,1)+dispY(j) < 0))
            dispY(j) = -dispY(j);
            posY(j,2) = posY(j,1)+dispY(j);
        % Approaching a contact from the right
        elseif ((posX(j,1) < 80e-9) && (posX(j,1) + dispX(j) > 80e-9) && ((posY(j,1) > 60e-9) || (posY(j,1) < 40e-9)))
            dispX(j) = -dispX(j);
            posX(j,2) = posX(j,1)+dispX(j);
        % Approaching a contact from the left
        elseif ((posX(j,1)) > 120e-9 && (posX(j,1) + dispX(j) < 120e-9) && ((posY(j,1) > 60e-9) || (posY(j,1) < 40e-9)))
            dispX(j) = -dispX(j);
            posX(j,2) = posX(j,1)+dispX(j);
        % Top and bottom
        elseif (((posX(j,1) > 80e-9) && (posX(j,1) < 120e-9)) && ((posY(j,1) > 40e-9) && (posY(j,1) < 60e-9)) && ((posY(j,1) + dispY(j) < 40e-9) || (posY(j,1) + dispY(j) > 60e-9)))
            dispY(j) = -dispY(j);
            posY(j,2) = posY(j,1)+dispY(j);
        % Left Corner
        elseif ((posX(j,1) + dispX(j) == 80e-9) && ((posY(j,1) + dispY(j) == 40e-9) || (posY(j,1) + dispY(j) == 60e-9)))
            dispY(j) = -dispY(j);
            dispX(j) = -dispX(j);
            posY(j,2) = posY(j,1)+dispY(j);
            posX(j,2) = posX(j,1)+dispX(j);
        else
            posY(j,2) = posY(j,1)+dispY(j);
        end
    end
    if (ic == 0)
       figure(3);
       scatter(posX(1:10, 2),posY(1:10,2),1,colors(1:10,1));
       hold on;
       plot([80e-9,80e-9],[0,40e-9],'k',[80e-9,120e-9],[40e-9,40e-9],'k',[120e-9,120e-9],[40e-9,0],'k',[80e-9,120e-9],[0,0],'k')
       plot([80e-9,80e-9],[100e-9,60e-9],'k',[80e-9,120e-9],[60e-9,60e-9],'k',[120e-9,120e-9],[60e-9,100e-9],'k',[80e-9,120e-9],[100e-9,100e-9],'k')
       title('2-D Particle pathing');
       xlabel('Length');
       ylabel('Width');
       xlim([0 L]);
       ylim([0 W]);
    elseif (ic < 1000)
       title('2-D Particle pathing');
       scatter(posX(1:10, 2),posY(1:10,2),1,colors(1:10,1));
       hold on;
    else
       title('2-D Particle pathing');
       scatter(posX(1:10, 2),posY(1:10,2),1,colors(1:10,1));
       hold off; 
    end
    pause(0.001);
    posX(:,1) = posX(:,2);
    posY(:,1) = posY(:,2);
    ic = ic + 1;
end
% The code works as intended for the most part, however there is sometimes
% an index out of bounds error when calculating Fx and Fy, the mesh spacing
% from parts a and b (assignment 2) needs to be manually adapted to fit the
% number of electrons and iterations that need to be done for this part.
% The number of electrons is reduced in the graph due to me not having
% enough RAM and/or computing power to simulate 1000 electrons for 1000 steps.

