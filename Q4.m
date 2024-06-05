%Basic initialization
nx = 200; %assumning the domain is rechtangular 
dx  = 2/nx; %assuming discretication is rechtangular
tfinal = 3.0;
C = 0.7;
T = 0;
Count = 1;

%set up intial conditions - note that its bad to let h be zero 
U0 = zeros(3,nx,nx);
for i = 1:nx
    for j = 1:nx
        if (abs(-1+i*dx)<= 0.5 & abs(-1+j*dx)<= 0.5)
            U0(1,i,j) = 2;
        else 
            U0(1,i,j) = 1;
        end
    end
end

%surf(linspace(-1,1,nx),linspace(-1,1,nx),squeeze(U0(1,:,:)))
%axis([-1 1 -1 1 1 3])

U = U0;
outputMovie3Dh = zeros(1,1,nx,nx);

outputChart1 = zeros(1,nx,nx);
outputChart2 = zeros(1,nx,nx);
outputChart3 = zeros(1,nx,nx); 

outputChart1 = U(1,:,:);

while T<tfinal
    
    Unew = U;
    LXG = getXMaxGlobal(U,nx);
    LYG = getYMaxGlobal(U,nx);
    dt = 0.5*C*min(dx/LXG,dx/LYG);  
    T = T + dt
    Count = Count + 1;
    
    %ok cool now I can acutally do LF!
    for i = 1:nx
        for j = 1:nx
            
            %function gets U for pretty much every location.
            %allows for easy handeling of boundary values **if the function
            %is correct**
            [Up, Ue, Uw,Un,Us] = getU(U,i,j,nx);
    
            %Get the maximum speeds 
            Le = getMaxLX(Ue,Up);
            Lw = getMaxLX(Uw,Up);
            Ln = getMaxLY(Un,Up);
            Ls = getMaxLY(Us,Up);

            %Get the LF fluxes 
            Fe = 0.5*(Fx(Up)+Fx(Ue)) - 0.5*Le*(Ue - Up); 
            Fw = 0.5*(Fx(Up)+Fx(Uw)) - 0.5*Lw*(Up - Uw);

            Fn = 0.5*(Fy(Up)+Fy(Un)) - 0.5*Le*(Un - Up);
            Fs = 0.5*(Fy(Up)+Fy(Us)) - 0.5*Lw*(Up - Us);

            %perform the update
            Unew(:,i,j) = Up - (dt/dx)*(Fe-Fw+Fn-Fs);   
        end 
    end 
  
    U = Unew;
    outputMovie3Dh(Count,:,:,:) = U(1,:,:); 

    if (round(T,2) == 1.35)
        outputChart2 = U(1,:,:);
    elseif (round(T,2) == 3)
        outputChart3 = U(1,:,:);
    end 
    
    %imagesc(squeeze(U(1,:,:)))
end   

surf(linspace(-1,1,nx),linspace(-1,1,nx),squeeze(U(1,:,:)))
axis([-1 1 -1 1 0 3])


%FLUX FUNCTIONS 
function [flux] = Fx(U)
    q1 = U(1);
    q2 = U(2);
    q3 = U(3);
    flux = [q2 ; (q2^2 / q1)+0.5*(q1^2); q2*q3/q1];
end 

function [flux] = Fy(U)
    q1 = U(1);
    q2 = U(2);
    q3 = U(3);
    flux = [q3 ; q2*q3/q1 ; (q3^2 / q1)+0.5*(q1^2)];
end 

%this function gets the maximum speed
function [L] = getMaxLX(Ul,Ur)
    U = [Ul,Ur];
    speeds = zeros(2,1);
    for i = 1:2
        h = U(1,i);
        hu = U(2,i);
        hv = U(3,i);

        speeds(i) =  abs((hu/h)) + sqrt(h);
    end 
    L = max(speeds);
end 

function [L] = getMaxLY(Un,Us)
    U = [Un,Us];
    speeds = zeros(2,1);
    for i = 1:2
        h = U(1,i);
        hu = U(2,i);
        hv = U(3,i);

        speeds(i) =  abs((hv/h)) + sqrt(h);
    end 
    L = max(speeds);
end 

%THESE FUNCTIONS GET THE GLOBAL MAXIMUM SPEEDS TO SET DT 
function [Lg] = getXMaxGlobal(U,nx)
    speeds = zeros(nx*nx,1);
    for i = 1:nx
        for j = 1:nx
            h = U(1,i,j);
            hu = U(2,i,j);
            hv = U(3,i,j);

            speeds((i-1)*nx + j) = abs((hu/h)) + sqrt(h);
        end 
    end 
    Lg = max(speeds);
end 

function [Lg] = getYMaxGlobal(U,nx)
    speeds = zeros(nx*nx,1);
    for i = 1:nx
        for j = 1:nx
            h = U(1,i,j);
            hu = U(2,i,j);
            hv = U(3,i,j);

            speeds((i-1)*nx + j) = abs((hv/h)) + sqrt(h);
        end 
    end 
    Lg = max(speeds);
end 
%making getting U a function because.... it sucks 
%THIS IS A MASSIVE POTENTIAL SOURCE OF ERROR. IF BOUNDARIES ARE MISBEHAVING
%THIS IS ALMOST CERTAINLY THE ISSUE. TRY FIXED BOUNDARIES IF THINGS DON'T
%WORK
function [Up, Ue, Uw,Un,Us] = getU(U,i,j,nx)
    Up = U(:,i,j); %simple one 
    if (i < nx)
        Ue = U(:,i+1,j);
    else 
        Ue = [1 0 0 ; 0 -1 0; 0 0 1] * U(:,i,j);
    end

     if (i > 1)
        Uw = U(:,i-1,j);
    else 
        Uw = [1 0 0 ; 0 -1 0; 0 0 1] * U(:,i,j);
    end
    
    %set north and south 
    if (j < nx)
        Us = U(:,i,j+1);
    else 
        Us = [1 0 0 ; 0 1 0; 0 0 -1] * U(:,i,j);
    end

     if (j > 1)
        Un = U(:,i,j-1);
    else 
        Un = [1 0 0 ; 0 1 0; 0 0 -1] * U(:,i,j);
    end
end 

