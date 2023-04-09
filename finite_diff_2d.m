clear
clc
close all

N_surface = 300;
E_surface = 4000;
S_surface = 250;
W_surface = 250;
q_dot = 0.25; %W/m^3
k_air = 33.33E-3;

for N = 80
    A = eye(N*N);
    T = [];
    b = [];

    %boundary conditions (b-vector)
    b(1,1) = 0.25*(W_surface)+0.25*(N_surface); %top left node
    b((N*N)-(N-1),1) = 0.25*(W_surface)+0.25*(S_surface); %bottom left node
    b(N,1) = 0.25*(N_surface)+0.25*(E_surface); %top right node
    b(N*N,1) = 0.25*(E_surface)+0.25*(S_surface); %bottom right node

    for i = N+1:N:(N*N)-((N-1)+N) %left side non-corner
        b(i,1) = 0.25*W_surface;
    end

    for i = 2:1:N-1 %top side non-corner
        b(i,1) = .25*N_surface;
    end

    for i = N*2:N:(N*N)-N %right side non-corner
        b(i,1) = 0.25*E_surface;
    end

    for i = (N*N)-N+2:1:(N*N)-1 %bottom side non-corner
        b(i,1) = 0.25*S_surface;
    end
% 
%     for i = N+1:1:N*2
%         b(i,1) = b(i,1) + q_dot/k_air;
%     end
% 
%     for i = (N*N)-N+1-(2*N):1:(N*N)-(2*N)
%         b(i,1) = b(i,1) + q_dot/k_air;
%     end
%         


    %%%NODE DEPENDENCIES%%%
    %%%interior nodes%%%
    interior = N+2:1:2*N-1;
    for i = 1:N-3
        interior = [interior (interior(1)+i*N):1:(interior(N-2)+i*N)];
    end

    for i = interior
        A(i,i+1) = -0.25;
        A(i,i+N) = -0.25;
        A(i,i-N) = -0.25;
        A(i,i-1) = -0.25;
    end

    %%%CORNER NODES%%%
    %%%top left node%%%
    A(1,2) = -0.25;
    A(1,N+1) = -0.25;

    %%%top right node%%%
    A(N,N-1) = -0.25;
    A(N, N*2) = -0.25;

    %%%bottom left node%%%
    A(1+(N*(N-1)),1+(N*(N-1))-N) = -0.25;
    A(1+(N*(N-1)),2+(N*(N-1))) = -0.25;

    %%%bottom right node%%%
    A(N*N,(N*N)-1) = -0.25;
    A(N*N,(N*N)-N) = -0.25;

    %%%OTHER NODES%%%
    %%%top row non-corner nodes%%%
    for i = 2:N-1
        A(i,i-1) = -0.25;
        A(i,i+1) = -0.25;
        A(i,i+N) = -0.25;
    end


    %%%bottom row non-corner nodes%%%
    for i = N*N-(N-2):(N*N)-1
        A(i,i-1) = -0.25;
        A(i,i+1) = -0.25;
        A(i,i-N) = -0.25;
    end

    %%%right side non-corner nodes%%%
    for i = 2*N:N:(N*N)-N
        A(i,i-N) = -0.25;
        A(i,i-1) = -0.25;
        A(i,i+N) = -0.25;
    end

    %%%left side non-corner nodes%%%
    for i = 1+N:N:(N*N)-N-(N-1)
        A(i,i+1) = -0.25;
        A(i,i-N) = -0.25;
        A(i,i+N) = -0.25;
    end

    T = A\b;
    T = reshape(T,N,N)';

    figure();
    heatmap(T,'Colormap',jet)
    title("Temperature(K) Heatmap - " + N + "x" + N + " mesh")
    xlabel('Node Column')
    ylabel('Node Row')

    surf(T)


end
