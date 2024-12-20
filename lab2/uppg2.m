% lab 2 upg 2

%-----------------------------------------------

% DELUPPGIFT b mha upg a

%-----------------------------------------------

clear all; close all; clc; clf;

L = 3.40; % [m]
T0 = 300; % [K]
TL = 450; % [K] 

% Söker x = 1.40 [m]

n = 34; % Nr of steps (adjusted to "hit" x = 1.40)

tolerance = 1e-4;

prev_prev_T_index = 2;
prev_T_index = 1;

% Iteration-loop for "Finita Differensmetoden (för Randvillkorsproblem)"

disp('For x = 1.40, T(x) = : ')

for i = 1:10
    % Step size
    h = L/n;
    
    Q  =@(x) 285 .* exp(- ((x - (L .* 0.5)).^2));
    
    % C1 =@(x) -h + 6 * (3 + x/6);                    % For T_i-1
    % C2 =@(x) -2 * 6 * (3 + x/6);     % For T_i
    % C3 =@(x) h + 6 * (3 + x/6);     % For T_i+1

    C1 =@(x) (x/6 + 3)/(h^2) - 1/(12*h);                    % For T_i-1
    C2 =@(x) (x/6 + 3)/(h^2)*(-2);     % For T_i
    C3 =@(x) (x/6 + 3)/(h^2) + 1/(12*h);     % For T_i+1
    
    
    
    % Creates x_i (importanto to have aftr functions)
    x = 0:h:L;
    
    A = zeros(n-1); % Autoprogrammed to be a squared of the size of ( ) w zeros. square matrix as standard.
    A = sparse(A);
    % b-vector
    b = zeros([n-1, 1]);
    
    % Create matrix A with diagonals accordingly
    for row = 1:n-1
        if row - 1 > 0
            A(row, row - 1) = C1(x(row + 1));
        end
        
        A(row, row)     = C2(x(row + 1));
        
        if row + 1 < n
            A(row, row + 1) = C3(x(row + 1));
        end
        
        % b(row) = -6 * h^2 * Q(x(row + 1));
        b(row) = -Q(x(row + 1));
    end
    
    % disp('Matrix A: ')
    % disp(A)
    
    b(1) = b(1) - C1(x(2)) .* T0;
    b(n-1) = b(n-1) - C3(x(n - 1)) .* TL;
    
    T = A\b;

    % disp('T: ')
    % disp(T)
    
    x_index = floor(n * 1.40 / L) + 1;
    % disp(T(x_index))
    T_index = T(x_index);
    E_trunk = abs(T_index - prev_T_index);
    convergence = (T_index - prev_T_index) / (prev_T_index - prev_prev_T_index);
    disp("   T: " + T_index + "   prev_T: " + prev_T_index + "   prev_prev_T:" + prev_prev_T_index + "   E_trunk: " + E_trunk + "   Convergence: " + convergence + "   step_size: " + h) % convergence should in Finita Differensmetoden ca equal to 2^p, p = 2. 
    prev_prev_T_index = prev_T_index;
    prev_T_index = T_index;

    % exit condition
    if E_trunk < tolerance
        if i >= 5
            disp("iterations: " + i + "   E_trunk: " + E_trunk);
            break
        end
    end

    n = n * 2; % doubles nr of steps per iteration; will always hit x = 1.40
    % disp(['Nr of steps n: ', n])
end

plot(x, [T0,T',TL], '-');
hold on
plot([1.4, 1.4], [min(T), max(T)])
hold on
x_index = floor(n * 0.5 * 1.40 / L) + 1;
plot([0, L], [T(x_index), T(x_index)]);