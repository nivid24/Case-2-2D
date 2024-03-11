% Parameters
L = 0.1;         % Length of the rod (meters)
T_initial = 100; % Initial temperature (degrees Celsius)
alpha = 1e-4;    % Thermal diffusivity

% Discretization
Nx = 100;        % Number of spatial points
Nt = 1000;       % Number of time steps
dx = L / (Nx - 1);
dt = 1e-2;

% Initial conditions
x = linspace(0, L, Nx)';
T = ones(Nx, 1) * T_initial;

% Time-stepping loop
for n = 1:Nt
    % Interior points (finite difference scheme)
    T(2:end-1) = T(2:end-1) + alpha * dt / dx^2 * (T(3:end) - 2*T(2:end-1) + T(1:end-2));
    
    % Boundary conditions (fixed temperature at both ends)
    T(1) = T_initial;
    T(end) = T_initial;
    
    % Plot the temperature distribution at every 100 time steps
    if mod(n, 100) == 0
        plot(x, T);
        title(['Time = ', num2str(n*dt), ' seconds']);
        xlabel('Position (meters)');
        ylabel('Temperature (degrees Celsius)');
        ylim([T_initial-10, T_initial+10]);
        pause(0.01);
    end
    end