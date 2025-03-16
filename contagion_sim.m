clear
close all

%% Initial conditions
S_sim.N = 1e3; % Number of susceptible individuals
I_sim.N = 1;   % Number of infected individuals
R_sim.N = 0;   % Number of recovered individuals

S = S_sim.N;
I = I_sim.N;
R = R_sim.N;

L = 20; % Simulation area size
length = 1.1; % Infection radius
infection_rate = 0.1; % Probability of infection
recovery_rate = 0.004; % Probability of recovery

% Initialize positions and velocities randomly
S_sim.X = L*(rand(S_sim.N,2)*2 - 1);
I_sim.X = L*(rand(I_sim.N,2)*2 - 1);
R_sim.X = L*(rand(R_sim.N,2)*2 - 1);

S_sim.V = rand(S_sim.N,2)*2 - 1;
I_sim.V = rand(I_sim.N,2)*2 - 1;
R_sim.V = rand(R_sim.N,2)*2 - 1;

%% Simulation settings
ts = 0;
tf = 200;
h = 0.1;
t = ts:h:tf;

S_sim.X_str = cell(1,numel(t));
I_sim.X_str = cell(1,numel(t));
R_sim.X_str = cell(1,numel(t));

S_sim.N_str = zeros(1,numel(t));
I_sim.N_str = zeros(1,numel(t));
R_sim.N_str = zeros(1,numel(t));

S_sim.X_str{1} = S_sim.X;
I_sim.X_str{1} = I_sim.X;
R_sim.X_str{1} = R_sim.X;

S_sim.N_str(1) = S_sim.N;
I_sim.N_str(1) = I_sim.N;
R_sim.N_str(1) = R_sim.N;

%% Particle simulation loop
for step = 1:numel(t)
    % Calculate distances between susceptible and infected individuals
    distance = zeros(S_sim.N,I_sim.N);
    for i = 1:S_sim.N
        for j = 1:I_sim.N
            distance(i,j) = sqrt(sum((S_sim.X(i,:) - I_sim.X(j,:)).^2));            
        end
    end

    % Determine infections
    [i_idx,j_idx] = find(distance < length);
    i_idx = sort(unique(i_idx),'descend');
    i_str = [];
    
    for k = 1:numel(i_idx)
        i = i_idx(k);
        if rand < infection_rate && ~any(i_str == i)
            I_sim.X = [I_sim.X; S_sim.X(i,:)];
            I_sim.V = [I_sim.V; S_sim.V(i,:)];
            S_sim.X(i,:) = [];
            S_sim.V(i,:) = [];
            S_sim.N = S_sim.N - 1;
            I_sim.N = I_sim.N + 1;
            i_str = [i_str i];
        end
    end
    
    % Determine recoveries
    N_tmp = sort(1:I_sim.N,'descend');
    for l = 1:numel(N_tmp)
        j = N_tmp(l);
        if rand < recovery_rate
            R_sim.X = [R_sim.X; I_sim.X(j,:)];
            R_sim.V = [R_sim.V; I_sim.V(j,:)];
            I_sim.X(j,:) = [];
            I_sim.V(j,:) = [];
            I_sim.N = I_sim.N - 1;
            R_sim.N = R_sim.N + 1;
        end
    end
    
    % Apply boundary conditions (reflection)
    S_sim.V(abs(S_sim.X) - L > 0) = -S_sim.V(abs(S_sim.X) - L > 0);
    I_sim.V(abs(I_sim.X) - L > 0) = -I_sim.V(abs(I_sim.X) - L > 0);
    R_sim.V(abs(R_sim.X) - L > 0) = -R_sim.V(abs(R_sim.X) - L > 0); 
    
    % Update positions
    S_sim.X = S_sim.X + S_sim.V*h;
    I_sim.X = I_sim.X + I_sim.V*h;
    R_sim.X = R_sim.X + R_sim.V*h;
    
    % Store results
    S_sim.X_str{step} = S_sim.X;
    I_sim.X_str{step} = I_sim.X;
    R_sim.X_str{step} = R_sim.X;
    S_sim.N_str(step) = S_sim.N;
    I_sim.N_str(step) = I_sim.N;
    R_sim.N_str(step) = R_sim.N;
    
    disp([step,S_sim.N,I_sim.N,R_sim.N])
end

%% Visualization settings
fig = figure('Position',[100 100 1400 600]);
reset(groot)
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(groot,'defaultTextInterpreter','latex')
set(groot,'defaultLegendInterpreter','latex')
set(groot,'defaultAxesFontSize',13)
set(groot,'defaultLineLineWidth',1.5)
set(gca, 'LooseInset', get(gca, 'TightInset'));

tiledlayout(1,2)

%% Visualization of simulation
nexttile(1)
S_sim.points = scatter(S_sim.X_str{1}(:,1),S_sim.X_str{1}(:,2),50,'b','filled');
hold on
I_sim.points = scatter(I_sim.X_str{1}(:,1),I_sim.X_str{1}(:,2),50,'r','filled');
R_sim.points = scatter(R_sim.X_str{1}(:,1),R_sim.X_str{1}(:,2),50,'g','filled');
axis([-1 1 -1 1]*L)
daspect([1 1 1])
box on
xticks([])
yticks([])

%% Comparison with SIR model
nexttile(2)
Shandle = plot(t(1),S_sim.N_str(1),'b-');
hold on
Ihandle = plot(t(1),I_sim.N_str(1),'r-');
Rhandle = plot(t(1),R_sim.N_str(1),'g-');

beta = 0.0003;
gamma = 0.05;
S_str = zeros(1,numel(t));
I_str = zeros(1,numel(t));
R_str = zeros(1,numel(t));
S_str(1) = S;
I_str(1) = I;
R_str(1) = R;

for i = 1:numel(t)
    dS = -beta*S*I;
    dI = beta*S*I - gamma*I;
    dR = gamma*I;
    S = S + dS*h;
    I = I + dI*h;
    R = R + dR*h;
    S_str(i) = S;
    I_str(i) = I;
    R_str(i) = R;
end

plot(t,S_str,'b--')
plot(t,I_str,'r--')
plot(t,R_str,'g--')
xlabel('step $t$')
ylabel('population')
legend("S","I","R")

%% Drawing Differential Equations
text(tf*0.65,600,'$\frac{dS}{dt} = -\beta S I$','FontName','latex','FontSize',20)
text(tf*0.65,450,'$\frac{dI}{dt} = \beta S I - \gamma I$','FontName','latex','FontSize',20)
text(tf*0.65,300,'$\frac{dR}{dt} = \gamma I$','FontName','latex','FontSize',20)

%% Create simulation video
video = VideoWriter("contagion_sim.avi",'Uncompressed AVI');
open(video)

for step = 1:numel(t)
    S_sim.points.XData = S_sim.X_str{step}(:,1);
    S_sim.points.YData = S_sim.X_str{step}(:,2);
    I_sim.points.XData = I_sim.X_str{step}(:,1);
    I_sim.points.YData = I_sim.X_str{step}(:,2);
    R_sim.points.XData = R_sim.X_str{step}(:,1);
    R_sim.points.YData = R_sim.X_str{step}(:,2);

    Shandle.XData = t(1:step);
    Shandle.YData = S_sim.N_str(1:step);
    Ihandle.XData = t(1:step);
    Ihandle.YData = I_sim.N_str(1:step);
    Rhandle.XData = t(1:step);
    Rhandle.YData = R_sim.N_str(1:step);

    frame = getframe(gcf);
    writeVideo(video,frame)
end

close(video)