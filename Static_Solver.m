function Static_Solver
% Function produces the plots for the Static domains (solve the PDE + numerical modes at large t=T)

%% Parameters

a=0.2;
b=1;
d_1 = 0.004;
d_2 = 0.1;
L=1.2; 
tf = 1001;

%% Solve PDE

function zvec = kenetics(w)
    zvec = [a-w(1) + w(1)^2 * w(2); b - w(1)^2 * w(2)];
end

function u0 = pdeic(x)
    u0 = [(a+b); b/(a+b)^2]+[3;0]*(x<0.1);%[(a+b); b/(a+b)^2]+0.05*(ra5nd(1)-0.5)*[1;1];
end

function [c,f,s] = pdefun(x,t,w,dw)

    c = [1;1];
    f = [d_1; d_2].* dw;
    s = kenetics(w);
end

function [pl,ql,pr,qr] = pdebc(xl,ul,xr,ur,t)
    pl = [0;0];
    pr = [0;0];
    ql = [1;1];
    qr = [1;1];
end

x = linspace(0,L,1001);
t = linspace(0,500,tf);

sol = pdepe(0, @pdefun, @pdeic, @pdebc, x, t);
u = sol(:,:,1);
v = sol(:,:,2);


%% Plot solution
[X,T] = meshgrid(x,t);

figure('color','white')
subplot(1,2,1)
pcolor(X,T,u)
xlabel('Position $x$', Interpreter='latex', FontSize=18)
ylabel('Time $t$', Interpreter='latex', FontSize=18)
title('System evolution', 'Interpreter','latex', FontSize=15)
shading interp

%% Plot modes at t = tf
u_coef = extract_gen_fourier_coff(u(end,:), x,15);
v_coef = extract_gen_fourier_coff(v(end,:), x,15);

rho_range = [d_2*(b-a)-d_1*(a+b)^3 - sqrt((d_2*(b-a)-d_1*(a+b)^3)^2-4*d_1*d_2*(a+b)^4),...
    d_2*(b-a)-d_1*(a+b)^3 + sqrt((d_2*(b-a)-d_1*(a+b)^3)^2-4*d_1*d_2*(a+b)^4)]/(2*d_1*d_2*(a+b));


k_range = sqrt(rho_range)*L/pi;

subplot(1,2,2)
plot_modes(u_coef, v_coef, k_range)

sgtitle(['System evolution and mode composition for  $\ell = ', num2str(L), '$'], Interpreter="latex", FontSize=20)

end

%% Helper function for plot of dominant modes + instabilities
function plot_modes(u_coef,v_coef, k_range)

magnitudes = sqrt(u_coef.^2 + v_coef.^2);
mode_labels = 1:length(magnitudes);
hold on
area([k_range(1), min([length(magnitudes), k_range(2)])],...
    [100,100]);
plot(mode_labels, magnitudes, '.', MarkerSize=40)
ylim([0, max([max(magnitudes)*1.1,0.001])])
xlabel('Wavenumber $k$',Interpreter="latex", FontSize=18)
ylabel('Magnitude $||\mathbf{W}_k||$',Interpreter="latex",FontSize=18)
xticks(mode_labels)
legend({"Unstable $k$ range", ''}, Interpreter="latex", fontSize = 15)
title('Large time mode composition', Interpreter="latex", FontSize=15);
hold off

end