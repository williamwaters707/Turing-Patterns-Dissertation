function turing_space_plot
% This function plots the Turing space for Schnakenburg kenitics given
% fixed diffusion parameters


%% Paramaters + a,b ranges
d1 = 0.004;
d2 = 0.1;

a=0:0.0005:0.5;
b=0:0.0005:2.5;

[A,B] = meshgrid(a,b);

%% Check conditions
Condt1 = ((A+B).^3-B+A)>0;
Condt2 = (d2*(B-A)-d1*(A+B).^3)>0;
Condt3 = ((d1*(A+B).^3-d2*(B-A)).^2 - 4*d1*d2*(A+B).^4)>0;

final = Condt1 .* Condt2 .* Condt3;

%% Plot
figure('color','white')
pcolor(A,B,1-double(final))
fontsize(gca, 13, 'points')
xlabel('$a$', Interpreter='latex', FontSize=30)
ylabel('$b$', Interpreter='latex', FontSize=30)

shading flat
end


