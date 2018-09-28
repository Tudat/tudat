unperturbed = importdata('unperturbedStateHistory.txt');
perturbed = importdata('perturbedStateHistory.txt');
plot(perturbed(:,1),(perturbed(:,2:4)-unperturbed(:,2:4))/1e3); grid on;
xlabel('Time [seconds since J2000]');
ylabel('Position Difference [km]');
legend('x','y','z','Location','NW');