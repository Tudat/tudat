results = importdata('stateHistory.txt');
plot(results(:,1),results(:,2:4)/1e3); grid on;
xlabel('Time [seconds since J2000]');
ylabel('Position [km]');
legend('x','y','z');
