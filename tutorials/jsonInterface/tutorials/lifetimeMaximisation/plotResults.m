%% READ
day = 1:365;
t0 = zeros(size(day));
dhpdt = t0;
for i = 1:length(t0)
    results = importdata(sprintf('outputs/day%i.txt',i));
    t0(i) = results(1,1);
    dhpdt(i) = (results(2,2) - results(1,2))/(results(2,1) - t0(i));
end

%% PLOT
figure;
plot(day,-dhpdt/1e3*86400);
grid on;
xlabel('Day initial epoch [year 2000]');
ylabel('Rate of decrease of the perigee altitude [km/day]');
