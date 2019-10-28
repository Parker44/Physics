tmax = 5;

m1 = 1;
m2 = 0.5;
d = 0.5;

[t6, r6] = mutual_orbit(tmax, 7, m1, m2, d);
[t7, r7] = mutual_orbit(tmax, 8, m1, m2, d);
[t8, r8] = mutual_orbit(tmax, 9, m1, m2, d);
[t9, r9] = mutual_orbit(tmax, 10, m1, m2, d);

r6 = r6(1, 1, :);
r6 = reshape(r6, [1, length(r6)]);
r7 = r7(1, 1, :);
r7 = reshape(r7, [1, length(r7)]);
r8 = r8(1, 1, :);
r8 = reshape(r8, [1, length(r8)]);
r9 = r9(1, 1, :);
r9 = reshape(r9, [1, length(r9)]);

clf; 
hold on;
plot(t6, r6, 'r-.o');
plot(t7, r7, 'g-.+');
plot(t8, r8, 'b-.*');

r7 = r7(1:2:end);
r8 = r8(1:4:end);
r9 = r9(1:8:end);

r67 = r6 - r7;
r78 = r7 - r8;
r89 = r8 - r9;

clf; 
hold on; 
plot(t6, r67, 'r-.o'); 
plot(t6, r78, 'g-.+');
plot(t6, r89, 'b-.*'); 

r78 = 4 * r78;
r89 = 16 * r89;

clf;
hold on; 
plot(t6, r67, 'r-.o'); 
plot(t6, r78, 'g-.+');
plot(t6, r89, 'b-.*'); 

