
tmax = 3;
level = 9;
cores_m = [1 1];
cores_r0 = [-0.20 0 0; 0.20 0 0];
cores_v0 = [0 sqrt(0.20)/0.4 0; 0 -sqrt(0.20)/0.4 0];
cores_ns = 3000;


[t, r] = toomre(tmax, level, cores_m, cores_r0, cores_v0, cores_ns);
