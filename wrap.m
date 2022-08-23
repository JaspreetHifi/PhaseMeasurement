function wang = wrap( ang )

wang = ang;

idx = find(wang > pi );
wang(idx) = wang(idx) - 2*pi;
idx = find(wang < -pi);
wang(idx) = wang(idx) + 2*pi;


