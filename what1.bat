
@echo off
echo p	N	T	S
for /L %%x in (1,1,5) do for /L %%p in (1, 1, 12) do for %%N in (1,2,4,7,10,22,47,100,220,470,1000,2200,4700,10000,22000,47000,100000,220000,470000,1000000,2200000,4700000,10000000,22000000,47000000,100000000,220000000,470000000,1000000000,2200000000,4700000000,10000000000) do mpiexec -n %%p App311.exe %%N