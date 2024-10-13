@echo off
echo p	t	s
for /L %%k in (1, 1, 10) do for /L %%i in (1, 1, 12) do mpiexec -n %%i App311.exe 100000000