function memUsed_linux
[~,pid] = system('pgrep MATLAB');
[~,mem_usage] = system(['cat /proc/' strtrim(pid) '/status | grep VmSize']);
fprintf("%i MB\n", round(str2num(strtrim(extractAfter(extractBefore(mem_usage, ' kB'), ':'))) / 1000));