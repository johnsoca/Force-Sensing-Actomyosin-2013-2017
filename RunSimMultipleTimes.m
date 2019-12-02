% RunSimMultipleTimes.m
%  4/9/16
% Callie J Miller
clear all
close all
runSims='/Users/mill29ca/Desktop/Code for Lance/Code for Martin Simulations';
E=100;
cd(runSims);
newFolder=runSims;
VariableTestingSimMulti(runSims,newFolder,E);
save data.mat

% E=10;
% cd(runSims);
% newFolder='/Users/calliemiller/Desktop/Adam work/MultiSim Anisotropic Circle E=100 kPa top E=10 kPa sides/1';
% VariableTestingSimMulti(runSims,newFolder,E);
% save data.mat
% clearvars -except runSims E
% close all
% 
% cd(runSims);
% newFolder='/Users/calliemiller/Desktop/Adam work/MultiSim Anisotropic Circle E=100 kPa top E=10 kPa sides/2';
% VariableTestingSimMulti(runSims,newFolder,E);
% save data.mat
% clearvars -except runSims E
% close all
% 
% cd(runSims);
% newFolder='/Users/calliemiller/Desktop/Adam work/MultiSim Anisotropic Circle E=100 kPa top E=10 kPa sides/3';
% VariableTestingSimMulti(runSims,newFolder,E);
% save data.mat
% clearvars -except runSims E
% close all
% 
% cd(runSims);
% newFolder='/Users/calliemiller/Desktop/Adam work/MultiSim Anisotropic Circle E=100 kPa top E=10 kPa sides/4';
% VariableTestingSimMulti(runSims,newFolder,E);
% save data.mat
% clearvars -except runSims E
% close all
% 
% cd(runSims);
% newFolder='/Users/calliemiller/Desktop/Adam work/MultiSim Anisotropic Circle E=100 kPa top E=10 kPa sides/5';
% VariableTestingSimMulti(runSims,newFolder,E);
% save data.mat
% clearvars -except runSims E
% close all
% 
% cd(runSims);
% newFolder='/Users/calliemiller/Desktop/Adam work/MultiSim Anisotropic Circle E=100 kPa top E=10 kPa sides/6';
% VariableTestingSimMulti(runSims,newFolder,E);
% save data.mat
% clearvars -except runSims E
% close all
% 
% cd(runSims);
% newFolder='/Users/calliemiller/Desktop/Adam work/MultiSim Anisotropic Circle E=100 kPa top E=10 kPa sides/7';
% VariableTestingSimMulti(runSims,newFolder,E);
% save data.mat
% clearvars -except runSims E
% close all
% 
% cd(runSims);
% newFolder='/Users/calliemiller/Desktop/Adam work/MultiSim Anisotropic Circle E=100 kPa top E=10 kPa sides/8';
% VariableTestingSimMulti(runSims,newFolder,E);
% save data.mat
% clearvars -except runSims E
% close all
% 
% cd(runSims);
% newFolder='/Users/calliemiller/Desktop/Adam work/MultiSim Anisotropic Circle E=100 kPa top E=10 kPa sides/9';
% VariableTestingSimMulti(runSims,newFolder,E);
% save data.mat
% clearvars -except runSims E
% close all
% 
% cd(runSims);
% newFolder='/Users/calliemiller/Desktop/Adam work/MultiSim Anisotropic Circle E=100 kPa top E=10 kPa sides/10';
% VariableTestingSimMulti(runSims,newFolder,E);
% save data.mat
% clearvars -except runSims E
% close all
% 
% E=5;
% cd(runSims);
% newFolder='/Users/calliemiller/Desktop/Adam work/MultiSim Anisotropic Circle E=100 top E=5 kPa sides/1';
% VariableTestingSimMulti(runSims,newFolder,E);
% save data.mat
% clearvars -except runSims E
% close all
% 
% cd(runSims);
% newFolder='/Users/calliemiller/Desktop/Adam work/MultiSim Anisotropic Circle E=100 top E=5 kPa sides/2';
% VariableTestingSimMulti(runSims,newFolder,E);
% save data.mat
% clearvars -except runSims E
% close all
% 
% cd(runSims);
% newFolder='/Users/calliemiller/Desktop/Adam work/MultiSim Anisotropic Circle E=100 top E=5 kPa sides/3';
% VariableTestingSimMulti(runSims,newFolder,E);
% save data.mat
% clearvars -except runSims E
% close all
% 
% cd(runSims);
% newFolder='/Users/calliemiller/Desktop/Adam work/MultiSim Anisotropic Circle E=100 top E=5 kPa sides/4';
% VariableTestingSimMulti(runSims,newFolder,E);
% save data.mat
% clearvars -except runSims E
% close all
% 
% cd(runSims);
% newFolder='/Users/calliemiller/Desktop/Adam work/MultiSim Anisotropic Circle E=100 top E=5 kPa sides/5';
% VariableTestingSimMulti(runSims,newFolder,E);
% save data.mat
% clearvars -except runSims E
% close all
% 
% cd(runSims);
% newFolder='/Users/calliemiller/Desktop/Adam work/MultiSim Anisotropic Circle E=100 top E=5 kPa sides/6';
% VariableTestingSimMulti(runSims,newFolder,E);
% save data.mat
% clearvars -except runSims E
% close all
% 
% cd(runSims);
% newFolder='/Users/calliemiller/Desktop/Adam work/MultiSim Anisotropic Circle E=100 top E=5 kPa sides/7';
% VariableTestingSimMulti(runSims,newFolder,E);
% save data.mat
% clearvars -except runSims E
% close all
% 
% cd(runSims);
% newFolder='/Users/calliemiller/Desktop/Adam work/MultiSim Anisotropic Circle E=100 top E=5 kPa sides/8';
% VariableTestingSimMulti(runSims,newFolder,E);
% save data.mat
% clearvars -except runSims E
% close all
% 
% cd(runSims);
% newFolder='/Users/calliemiller/Desktop/Adam work/MultiSim Anisotropic Circle E=100 top E=5 kPa sides/9';
% VariableTestingSimMulti(runSims,newFolder,E);
% save data.mat
% clearvars -except runSims E
% close all
% 
% cd(runSims);
% newFolder='/Users/calliemiller/Desktop/Adam work/MultiSim Anisotropic Circle E=100 top E=5 kPa sides/10';
% VariableTestingSimMulti(runSims,newFolder,E);
% save data.mat
% clearvars -except runSims E
% close all
% % 
% % E=0;
% % cd(runSims);
% % newFolder='/Users/calliemiller/Desktop/Adam work/MultiSim Isotropic Circle E=0 kPa/1';
% % VariableTestingSimMulti(runSims,newFolder,E);
% % save data.mat
% % clearvars -except runSims E
% % close all
% % 
% % cd(runSims);
% % newFolder='/Users/calliemiller/Desktop/Adam work/MultiSim Isotropic Circle E=0 kPa/2';
% % VariableTestingSimMulti(runSims,newFolder,E);
% % save data.mat
% % clearvars -except runSims E
% % close all
% % 
% % cd(runSims);
% % newFolder='/Users/calliemiller/Desktop/Adam work/MultiSim Isotropic Circle E=0 kPa/3';
% % VariableTestingSimMulti(runSims,newFolder,E);
% % save data.mat
% % clearvars -except runSims E
% % close all
% % 
% % cd(runSims);
% % newFolder='/Users/calliemiller/Desktop/Adam work/MultiSim Isotropic Circle E=0 kPa/4';
% % VariableTestingSimMulti(runSims,newFolder,E);
% % save data.mat
% % clearvars -except runSims E
% % close all
% % 
% % cd(runSims);
% % newFolder='/Users/calliemiller/Desktop/Adam work/MultiSim Isotropic Circle E=0 kPa/5';
% % VariableTestingSimMulti(runSims,newFolder,E);
% % save data.mat
% % clearvars -except runSims E
% % close all
% % 
% % cd(runSims);
% % newFolder='/Users/calliemiller/Desktop/Adam work/MultiSim Isotropic Circle E=0 kPa/6';
% % VariableTestingSimMulti(runSims,newFolder,E);
% % save data.mat
% % clearvars -except runSims E
% % close all
% % 
% % cd(runSims);
% % newFolder='/Users/calliemiller/Desktop/Adam work/MultiSim Isotropic Circle E=0 kPa/7';
% % VariableTestingSimMulti(runSims,newFolder,E);
% % save data.mat
% % clearvars -except runSims E
% % close all
% % 
% % cd(runSims);
% % newFolder='/Users/calliemiller/Desktop/Adam work/MultiSim Isotropic Circle E=0 kPa/8';
% % VariableTestingSimMulti(runSims,newFolder,E);
% % save data.mat
% % clearvars -except runSims E
% % close all
% % 
% % cd(runSims);
% % newFolder='/Users/calliemiller/Desktop/Adam work/MultiSim Isotropic Circle E=0 kPa/9';
% % VariableTestingSimMulti(runSims,newFolder,E);
% % save data.mat
% % clearvars -except runSims E
% % close all
% % 
% % cd(runSims);
% % newFolder='/Users/calliemiller/Desktop/Adam work/MultiSim Isotropic Circle E=0 kPa/10';
% % VariableTestingSimMulti(runSims,newFolder,E);
% % save data.mat
