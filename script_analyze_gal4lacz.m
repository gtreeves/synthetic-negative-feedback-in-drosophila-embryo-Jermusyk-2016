% script_analyze_gal4lacz
%
% This script is designed to run the pipeline code on Gal4/LacZ embryos.


pthbase = 'D:\GDrive\Students\Etika Goyal\Confocal\2018-10-05 first try\';
pths = {'lacZ'
	'MS2'
	'Zld_MS2'};
gtypes = {'Gal4x4,UAS-lacZ'
	'Gal4x4,UAS-MS2-lacZ'
	'Gal4x4,UAS-Zld-MS2-lacZ'};

Soln = [];
for i = 1:1%length(pths)
	soln = run_analyze_AP([pthbase,pths{i}],[1 2],{'h3','GCN4x4_lacZ'},...
		'Gal4x4,UAS-Zld-MS2-lacZ');
	
	Soln = [Soln;soln];
end