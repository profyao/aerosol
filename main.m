addpath(genpath('~/projects/aerosol/src'));
config

preprocess_batch('Baltimore',const,4400) %r
par_aod_retri_batch('CD','Baltimore',4400,4,const); %r,core
%time = get_time('2011.06.02',16,60934,59,const);

%plot_result('scatter_4band',4400,const,'LA','CD-random')
plot_result('scatter',4400,const,'Baltimore','CD')
%legend({'MAP','MCMC','MISR'})
plot_result('overlay',1100,const,'2011.06.02',16,60934,59,'Baltimore','CD',jet(256))
%plot_result('overlay',4400,const,'2011.06.04',14,60963,59,'Baltimore','CD',hot(256))
%plot_result('overlay',4400,const,'2011.07.20',16,61633,59,'Baltimore','CD',hot(256))
%plot_result('overlay',4400,const,'2011.07.22',14,61662,59,'Baltimore','MISR',hot(256))
%plot_result('overlay',4400,const,'2011.07.29',15,61764,59,'Baltimore','CD-random',hot(256))

%plot_result('resid',4400,const,'2011.06.02',16,60934,59,'CD-random',1)
%plot_result('stability',4400, const,'Baltimore',false,'CD-random')
%plot_result('post_dist',const,'Baltimore',1)
%plot_result('stability',const,'Baltimore',1,'CD-random')
%plot_result('theta',4400,const,'2011.06.02',16,60934,59,'MCMC',jet)
%plot_result('theta',4400,const,'2011.06.04',14,60963,59,'CD',jet(256))
%plot_result('theta',4400,const,'2011.07.20',16,61633,59,'MCMC',jet(256))
%plot_result('theta',4400,const,'2011.07.22',14,61662,59,'CD-random',jet(256))
%plot_result('theta',4400,const,'2011.07.29',15,61764,59,'CD-random',jet(256))

%plot_result('trace',4400,const,'Baltimore','raw')