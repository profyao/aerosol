Code developer: Shijing Yao, yaoshijing@gmail.com

This library contains all matlab functions needed to perform aerosol retrieval using MCMC and MAP methods. The hierarchical Bayesian model is described in this paper: 

https://www.stat.berkeley.edu/~binyu/ps/papers2012/WangJYJ12.pdf

To use this library, you may want to first create a project folder on your machine, say ‘aerosol/’, and add this library into that folder as ‘aerosol/src/‘. Then ‘mv main.m ..’ to perform various jobs.

1. Download and preprocess MISR data, which includes data aggregation, EOF computing, and SMART interpolation. The downloaded files will be automatically stored at ‘aerosol/cache/data’

run preprocess_batch(Location,const,r) in ‘main.m’

You need to specify Location (a string, currently support ‘Baltimore’, ‘Beijing’ and ’LA’). r by default is 4400(m), but you can also specify it as 1100, 17600 or other MISR resolutions.

2. Perform MCMC or MAP AOD retrievals

run par_aod_retri_batch(Method,Location,r,core,const) to perform AOD retrieval.

Method: {‘CD’, ‘CD-noprior’, ‘CD-random’, ‘CD-random-noprior’, ‘MCMC’}, ‘CD’ means coordinate descent, and ‘CD-random’ means coordinate-wise coordinate descent. ‘-noprior’ means Gaussian Markov random field prior and Dirichlet prior is turned off.

core: number of cores, usually set to 4 for a single laptop.

The results will be stored at ‘aerosol/cache/result’

3. Visualize results and validate with AERONET measurement
Prepare AERONET data using scripts aeronet.R. You need to download AERONET data of particular location from their website ‘http://aeronet.gsfc.nasa.gov/cgi-bin/type_piece_of_map_opera_v2_new?level=3' to ‘aerosol/aeronet/raw’. Or please email code developer to request processed AERONET data for fast demo.

To get scatter plot of MAP using CD with 4400m resolution, run plot_result('scatter',4400,const,'Baltimore','CD')

To plot out AOD retrievals on google map with 1100m resolution , run plot_result('overlay',1100,const,'2011.06.02',16,60934,59,'Baltimore','CD',jet(256))

To learn other functions or have questions about the code, please email code developer for details.

