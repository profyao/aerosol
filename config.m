% file header
const.header_MI1B2T_url = 'ftp://l5eil01.larc.nasa.gov/MISR/MI1B2T.003/';
const.header_MI1B2T_filename = 'MISR_AM1_GRP_TERRAIN_GM_P';
const.header_MIL2ASAE_url = 'ftp://l5eil01.larc.nasa.gov/MISR/MIL2ASAE.002/';
const.header_MIL2ASAE_filename = 'MISR_AM1_AS_AEROSOL_P';
const.header_MIL2ASAF = 'ftp://l5eil01.larc.nasa.gov/MISR/MIL2ASAF.001/';
const.MIANSMT_SS_filename = 'MISR_AM1_SMART_TOA_RHO_ATM_SS_F02_0009.hdf';
const.MIANSMT_MS_filename = 'MISR_AM1_SMART_TOA_RHO_ATM_MS_F02_0009.hdf';
const.MIANSMT_TDIFF_filename = 'MISR_AM1_SMART_TDIFF_F02_0009.hdf';
const.MIANSMT_EDIFF_filename = 'MISR_AM1_SMART_BOA_EDIFF_F02_0009.hdf';
const.header_MIANCAGP_url1 = 'ftp://l5eil01.larc.nasa.gov/MISR/MIANCAGP.001/1999.11.07/';
const.header_MIANCAGP_url2 = 'ftp://l5eil01.larc.nasa.gov/MISR/MIANCAGP.001/1999.11.08/';
const.header_MIANCAGP_filename = 'MISR_AM1_AGP_P';

% MISR cameras parameters
const.Cam_DF = 1;
const.Cam_CF = 2;
const.Cam_BF = 3;
const.Cam_AF = 4;
const.Cam_AN = 5;
const.Cam_AA = 6;
const.Cam_BA = 7;
const.Cam_CA = 8;
const.Cam_DA = 9;
const.Cam_Dim = 9;
const.Cam_Name = {'DF','CF','BF','AF','AN','AA','BA','CA','DA'};
const.Cam_Dim = length(const.Cam_Name);

% MISR spatial resolutions
const.r275 = 275;
const.r1100 = 1100;
const.r4400 = 4400;
const.r8800 = 8800;
const.r17600 = 17600;
const.XDim_r1100 = 128;
const.YDim_r1100 = 512;
const.XDim_r4400 = 32;
const.YDim_r4400 = 128;
const.XDim_r8800 = 16;
const.YDim_r8800 = 64;
const.XDim_r17600 = 8;
const.YDim_r17600 = 32;
const.r = const.r4400; % default resolution for retrieval
const.XDim_r = const.XDim_r4400; % default X dimension
const.YDim_r = const.YDim_r4400; % default Y dimension

% Number of subregions in a region
const.RegSize = const.r / const.r1100;
% Scale factor to the 17.6-km standard region
const.RegScale = const.r17600 / const.r;

% XDim is the number of rows in a block, depending on the resolution
% YDim is the number of columns in a block,  depending on the resolution

% MISR bands parameters
const.Band_Blue = 1;
const.Band_Green = 2;
const.Band_Red = 3;
const.Band_NIR = 4;
const.Band_Dim = 4;
const.NChannel = const.Band_Dim*const.Cam_Dim;
const.Band_Name = {'BlueBand','GreenBand','RedBand','NIRBand'};
const.Band_Used = [1,1,1,1];
const.Channel_Used = logical(kron(const.Band_Used',ones(const.Cam_Dim,1)));


const.Band_Radiance = {'Blue Radiance/RDQI','Green Radiance/RDQI','Red Radiance/RDQI','NIR Radiance/RDQI'};
const.Config_rdqi1 = 1;
const.Config_c_lambda = [5.67e-6, 1.04e-4, 4.89e-5,3.94e-6];
const.Config_spectral_corr_matrix = [1.0106,-0.0057,-0.0038,-0.0011;
    -0.0080,1.0200,-0.0086,-0.0034;
    -0.0060,-0.0048,1.0145,-0.0036;
    -0.0048,-0.0033,-0.0136,1.0217];

const.sample_size = const.RegSize*const.RegSize;
const.Config_min_het_subr_thresh = const.sample_size/4;
const.min_cam_used=2;
const.Config_first_eigenvalue_for_eofs = 1;
const.Config_eigenvector_variance_thresh = 0.95;

% MISR aerosol model parameters
const.Model_ComponentDim = 21;
const.Model_MixtureDim = 74;
const.Model_NumComponent = 3;
const.Model_Pressure = [607.95,1013.25];

const.Model_mu0Grid = 0.2:0.01:1.0;
const.Model_muGrid = [0.31,0.32,0.33,0.34,0.35,0.47,0.48,0.49,0.5,0.51,0.66,0.67,0.68,0.69,0.7,0.71,0.84,0.85,0.86,0.87,0.88,0.89,0.9,0.95,0.96,0.97,0.98,0.99,1];
const.Model_ScatterAngleGrid = [-1,0,2.5,5,7.5,10,12.5,15,17.5,20,22.5,25,27.5,30,32.5,35,37.5,40,42.5,45,47.5,50,52.5,55,57.5,60,62.5,65,67.5,70,72.5,75,77.5,80,82.5,...
    85,87.5,90,92.5,95,97.5,100,102.5,105,107.5,110,112.5,115,117.5,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,...
    145,146,147,148,149,150,152.5,155,157.5,160,162.5,165,167.5,170,172.5,175,176,177,178,179,180,181];
const.Model_OpticalDepthGrid = [0,0.05,0.1,0.2,0.4,0.6,0.8,1,1.5,2,3,4,6];
const.Model_OpticalDepthLen = length(const.Model_OpticalDepthGrid);

const.Model_AOTGridGap = 0.025;
const.Model_OpticalDepthFinerGrid = 0:const.Model_AOTGridGap:3;
const.Model_OpticalDepthFinerGridLen = length(const.Model_OpticalDepthFinerGrid);

const.Component_Particle = [1,2,3,6,8,14,19,21];
const.Component_Num = length(const.Component_Particle);

const.Config_albedo_thresh_land = 0.015;

const.cols = [0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.6350    0.0780    0.1840];

const.str_kf = {'','_k-factor'};
const.str_dy = {'','_dynamic_component'};
const.str_par = {'','_par'};
