clear all;
close all;
clc;
%load('young_old_RIOts68.mat');
load('AC_20120917_fMRI_new.mat');
%load('AC_20120917_ROIts_68_Bold');
%load('/home/shruti/Desktop/Himadri/SC-FC_fMRI/AC_20120917_ROIts_68.mat');
sig_sum=zeros(10,10);
%%%%% %% Default Mode Network %% %%%%%%%%

%for k=26:49

%% PCC/Pcu, Rostral Middle Frontal, IPC
% gc_lprec=gc_young(:,24,k);
% gc_lprec=reshape(gc_lprec,661,1);
% 
% gc_rprec=gc_young(:,58,k);
% gc_rprec=reshape(gc_rprec,661,1);
% 
% gc_lsfc=gc_young(:,27,k);
% gc_lsfc=reshape(gc_lsfc,661,1);
% 
% gc_rsfc=gc_young(:,61,k);
% gc_rsfc=reshape(gc_rsfc,661,1);
% 
% gc_lmfc=gc_young(:,26,k);
% gc_lmfc=reshape(gc_lmfc,661,1);
% 
% gc_rmfc=gc_young(:,60,k);
% gc_rmfc=reshape(gc_rmfc,661,1);
% 
% gc_lipc=gc_young(:,7,k);
% gc_lipc=reshape(gc_lipc,661,1);
% 
% gc_ripc = gc_young(:,41,k);
% gc_ripc=reshape(gc_ripc,661,1);
% 
% gc_lpcc= gc_young(:,22,k);
% gc_lpcc=reshape(gc_lpcc,661,1);
% 
% gc_rpcc=gc_young(:,56,k);
% gc_rpcc=reshape(gc_rpcc,661,1);
% 
% 
% % gc_lst=gc(:,63,:);
% % gc_lst=reshape(gc_lst,661,49);
% % 
% % gc_lmt=gc(:,48,:);
% % gc_lmt=reshape(gc_lmt,661,49);
% 
% 
% X(1,:)=gc_lpcc;
% X(2,:)=gc_rpcc;
% X(3,:)=gc_lmfc;
% X(4,:)=gc_rmfc;
% X(5,:)=gc_lipc;
% X(6,:)=gc_ripc;
% % X(9,:)=gc_lpcc;
% % X(10,:)=gc_rpcc;


%%% Salience Network %%%%
gc_lracc=AC_20120917_ROIts_68(:,25);
gc_lracc=reshape(gc_lracc,661,1);

gc_rracc=AC_20120917_ROIts_68(:,59);
gc_rracc=reshape(gc_rracc,661,1);

gc_lcacc=AC_20120917_ROIts_68(:,2);
gc_lcacc=reshape(gc_lcacc,661,1);

gc_rcacc=AC_20120917_ROIts_68(:,36);
gc_rcacc=reshape(gc_rcacc,661,1);

gc_linsula=AC_20120917_ROIts_68(:,34);
gc_linsula=reshape(gc_linsula,661,1);

gc_rinsula=AC_20120917_ROIts_68(:,68);
gc_rinsula=reshape(gc_rinsula,661,1);



% X(1,:)=gc_lracc;
% X(2,:)=gc_rracc;
% X(3,:)=gc_lcacc;
% X(4,:)=gc_rcacc;
% X(5,:)=gc_linsula;
% X(6,:)=gc_rinsula;
% X(9,:)=gc_lpcc;
% % X(10,:)=gc_rpcc;


%%% Central Executive Network %%%%
gc_lrmfc=AC_20120917_ROIts_68(:,22);
gc_lrmfc=reshape(gc_lrmfc,661,1);

gc_rrmfc=AC_20120917_ROIts_68(:,56);
gc_rrmfc=reshape(gc_rrmfc,661,1);

gc_lcmfc=AC_20120917_ROIts_68(:,27);
gc_lcmfc=reshape(gc_lcmfc,661,1);


gc_rcmfc=AC_20120917_ROIts_68(:,61);
gc_rcmfc=reshape(gc_rcmfc,661,1);

gc_lspc=AC_20120917_ROIts_68(:,28);
gc_lspc=reshape(gc_lspc,661,1);


gc_rspc=AC_20120917_ROIts_68(:,7);
gc_rspc=reshape(gc_rspc,661,1);


gc_lec=AC_20120917_ROIts_68(:,5);
gc_lec=reshape(gc_lec,661,1);


gc_rec=AC_20120917_ROIts_68(:,39);
gc_rec=reshape(gc_rec,661,1);

gc_lacc=AC_20120917_ROIts_68(:,25);
gc_lacc=reshape(gc_lacc,661,1);


gc_racc=AC_20120917_ROIts_68(:,59);
gc_racc=reshape(gc_racc,661,1);



X(1,:)=gc_lrmfc;
X(2,:)=gc_rrmfc;
X(3,:)=gc_lcmfc;
X(4,:)=gc_rcmfc;
X(5,:)=gc_lspc;
X(6,:)=gc_rspc;
X(7,:)=gc_lec;
X(8,:)=gc_rec;
X(9,:)=gc_lacc;
X(10,:)=gc_racc;


%% MVGC demo
%
% Demonstrates typical usage of the MVGC toolbox on generated VAR data for a
% 5-node network with known causal structure (see <var5_test.html |var5_test|>).
% Estimates a VAR model and calculates time- and frequency-domain
% pairwise-conditional Granger causalities (also known as the "causal graph").
% Also calculates Seth's causal density measure [2].
%
% This script is a good starting point for learning the MVGC approach to
% Granger-causal estimation and statistical inference. It may serve as a useful
% template for your own code. The computational approach demonstrated here will
% make a lot more sense alongside the reference document> [1], which we
% _strongly recommend_ you consult, particularly Section 3 on design principles
% of the toolbox. You might also like to refer to the <mvgc_schema.html schema>
% of MVGC computational pathways - <mvgc_schema.html#3 algorithms> |A<n>| in
% this demo refer to the algorithm labels listed there - and the
% <mvgchelp.html#4 Common variable names and data structures> section of the
% Help documentation.
%
% *_FAQ:_* _Why do the spectral causalities look so smooth?_ 
%
% This is because spectral quantities are calculated from the estimated VAR,
% rather than sampled directly. This is in accordance with the MVGC design
% principle that all causal estimates be based on the <mvgc_demo.html#6
% estimated VAR model> for your data, and guarantees that spectral causalities
% <mvgc_demo.html#10 integrate correctly> to time-domain causality as theory
% requires. See [1] for details.
% 
% *_Note_*: Do _not_ pre-filter your data prior to GC estimation, _except_
% possibly to improve stationarity (e.g notch-filtering to eliminate line noise
% or high-pass filtering to suppress low-frequency transients). Pre-filtering
% (of stationary data) may seriously degrade Granger-causal inference! If you
% want (time-domain) GC over a limited frequency range, rather calculate
% "band-limited" GC; to do this, calculate frequency-domain GCs over the full
% frequency range, then integrate over the desired frequency band [3]; see
% <smvgc_to_mvgc.html |smvgc_to_mvgc|>.
%
%% References
%
% [1] L. Barnett and A. K. Seth,
% <http://www.sciencedirect.com/science/article/pii/S0165027013003701 The MVGC
%     Multivariate Granger Causality Toolbox: A New Approach to Granger-causal
% Inference>, _J. Neurosci. Methods_ 223, 2014
% [ <matlab:open('mvgc_preprint.pdf') preprint> ].
%
% [2] A. B. Barrett, L. Barnett and A. K. Seth, "Multivariate Granger causality
% and generalized variance", _Phys. Rev. E_ 81(4), 2010.
%
% [3] L. Barnett and A. K. Seth, "Behaviour of Granger causality under
% filtering: Theoretical invariance and practical application", _J. Neurosci.
% Methods_ 201(2), 2011.
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%% Parameters

ntrials   = 1;     % number of trials
nobs      = 661;   % number of observations per trial

regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)

morder    = 'AIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
momax     = 10;     % maximum model order for model order estimation

acmaxlags = 661;   % maximum autocovariance lags (empty for automatic calculation)

tstat     = 'chi2';     % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
alpha     = 0.05;   % significance level for significance test
mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')

fs        = 200;    % sample rate (Hz)
fres      = [];     % frequency resolution (empty for automatic calculation)

seed      = 0;      % random seed (0 for unseeded)

%% Generate VAR test data (<mvgc_schema.html#3 |A3|>)
%
% _*Note:*_ This is where you would read in your own time series data; it should
% be assigned to the variable |X| (see below and <mvgchelp.html#4 Common
% variable names and data structures>).

% Seed random number generator.

rng_seed(seed);

% Get VAR coefficients for 5-node test network.

%AT = var5_test;
nvars = 10; %size(AT,1); % number of variables

% Residuals covariance matrix.

SIGT = eye(nvars);

% Generate multi-trial VAR time series data with normally distributed residuals
% for specified coefficients and covariance matrix.

% ptic('\n*** var_to_tsdata... ');
% X = var_to_tsdata(AT,SIGT,nobs,ntrials);
% ptoc;

%% Model order estimation (<mvgc_schema.html#3 |A2|>)

% Calculate information criteria up to specified maximum model order.

ptic('\n*** tsdata_to_infocrit\n');
[AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(X,momax,icregmode);
ptoc('*** tsdata_to_infocrit took ');

% Plot information criteria.

% figure(1); clf;
% plot_tsdata([AIC BIC]',{'AIC','BIC'},1/fs);
% title('Model order estimation');

%amo = size(AT,3); % actual model order

fprintf('\nbest model order (AIC) = %d\n',moAIC);
fprintf('best model order (BIC) = %d\n',moBIC);
%fprintf('actual model order     = %d\n',amo);

% Select model order.

if     strcmpi(morder,'actual')
 %   morder = amo;
  %  fprintf('\nusing actual model order = %d\n',morder);
elseif strcmpi(morder,'AIC')
    morder = moAIC;
    fprintf('\nusing AIC best model order = %d\n',morder);
elseif strcmpi(morder,'BIC')
    morder = moBIC;
    fprintf('\nusing BIC best model order = %d\n',morder);
else
    fprintf('\nusing specified model order = %d\n',morder);
end

%% VAR model estimation (<mvgc_schema.html#3 |A2|>)

% Estimate VAR model of selected order from data.

ptic('\n*** tsdata_to_var... ');
[A,SIG] = tsdata_to_var(X,morder,regmode);
ptoc;

% Check for failed regression

assert(~isbad(A),'VAR estimation failed');

% NOTE: at this point we have a model and are finished with the data! - all
% subsequent calculations work from the estimated VAR parameters A and SIG.

%% Autocovariance calculation (<mvgc_schema.html#3 |A5|>)

% The autocovariance sequence drives many Granger causality calculations (see
% next section). Now we calculate the autocovariance sequence G according to the
% VAR model, to as many lags as it takes to decay to below the numerical
% tolerance level, or to acmaxlags lags if specified (i.e. non-empty).

ptic('*** var_to_autocov... ');
[G,info] = var_to_autocov(A,SIG,acmaxlags);
ptoc;

% The above routine does a LOT of error checking and issues useful diagnostics.
% If there are problems with your data (e.g. non-stationarity, colinearity,
% etc.) there's a good chance it'll show up at this point - and the diagnostics
% may supply useful information as to what went wrong. It is thus essential to
% report and check for errors here.

var_info(info,true); % report results (and bail out on error)

%% Granger causality calculation: time domain  (<mvgc_schema.html#3 |A13|>)

% Calculate time-domain pairwise-conditional causalities - this just requires
% the autocovariance sequence.

ptic('*** autocov_to_pwcgc... ');
F = autocov_to_pwcgc(G);
ptoc;

% Check for failed GC calculation

assert(~isbad(F,false),'GC calculation failed');

% Significance test using theoretical null distribution, adjusting for multiple
% hypotheses.

pval = mvgc_pval(F,morder,nobs,ntrials,1,1,nvars-2,tstat); % take careful note of arguments!
sig  = significance(pval,alpha,mhtc);
%sig(isnan(sig))=0
subject_sig(:,:)=sig;
% Plot time-domain causal graph, p-values and significance.

figure(2); clf;
subplot(1,3,1);
plot_pw(F);
title('Pairwise-conditional GC');
subplot(1,3,2);
plot_pw(pval);
title('p-values');
subplot(1,3,3);
plot_pw(sig);
title(['Significant at p = ' num2str(alpha)]);
colormap('jet');colorbar();

% For good measure we calculate Seth's causal density (cd) measure - the mean
% pairwise-conditional causality. We don't have a theoretical sampling
% distribution for this.

cd = mean(F(~isnan(F)));

fprintf('\ncausal density = %f\n',cd);

%% Granger causality calculation: frequency domain  (<mvgc_schema.html#3 |A14|>)

% Calculate spectral pairwise-conditional causalities at given frequency
% resolution - again, this only requires the autocovariance sequence.

ptic('\n*** autocov_to_spwcgc... ');
f = autocov_to_spwcgc(G,fres);
ptoc;
% 
% % Check for failed spectral GC calculation
% 
assert(~isbad(f,false),'spectral GC calculation failed');
% 
% % Plot spectral causal graph.
% 
 figure(3); clf;
 plot_spw(f,fs);

%% Granger causality calculation: frequency domain -> time-domain  (<mvgc_schema.html#3 |A15|>)

% Check that spectral causalities average (integrate) to time-domain
% causalities, as they should according to theory.

fprintf('\nchecking that frequency-domain GC integrates to time-domain GC... \n');
Fint = smvgc_to_mvgc(f); % integrate spectral MVGCs
mad = maxabs(F-Fint);
madthreshold = 1e-5;
if mad < madthreshold
    fprintf('maximum absolute difference OK: = %.2e (< %.2e)\n',mad,madthreshold);
else
    fprintf(2,'WARNING: high maximum absolute difference = %e.2 (> %.2e)\n',mad,madthreshold);
end

%%
% <mvgc_demo.html back to top>
sig_sum=sig_sum+sig;

%end

