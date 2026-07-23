# Coda-Q-Inversion-USArray
MATLAB codes for coda-wave attenuation analysis based on the coda-Q inversion framework of Wang and Shearer (2019), with adaptations for large-scale USArray waveform processing and crustal attenuation imaging.

## Code Description
1. `process_coda_waveforms.m`  
   Processes SAC waveforms, calculates coda-wave energy envelopes, applies quality control, and generates inversion input files.

2. `plot_coda_trace.m`  
   Generates quality-control plots of waveform windows, coda envelopes, and fitting results.

3. `run_inv_coda_fixed_window.m`  
   Performs coda-Q inversion using a fixed coda-window length.

4. `run_inv_coda_fixed_tmax.m`  
   Performs coda-Q inversion using a fixed maximum lapse-time window.

5. `public_run_coda_inversion.m`  
   Controls the main inversion workflow and generates inversion results.

6. `reorder_coda_for_inversion.m`  
   Reorganizes event and station information for inversion.

7. `filter_event_station_pairs.m`  
   Filters event-station pairs according to data availability criteria.

8. `coda_q_inversion.m`  
   Core function for estimating coda attenuation parameters.

9. `hann_taper.m`  
   Applies a half-Hann taper to waveform data.

10. `readsac.m`  
   Reads SAC waveform files.

## Reference
Wang, W., & Shearer, P. M. (2019). An Improved Method to Determine Coda-Q, Earthquake Magnitude, and Site Amplification: Theory and Application to Southern California. Journal of Geophysical Research: Solid Earth, 124(1), 578–598. https://doi.org/10.1029/2018JB015961
