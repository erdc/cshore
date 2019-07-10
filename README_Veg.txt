1. The CSHORE source code remains the same when IVEG=0-2. The vegetation advances developed by Dr. Qin J. Chen's research team are activated only when IVEG=3. A new makeinfile matlab script has been added: /usace_distribute_bundle/mfiles/makeinfile_usace_vegfeature.m


2. When IVEG = 3, new parameters (NFR) and new variables (IDISS, IFV,FREQMIN, FREQMAX, FREQNUM, FREQMINBC, VEGCDM, TZ, VMEASOMEG, VMEASSE, VMEASWNUM, NMEASSPEC, GAMJONSWAP) are added.
   NFR       = maximum number of frequency beams for JONSWAP spectrum
   IDISS     = controls energy dissipation model (due to vegetation) 
   IFV       = controls phase-averaged depth-integrated drag model
   FREQMIN   = the minimum cutoff frequency 
   FREQMAX   = the maximum cutoff frequency
   FREQNUM   = the number of frequency components 
   FREQMINBC = the cutoff frequency at the offshore BC 
   GAMJONSWAP= the gamma parameter in JONSWAP spectrum
   VEGCDM    = the second set of drag coefficient if IFV = 2
   VMEASOMEG = the frequencies of the measured wave spectrum
   VMEASSE   = the spectral density of the measured wave spectrum
   VMEASWNUM = the wave numbers corresponding to the frequencies (VMEASOMEG) 
   NMEASSPEC = the number of frequency components in the measured spectrum


3. When IVEG = 3, new subroutines are added:
   subroutine DVEG:  compute the energy dissipation due to vegetation (DVEGSTA). Three different models are used: 
                     IDISS=1: Mendez and Losada (2004)  
                     IDISS=2: Chen and Zhao (2012) with JONSWAP spectrum
                     IDISS=3: Chen and Zhao (2012) with measured spectrum    

   subroutine INTERP1: doing interpolation.
   subroutine FINDHV2HTOMTABLE: compute the parameter (HV2HTOM) used in subroutine PHASEAVEFV().
   subroutine DISPERSION: compute wavenumber with given water depth and wave period.
   subroutine PHASEAVEFV: compute the phase-averaged depth-integrated drag (STREAMSTRESSSTA) due to vegetation. Three different models are used:
                     IFV = 1: original CSHORE 
                     IFV = 2: parametric model. This model follows the parametric model in Zhu and Chen (2019) for both submerged and emergent vegetation. For waves with large nonlinearity, a different Cdm value may be needed to compensate for the overestimated mean current (Stokes drift).
                     IFV = 3: hybrid model. The submerged part is computed with parametric model and the canopy part is computed with the formula based on linear wave theory (Dean and Bender 2006, a regular wave version). About determining the submerged part in random waves, we found that 95% of the water depth works for all 7 USDA flume experiments. No second set of Cdm is needed. 


4. If IDISS = 3 (use energy dissipation model by Chen and Zhao (2012) together with measured wave spectrum), the measured wave spectrum is read in the main code and stored in VMEASOMEG, VMEASSE and VMEASWNUM. The format of the measured spectrum file is: three columns representing frequency, spectral density, wave number, respectively.

5. Two new models of the phase-averaged depth-integrated drag (STREAMSTRESSSTA) are implemented in the cross-shore momentum equation. The lines for calculating "WSETUP" are revised.

6. Four application cases are added in /usace_distribute_bundle/applications

