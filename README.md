# IW_Corr_Length (Date: 2021/02) 
 Calculating correlation length scales of internal waves based on a 3D matrix from the model results

 Creator: Yankun Gong (University of Western Australia)
 
 Email:  Yankun.Gong@research.uwa.edu.au

===========================================================================

1. Dataset in the folder

     (a) RowleyShoals.topo.3d:
	 A binary file of the bathymetry around the Rowley Shoals on the Australian NorthWest Shelf in the model domain.
	 			
     (b) Xgrid_meters.txt and Ygrid_meters.txt:
	 Model grid setup in x and y directions.
	   
     (c) Corr_ellipse_params.mat:
   An example of internal wave correlation ellipse from the MITgcm model. It includes semimajor, semiminor axes and the ellipse orientations.  
   
===========================================================================

2. Scripts and functions

     (a) get1_Var_correlation: (Mazloff et al. 2018)
         Based on the the 3D matrix of a variable from the model, we can find the correlation of wave properties between neighboring cells.
         This script is to calculateparameters of the correlation ellipse for internal waves, including crest widths, crest lengths and propagation directions.
	   
     (b) plot1_Var_correlation:
         Take the Rowley Shoals MITgcm model as an example, plot the internal wave properties via the correlation ellipse calculation. 
         (Gong et al. 2021).
		 
===========================================================================

3. Images in the folder
		 
     (a) Corr_ellipse_IWs.png:
         Correlation ellipses of internal waves around the Rowley Shoals within the entire model domain.
		 
     (b) Corr_lengthscale_IWs.png:
         A four-panel plot, characterizing internal waves properties (i.e. crest lengths, crest widths, and propagation directions) in a small domain.
	 				 
===========================================================================

4. References

    Mazloff, M.R., Cornuelle, B.D., Gille, S.T. and Verdy, A., 2018. Correlation lengths for estimating the large‐scale carbon and heat content of the Southern Ocean. Journal of Geophysical Research: Oceans, 123(2), pp.883-901.
  
    Gong Y., 2021. The Effects of Remotely-generated Internal Tides on Internal Waves on a Continental Shelf (Doctoral dissertation, The University of Western Australia).
