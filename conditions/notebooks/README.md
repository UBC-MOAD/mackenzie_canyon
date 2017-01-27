The Jupyter Notebooks in this directory are used to explore and develop the boundary conditions for the Mackenzie Canyon model configuration.

The links below are to static renderings of the notebooks via
[nbviewer.jupyter.org](http://nbviewer.jupyter.org/).
Descriptions under the links below are from the first cell of the notebooks
(if that cell contains Markdown or raw text).

##Author: Idalia Machuca, University of British Columbia

* ##[prepare_stratification.ipynb](http://nbviewer.jupyter.org/urls/bitbucket.org/CanyonsUBC/mackenzie_canyon/raw/tip/conditions/prepare_stratification.ipynb)  
    
* ##[PlotANHA.ipynb](http://nbviewer.jupyter.org/urls/bitbucket.org/CanyonsUBC/mackenzie_canyon/raw/tip/conditions/PlotANHA.ipynb)  
    
    January, 2016: Getting familiarized with ANHA model output.  

* ##[room_of_requirement_2.ipynb](http://nbviewer.jupyter.org/urls/bitbucket.org/CanyonsUBC/mackenzie_canyon/raw/tip/conditions/room_of_requirement_2.ipynb)  
    
    **Notes: NEMO's analytical function**  
      
    **[Page 57]** In the vertical, the model mesh is determined by four things: (1) the bathymetry given in meters (2) the number of levels of the model (jpk) (3) the analytical transformation z(i, j, k) and the vertical scale factors (4) the masking system (the number of wet model levels at each (i, j) column of points).  
      
    **[Page 63]** In z-coordinate partial step (ln zps=.true.), the depths of the model levels are defined by the reference analytical function z0(k), except in the bottom layer. The thickness of the bottom layer is allowed to vary as a function of geographical location (lon, lat) to allow a better representation of the bathymetry.  
      
    **[Page 60]** The reference coordinate transformation z0(k) defines the arrays gdept0 and gdepw0 for t- and w-points, respectively. The vertical location of w- and t-levels is defined from the analytic expression of the depth z0(k) whose analytical derivative with respect to k provides the vertical scale factors.   
      
    It is possible to define a simple regular vertical grid by giving zero stretching (ppacr=0). It is often desirable to concentrate the vertical resolution near the ocean surface. The following function is proposed as a standard for  
    a z-coordinate (with either full or partial steps):  
      
    $z_0(k) = h_{sur} - h_0 k - h_1 log [ cosh ((k - h_{th}) / h_{cr})]$  
      
    With $h_{cr}$ and jpk, the four parameters (h...) have been determined such that these are satisfied, through an optimisation procedure using a bisection method.  
    $e_3(1 + 1/2) = 10$ and $e_3(jpk - 1/2) = 500$ and $z(1) = 0$ and $z(jpk) = -500$  
      
    Rather than entering the parameters (h...) directly, set ppsur=ppa0=ppa1=999999 in namcfg and specify the following:  
    * ppacr=hcr : stretching factor   
    * ppkth=hth : level at which maximum stretching  
    * ppdzmin : top layer minimum thickness  
    * pphmax : total depth of the ocean  
      
    **[Page 63]** With partial steps, layers from 1 to jpk-2 can have a thickness smaller than e3t(jk). The model deepest layer (jpk-1) is allowed to have either a smaller or larger thickness than e3t(jpk): the maximum thickness allowed is 2 * e3t(jpk -1)  
      
    Two variables in the namdom namelist are used to define the partial step vertical grid. The mimimum water thickness (in meters) allowed for a cell partially filled with bathymetry at level jk is the minimum of rn_e3zps_min (thickness in meters, usually 20 m) or e3t(jk) * rn_e3zps_rat (a fraction, usually 10%, of the default thickness e3t(jk)).  
      
      
    * phy_cst : initialization of ocean parameters and constants  
    * dom_nam  : domain initialization through namelist read (Namelist namdom : space & time domain)     
    * zgr_z   : Reference vertical z-coordinates  
    * dta_tsd_init : Temperature & Salinity data   

* ##[NEMOgcm_DataStructure.ipynb](http://nbviewer.jupyter.org/urls/bitbucket.org/CanyonsUBC/mackenzie_canyon/raw/tip/conditions/NEMOgcm_DataStructure.ipynb)  
    
    January, 2016: Getting familiarized with ANHA model output.  


##License

These notebooks and files are licensed under the Apache License, Version 2.0.
http://www.apache.org/licenses/LICENSE-2.0 by the CanyonsUBC group.
Please see the LICENSE file for details of the license.
