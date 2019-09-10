# SLAC D/P cross section ratio 

(The original data was presented as F2 ratio with the assumption Rd=Rp. So it is actually cross section ratio) 

## Data files: 
  * sigr (F2) d/p    : [xlsx](../data/JAM/10034.xlsx), [csv](../data/JAM/csv/10034.csv)   

## Reference:  
L.W.Whitlow, SLAC-Report-357, Ph.D. Thesis, Stanford University, March 1990.  

## Source: 
http://www.slac.stanford.edu/exp/e140/F2.DPRAT357


##HELP DOCUMENT

Source: http://www.slac.stanford.edu/exp/e140/HELP.DOC_357

File E.6:   F2.DPRATIO                                                           
                                                                                
File contains deuterium/hydrogen F_2 ratios extracted from the cross            
sections of Files E.2 and E.3, assuming that R¢d=R¢p (see Reference).           
Note, given this assumption, this file also presents deuterium/hydrogen         
cross section ratios.  Deuterium/hydrogen structure function ratios are         
presented wherever cross sections have been measured on both targets at         
identical or nearly identical kinematics.  Thus, the counting index, I,         
used here is NOT related the counting index used previously.                    
                                                                                
See documentation to Files E.4 and E.5 for a full description of the            
propagational properties of the errors dFSR, dFNM1, and dFNM2.                  
                                                                                
                                                                                
I       = running counting index provided for convenience and                   
          correspondence.                       

J       = experiment number, as ordered in Table 1.1 of Reference:              
          (1=E49a,2=E49b,3=E61,4=E87,5=E89a,6=E89b)                             
          Note that there are no hydrogen data from E139 and E140.   

x       = Bjorken scaling variable.

Q¢2     = 4-momentum transfer squared (GeV¢2).   

F2d/F2p = final extracted structure function ratio.  Units are [per             
          nucleon] for deuterium, and so, ratios are typically < 1. 

dST     = fractional statistical counting uncertainty in F_2 ratios. 

dSR     = fractional random systematic uncertainty in F_2 ratios. 

dN1     = fractional statistical uncertainty in F_2 ratios due to the           
          relative normalizations of the experiments.           

dN2     = fractional systematic uncertainty in F_2 due to the relative          
          normalizations of the experiments.    

stat    = total fractional random error, uncorrelated point to point.           
          stat=sqrt(dST¢2+dSR¢2).      
                                                   
syst    = total fractional experimental systematic error, correlated            
          between neighboring points.  syst=sqrt(dN1¢2+dN2¢2).                  
                                                                                
General remarks:                                                                
^^^^^^^^^^^^^^^^                                                                
stat and syst should be used primarily for making plots of the data.            
For sensitive analyses, the user should use the specific                        
(dST,dSR,dN1,dN2,dNM) error vector.  Note that errors of type dSY,              
dRC, dSZ become negligible in the structure function (or cross section)         
ratio, and are thus ignored.  Similarly, errors of type dSE are                 
identically zero in the ratio.  Additionally, there is an overall               
normalization uncertainty, denoted by dNM, of size                              
                                                                                
     dNM = .010 ,                                                               
                                                                                
also defined fractionally.                                                      
                                             
