	  �  +   k820309              12.0        ��Y                                                                                                           
       D:\workf90_1\MAIN_ITF\main_euler\PHYS_FLX\weno\flux_function_weno2.f90 FLUX_FUNCTIONS              EPWENO                                                    
                                                          
       %         @                                                   
       #FLUX_IN_X%DABS    #J    #W    #G                                                     DABS           
                                                       
                                                    
    p          p            p                                    
                                      
      %         @                                                    
       #W 	   #GAMMA 
             
                                 	                   
    p          p            p                                    
                                 
     
      %         @                                                    
       #W              
                                                    
    p          p            p                          %         @                                                    
       #W              
                                                    
    p          p            p                          %         @                                                   
       #FLUX_IN_Y%DABS    #J    #W    #G                                                     DABS           
                                                       
                                                    
    p          p            p                                    
                                      
                                                                     
      p          p            p                          #         @                                                    #WFLUX_X%DSQRT    #WU    #FLUX    #GAMMA                                                    DSQRT           
  @   �                                              
    p          & p ��������p          p            p          p                                    D @                                                 
     p          p            p                                    
  @                                   
      #         @                                                    #WFLUX_Y%DSQRT    #WU    #FLUX    #GAMMA                                                    DSQRT           
  @   �                                              
    p          & p ��������p          p            p          p                                    D @                                                 
     p          p            p                                    
  @                                   
      #         @                                                 
   #EVL     #EVR !   #EVLL "   #EVRR #   #F $   #WUU %   #FFF &   #WU '   #XY_DIR (   #FLUX )             
                                                     
    p          p          p            p          p                                    
                                 !                   
    p          p          p            p          p                                    
                                 "                   
     p          p          p            p          p                                    
                                 #                   
 !   p          p          p            p          p                                    
      �                           $                   
 "   p          & p ��������p          p            p          p                                    
      �                           %                   
 #   p          & p ��������p          p            p          p                                    
      �                           &                   
 $   p          & p ��������p          p            p          p                                    
      �                           '                   
    p          & p ��������p          p            p          p                                    
                                  (                     D                                )                   
     p          p            p                             �   ^      fn#fn $   �      b   uapp(FLUX_FUNCTIONS      @   J   GRID     U  @   J   EULER_FUNCTIONS *   �  y       FLUX_IN_X+EULER_FUNCTIONS /     =      FLUX_IN_X%DABS+EULER_FUNCTIONS ,   K  @   a   FLUX_IN_X%J+EULER_FUNCTIONS ,   �  �   a   FLUX_IN_X%W+EULER_FUNCTIONS ,     @   a   FLUX_IN_X%G+EULER_FUNCTIONS #   _  b       PF+EULER_FUNCTIONS %   �  �   a   PF%W+EULER_FUNCTIONS )   U  @   a   PF%GAMMA+EULER_FUNCTIONS #   �  W       UF+EULER_FUNCTIONS %   �  �   a   UF%W+EULER_FUNCTIONS #   �  W       VF+EULER_FUNCTIONS %   �  �   a   VF%W+EULER_FUNCTIONS *   k  y       FLUX_IN_Y+EULER_FUNCTIONS /   �  =      FLUX_IN_Y%DABS+EULER_FUNCTIONS ,   !  @   a   FLUX_IN_Y%J+EULER_FUNCTIONS ,   a  �   a   FLUX_IN_Y%W+EULER_FUNCTIONS ,   �  @   a   FLUX_IN_Y%G+EULER_FUNCTIONS    5  �       EMXY    �  x       WFLUX_X    A	  >      WFLUX_X%DSQRT    	  �   a   WFLUX_X%WU    C
  �   a   WFLUX_X%FLUX    �
  @   a   WFLUX_X%GAMMA      x       WFLUX_Y    �  >      WFLUX_Y%DSQRT    �  �   a   WFLUX_Y%WU    �  �   a   WFLUX_Y%FLUX    %  @   a   WFLUX_Y%GAMMA    e  �       WENO    
  �   a   WENO%EVL    �  �   a   WENO%EVR    r  �   a   WENO%EVLL    &  �   a   WENO%EVRR    �  �   a   WENO%F    �  �   a   WENO%WUU    b  �   a   WENO%FFF    &  �   a   WENO%WU    �  @   a   WENO%XY_DIR    *  �   a   WENO%FLUX 