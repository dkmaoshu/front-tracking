	  �1  �   k820309              12.0        XP�V                                                                                                           
       D:\workf90_1\MAIN_ITF\main_euler\EULER\eulernew.f90 PHYSICAL_STATE              ADD SUBTRACT MULTIPLY DIVIDE_R FLUX_IN_X FLUX_IN_Y ALARM_IN_DENSITY TRIEM UF VF PF CF HF INTENG ASSIGN_SCALAR ASSIGN_EULER WEIGHTED_DIVISION_EULER INNER_PRODUCT                                                    
                            @                             
                                                                 |  #ASSIGN_SCALAR    #ASSIGN_EULER    #         @     @X                                              #A    #VALUE              D                                      (               #STATE              
                                      
      #         @     @X                                              #A    #W 	             D                                      (               #STATE              
                                  	     (              #STATE                                                               #IF_GREATER 
                                                              #IF_LESS                                                               #ADD    &         @    @X                                (                      #A    #B    #STATE              
                                       (              #STATE              
                                       (              #STATE                                                               #SUBTRACT    &         @    @X                                (                      #A    #B    #STATE              
                                       (              #STATE              
                                       (              #STATE                                                               #MULTIPLY    #INNER_PRODUCT    &         @    @X                                (                      #A    #B    #STATE              
                                      
                
                                       (              #STATE    %         @    @X                                                
       #A    #B              
                                       (              #STATE              
                                       (              #STATE                                                               #DIVIDE_R    &         @    @X                                 (                      #B    #A    #STATE              
                                       (              #STATE              
                                      
                                                             u #WEIGHTED_DIVISION_EULER    %         @    @X                                               
       #WEIGHTED_DIVISION_EULER%DMIN1    #YY    #MAIN_JUMP    #LEFT_DIFFERENCE    #RIGHT_DIFFERENCE     #NORMAL !   #LEFT_MAGNITUDE "   #RIGHT_MAGNITUDE #   #WAVE_NUMBER $                                                   DMIN1           
                                       (              #STATE              
                                       (              #STATE              
                                       (              #STATE              
                                        (              #STATE              
                                 !                   
    p          p            p                                    
                                  "     (              #STATE              
                                  #     (              #STATE              
                                  $                          �                                      u #DATA_FILL_P %                                                         u #ABSOLUTE_VALUE &   #         @                                 '                    #         @                                 (                   #RHO )   #U *   #V +   #P ,   #W -   #GAMMA .             
                                 )     
                
                                 *     
                
                                 +     
                
                                 ,     
                                                -                   
 	    p          p            p                                    
                                 .     
                                                  /     
                   
                       ���                          @                           0     '                    #FACING 1   #STATUS 2   #NUMBER 3                �                               1                                        �                               2                                       �                               3                                                              4     
                 
                 �������?        0.05D0                                            5     
                 
                       �?        0.25D0                                            6     
                 
                       �?        0.25D0                                             7                                                      4                                             8                                                      3                  @                                '(                    #VALUE 9   #GAMMA :                �                              9                              
  p          p            p                                       �                              :                
   %         @     X                            
                           #A ;   #VALUE <             
                                  ;     (              #STATE              
                                 <     
      %         @     X                                                       #A =   #VALUE >             
                                  =     (              #STATE              
                                 >     
      #         @      X                            %                   #UXY ?   #HIGHEST_ORDER @   #RANGE_SLIDE A             
D @   �                            ?            (             p ��������  & p ��������p            p                          #STATE              
                                  @                     
                                  A                                                                       B     DABS %         @     X                            &                   
       #ABSOLUTE_VALUE%DMAX1 C   #A D                                              C     DMAX1           
                                  D     (              #STATE    &         @                                E     (                      #U F   #STATE              
                                  F     (              #STATE    &         @                                G     (                      #U H   #STATE              
                                  H     (              #STATE    #         @                                  I                
   #WL J   #WR K   #NORMAL L   #NUM M   #X N   #Y O   #C1 P   #C2 Q   #VELOCITY R   #WAVE_TYPE S             
                                  J     (              #STATE              
                                  K     (              #STATE              
  @                              L                   
    p          p            p                                    
                                  M                     
                                 N     
                
                                 O     
                D @                               P     (               #STATE              D @                               Q     (               #STATE              D                                R     
                 D                                 S                            #         @                                  T                   #STATE_1 U   #STATE_2 V   #STATE_OR W   #WAVE_NUMBER X   #DIFFERENCE Y             
                                  U     (              #STATE              
                                  V     (              #STATE              
  @                               W     (              #STATE              
                                  X                     D @                               Y     (               #STATE    #         @                                 Z                   #U [   #NORMAL \   #UU ]             
                                  [     (              #STATE              
                                 \                   
    p          p            p                                    D                                 ]     (               #STATE    #         @                                  ^                  #COURANT_NUMBER%DMAX1 _   #A `   #X a   #Y b   #CH_SP c   #CRN d                                              _     DMAX1           
                                  `     (              #STATE              
                                 a     
                
                                 b     
                D                                 c     (               #STATE              D                                d     
       %         @                                 e                           #ST f             
                                  f     (              #STATE    %         @                                 g                    
       #A h             
                                  h     (              #STATE    %         @                                 i                    
       #A j             
                                  j     (              #STATE    %         @                                 k                    
       #A l             
                                  l     (              #STATE    %         @                                 m                    
       #A n             
                                  n     (              #STATE    %         @                                 o                    
       #A p             
                                  p     (              #STATE    %         @                                 q                    
       #A r             
                                  r     (              #STATE    %         @                                 s                    
       #A t             
                                  t     (              #STATE    #         @                                  u                   #UXY v   #HIGHEST_ORDER w             
D     �                            v            (             p ��������  & p ��������p            p                          #STATE              
@ @                               w           %         @                                 x                           #WAVE y             
                                  y           #         @                                 z                   #UU {   #NORMAL |   #U }             
                                  {     (              #STATE              
                                 |                   
    p          p            p                                    D                                 }     (               #STATE    %         @                                 ~                           #U              
                                       (              #STATE       �   K      fn#fn $   �   �   b   uapp(PHYSICAL_STATE     �  @   J   EULER_FUNCTIONS #   �  @   J   DATA_EXTRAPOLATION      e      i@|    �  Z      ASSIGN_SCALAR     �  S   a   ASSIGN_SCALAR%A $   .  @   a   ASSIGN_SCALAR%VALUE    n  V      ASSIGN_EULER    �  S   a   ASSIGN_EULER%A      S   a   ASSIGN_EULER%W    j  P      i@    �  M      i@      I      i@    P  i      ADD    �  S   a   ADD%A      S   a   ADD%B    _  N      i@    �  i      SUBTRACT      S   a   SUBTRACT%A    i  S   a   SUBTRACT%B    �  a      i@      i      MULTIPLY    �  @   a   MULTIPLY%A    �  S   a   MULTIPLY%B    	  ^      INNER_PRODUCT     w	  S   a   INNER_PRODUCT%A     �	  S   a   INNER_PRODUCT%B    
  N      i@    k
  i      DIVIDE_R    �
  S   a   DIVIDE_R%B    '  @   a   DIVIDE_R%A &   g  ]       gen@WEIGHTED_DIVISION (   �  �      WEIGHTED_DIVISION_EULER .   �  >      WEIGHTED_DIVISION_EULER%DMIN1 +   �  S   a   WEIGHTED_DIVISION_EULER%YY 2   P  S   a   WEIGHTED_DIVISION_EULER%MAIN_JUMP 8   �  S   a   WEIGHTED_DIVISION_EULER%LEFT_DIFFERENCE 9   �  S   a   WEIGHTED_DIVISION_EULER%RIGHT_DIFFERENCE /   I  �   a   WEIGHTED_DIVISION_EULER%NORMAL 7   �  S   a   WEIGHTED_DIVISION_EULER%LEFT_MAGNITUDE 8   0  S   a   WEIGHTED_DIVISION_EULER%RIGHT_MAGNITUDE 4   �  @   a   WEIGHTED_DIVISION_EULER%WAVE_NUMBER    �  Q       gen@DATA_FILL      T       gen@DABS #   h  H       ERROR_MESSAGE+GRID &   �  x       TRAN3+EULER_FUNCTIONS *   (  @   a   TRAN3%RHO+EULER_FUNCTIONS (   h  @   a   TRAN3%U+EULER_FUNCTIONS (   �  @   a   TRAN3%V+EULER_FUNCTIONS (   �  @   a   TRAN3%P+EULER_FUNCTIONS (   (  �   a   TRAN3%W+EULER_FUNCTIONS ,   �  @   a   TRAN3%GAMMA+EULER_FUNCTIONS     �  p       ERROR_DATA+GRID    l  t       CURVE_MAP !   �  P   a   CURVE_MAP%FACING !   0  P   a   CURVE_MAP%STATUS !   �  H   a   CURVE_MAP%NUMBER    �  v       MASS_DIFFUSION    >  v       VISCOSITY     �  v       HEAT_CONDUCTION    *  q       STATE_NUMBER    �  q       WAVE_NUMBER      f       STATE    r  �   a   STATE%VALUE      H   a   STATE%GAMMA    V  b       IF_GREATER    �  S   a   IF_GREATER%A !     @   a   IF_GREATER%VALUE    K  b       IF_LESS    �  S   a   IF_LESS%A       @   a   IF_LESS%VALUE    @  u       DATA_FILL_P     �  �   a   DATA_FILL_P%UXY *   d  @   a   DATA_FILL_P%HIGHEST_ORDER (   �  P   a   DATA_FILL_P%RANGE_SLIDE    �  =       DABS    1  q       ABSOLUTE_VALUE %   �  >      ABSOLUTE_VALUE%DMAX1 !   �  S   a   ABSOLUTE_VALUE%A    3  b       F    �  S   a   F%U    �  b       G    J  S   a   G%U    �  �       RIEMANN    E  S   a   RIEMANN%WL    �  S   a   RIEMANN%WR    �  �   a   RIEMANN%NORMAL       @   a   RIEMANN%NUM    �   @   a   RIEMANN%X    �   @   a   RIEMANN%Y    ?!  S   a   RIEMANN%C1    �!  S   a   RIEMANN%C2 !   �!  @   a   RIEMANN%VELOCITY "   %"  P   a   RIEMANN%WAVE_TYPE    u"  �       SIDE_CUT !   #  S   a   SIDE_CUT%STATE_1 !   Y#  S   a   SIDE_CUT%STATE_2 "   �#  S   a   SIDE_CUT%STATE_OR %   �#  @   a   SIDE_CUT%WAVE_NUMBER $   ?$  S   a   SIDE_CUT%DIFFERENCE $   �$  c       FIND_PHYSICAL_STATE &   �$  S   a   FIND_PHYSICAL_STATE%U +   H%  �   a   FIND_PHYSICAL_STATE%NORMAL '   �%  S   a   FIND_PHYSICAL_STATE%UU    /&  �       COURANT_NUMBER %   �&  >      COURANT_NUMBER%DMAX1 !   �&  S   a   COURANT_NUMBER%A !   K'  @   a   COURANT_NUMBER%X !   �'  @   a   COURANT_NUMBER%Y %   �'  S   a   COURANT_NUMBER%CH_SP #   (  @   a   COURANT_NUMBER%CRN    ^(  X       IF_PHYSICAL    �(  S   a   IF_PHYSICAL%ST    	)  W       DENSITY    `)  S   a   DENSITY%A    �)  W       X_VELOCITY    
*  S   a   X_VELOCITY%A    ]*  W       Y_VELOCITY    �*  S   a   Y_VELOCITY%A    +  W       PRESSURE    ^+  S   a   PRESSURE%A     �+  W       INTERNAL_ENERGY "   ,  S   a   INTERNAL_ENERGY%A    [,  W       ENTHALPY    �,  S   a   ENTHALPY%A    -  W       SOUND_SPEED    \-  S   a   SOUND_SPEED%A    �-  d       DATA_FILL_FLUX #   .  �   a   DATA_FILL_FLUX%UXY -   �.  @   a   DATA_FILL_FLUX%HIGHEST_ORDER    /  Z       WAVE_REVERS !   \/  @   a   WAVE_REVERS%WAVE (   �/  c       FIND_CONSERVATIVE_STATE +   �/  S   a   FIND_CONSERVATIVE_STATE%UU /   R0  �   a   FIND_CONSERVATIVE_STATE%NORMAL *   �0  S   a   FIND_CONSERVATIVE_STATE%U !   91  W       WHETHER_STATE_OK #   �1  S   a   WHETHER_STATE_OK%U 