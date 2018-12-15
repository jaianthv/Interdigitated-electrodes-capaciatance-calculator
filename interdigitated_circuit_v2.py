import numpy as np
import sys
import os
import math



#No of tooth

N = 200;

#er of infinite layer

#er_infi = xx;


# lenth of electrodes

L = 2 * 10**-3;

#width of the tooth
w = 10 * 10**-6;

#Distance between opposite tooth +v and - V

G = 5 * 10**-6;

#Distance between two electrodes

lamb = 2 * (w+G) ; #5 * 10**-6;

#thickness of dielectric

h = 50 * 10**-9;


# r = h/lamb

r = h/lamb;

# eta = 2w/lamb

eta = 2*w/lamb;

# Permitivitty of vaccum

e0 = 8.85 * 10**(-12);

# Permitivitty of the dielectric

er = 8.9 ;

#permitivitty of substrate
es = 3.8;


######################################################################################

# Interior capacitance

# The value of q = exp(-4.pi.r)

q = np.exp(-4*math.pi*r)

# To obtain the value of 'k' with Jacobi elliptic integral of first kind v2(0,q)

v2 = (2* math.pow(q,(1/4))) + (2 * math.pow(q,(9/4))) + (2 * math.pow(q,(25/4)));

# To obtain the value of 'k' with Jacobi elliptic integral of first kind v3(0,q)

v3 = 1 + 2*q + 2*q**4 + 2*q**9 ;

# Finally K is

k = (v2/v3)**2;


# t4 is

t4 = 1/k;

# Last one t2 = sn(K(k)eta, k) This is a long one :(

## K(k)

kk = (math.pi/2)*(1 + ((1/2)**2 * k**2) + ( ( (1*3)/(2*4) )**2 * k**4 )) * eta;

## Now sn(u,k) --> u = kk(with eta) and k is still k

u = kk;

sn =  u - ( (1/6) * (1+k**2) * u**3 ) + ( (1/120) * (1 + math.pow(14,math.pow(k,2)) + k**4 ) * u**5 ) ;

t2 = sn ;

# ki

ki = t2 * math.sqrt(( ( (t4)**2 - 1 )/( (t4)**2 - (t2)**2) ));

# Modulus of ki

kim = math.sqrt(1 - (ki**2));

# to find K(Ki) and K(kim)

kk_ki = (math.pi/2)*(1 + ((1/2)**2 * ki**2) + ( ( (1*3)/(2*4) )**2 * ki**4 )) ;

kk_kim = (math.pi/2)*(1 + ((1/2)**2 * kim**2) + ( ( (1*3)/(2*4) )**2 * kim**4 )) ;

# internal capacitance is

Cint = e0 * er * L * ( kk_ki / kk_kim );

##############################################################################

# Exterior capacitance

x_t3= (math.pi*(1-eta)/8*r);
x_t4= (math.pi*(eta+1)/8*r);

t3 = 1+ (x_t3**2/2) + (x_t3**4/(4*3*2*1)) ;   #math.acosh((x_t3));
t4 = 1+ (x_t4**2/2) + (x_t4**4/(4*3*2*1)) ;   #math.acosh((x_t4));

# ke now

ke = 1/t3 * math.sqrt((t4**2 - t3**2)/(t4**2 - 1));

# ke modulus

kem = math.sqrt(1-ke**2)

# to find K(ke) and K(kem)

kk_ke = (math.pi/2)*(1 + ((1/2)**2 * ke**2) + ( ( (1*3)/(2*4) )**2 * ke**4 )) ;

kk_kem = (math.pi/2)*(1 + ((1/2)**2 * kem**2) + ( ( (1*3)/(2*4) )**2 * kem**4 )) ;

# external capacitance is

Cext = e0 * er * L * (kk_ke/kk_kem);


###############################################################################

# Infinite layer  capacitance


## internal electrodes


#Kinfinity

k_inf = math.sin((math.pi/2)* eta);

k_infm = math.sqrt(1-(k_inf)**2);


# K(kinfi) and K(kinfm)


kk_inf = (math.pi/2)*(1 + ((1/2)**2 * k_inf**2) + ( ( (1*3)/(2*4) )**2 * k_inf**4 )) ;

kk_infm = (math.pi/2)*(1 + ((1/2)**2 * k_infm**2) + ( ( (1*3)/(2*4) )**2 * k_infm**4 )) ;


# Internal capacitance with infinite dielectric thickness


#C_I_infi = e0* er_infi *L * (kk_infi/kk_infim);

C_I_infi = e0* er *L * (kk_inf/kk_infm);


## external electrodes


#K infinity  
ke_inf = ((2 * math.sqrt(eta)) / (1+ eta)) ;

q = 1-(ke_inf)**2;

ke_infm = math.sqrt(q);

# K(k infi) and K (k infim)


kke_inf = (math.pi/2)*(1 + ((1/2)**2 * ke_inf**2) + ( ( (1*3)/(2*4) )**2 * ke_inf**4 )) ;

kke_infm = (math.pi/2)*(1 + ((1/2)**2 * ke_infm**2) + ( ( (1*3)/(2*4) )**2 * ke_infm**4 )) ;

# External capacitance with infite dielectric thickness

#C_E_infi = e0 * er_infi * L (kke_inf/kke_infm);

C_E_infi = e0 * er * L * (kke_inf/kke_infm);
##### Total intermal capacitance


Term1 = kk_inf/kk_infm;

Term2 = (er - 1) * kk_ki / kk_kim;

Term3 = es * kk_ki / kk_kim;

C_Internal = e0 * L * Term1 * Term2 * Term3;


##### Total external capacitance


Term1_e = kke_inf/kke_infm;

Term2_e = (er - 1) * kke_inf / kke_infm;

Term3_e = es * kke_inf / kke_infm;

C_external = e0 * L * Term1_e * Term2_e * Term3_e;


#####Final capacitance For embedded electrodes



Final = (N-3) * (C_Internal/2)+ 2 * ((C_Internal *C_external)/(C_Internal  + C_external));




######## For electrodes on top


##### Total intermal capacitance


Term1_eletop = kk_inf/kk_infm;



Term3_eletop = es * kk_ki / kk_kim;

C_Internal_eletop = e0 * L * Term1_eletop  * Term3_eletop;


##### Total external capacitance


Term1__e_eletop = kke_inf/kke_infm;



Term3_e_eletop = es * kke_inf / kke_infm;

C_external_eletop = e0 * L * Term1__e_eletop  * Term3_e_eletop;

#####Final capacitance For NOT embedded electrodes



Final_NOTembedded = (N-3) * (C_Internal_eletop/2)+ 2 * ((C_Internal_eletop *C_external_eletop)/(C_Internal_eletop  + C_external_eletop));




print ("Electrodes embedded");
print (Final);
print ("Electrodes not embedded");
print  (Final_NOTembedded);














