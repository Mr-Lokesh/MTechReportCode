import math
import numpy as np

# Conversion of general orthotropic elastic constant into special orthotropy constants.
#
long_elastic_mod = float(input("Enter longitudinal elastic modulus:  "))
trans_elastic_mod = float(input("Enter transverse elastic modulus:  "))
shear_mod = float(input("Enter shear modulus:  "))
poisson_LT = float(input("Enter major poisson ratio:  "))
poisson_TL = float(input("Enter minor poisson ratio:  "))
theta = float(input("\nNote: Angle should be measure in anticlockwise direction from 'x' axis.\n"
                    "Enter angle of orientation of fiber:  "))


var1 = (((math.cos(math.radians(theta))**4)/long_elastic_mod) + ((math.sin(math.radians(theta))**4)/trans_elastic_mod))
var2 = 0.25 * ((1/shear_mod) - ((2*poisson_LT)/long_elastic_mod)) * (math.sin(math.radians(2*theta))**2)
long_mod_x = 1/ (var1 + var2)

var3 = ((math.cos(math.radians(theta))**4)/trans_elastic_mod) + ((math.sin(math.radians(theta))**4)/long_elastic_mod)
trans_mod_y = 1 / (var3 + var2)

var4 = (0.25 * ((1/long_elastic_mod)+((2*poisson_LT)/long_elastic_mod)+(1/trans_elastic_mod)-(1/shear_mod)))
var5 = (poisson_LT/long_elastic_mod) - (var4 * (math.sin(math.radians(2*theta))**2))
poisson_XY = long_mod_x * var5

var6 = (poisson_TL/trans_elastic_mod) - (var4 * (math.sin(math.radians(2*theta))**2))
poisson_YX = trans_mod_y * var6

var7 = ((1/long_elastic_mod)+(1/trans_elastic_mod)+((2*poisson_LT)/long_elastic_mod))
var8 = ((1/long_elastic_mod)+((2*poisson_LT)/long_elastic_mod)+(1/trans_elastic_mod)-(1/shear_mod))
shear_mod_XY = (1/(var7 - (var8 * (math.cos(math.radians(2*theta))**2))))

var9 = (poisson_LT + (long_elastic_mod / trans_elastic_mod) - (long_elastic_mod / (2*shear_mod)))
var10 = (1 + (long_elastic_mod/trans_elastic_mod) + (2*poisson_LT) - (long_elastic_mod/shear_mod))
coff_Mx = math.sin(math.radians(2*theta)) * (var9 - (var10 * (math.cos(math.radians(theta))**2)))

coff_My = math.sin(math.radians(2*theta)) * (var9 - (var10 * (math.sin(math.radians(theta))**2)))

print("\n------------------------------------------ RESULTS --------------------------------------")
print("\nLongitudinal Modulus in 'X' direction:  ", long_mod_x)
print("Modulus in 'Y' direction:   ", trans_mod_y)
print("Shear Modulus (XY):   ", shear_mod_XY)
print("Major Poisson Ratio (XY):   ", poisson_XY)
print("Minor Poisson's Ratio (YX):   ", poisson_YX)
print("Cross Coefficient value for 'X' direction:   ", coff_Mx)
print("Cross Coefficient value for 'Y' direction:   ", coff_My)
print("\n-----------------------------------------------------------------------------------------")
#
print("\n==================================== Strain Calculation ===================================\n")
sigma_x = float(input("Enter normal stress value in 'X' direction:   "))
sigma_y = float(input("Enter normal stress value in 'Y' direction:   "))
tau_xy = float(input("Enter shear stress value in 'XY' direction:   "))

if (theta == 0):     # CASE OF SPECIAL ORTHOTROPY
    strain_L = (sigma_x / long_elastic_mod) - (poisson_TL * (sigma_y / trans_elastic_mod))
    strain_T = (sigma_y / trans_elastic_mod) - (poisson_LT * (sigma_x / long_elastic_mod))
    shear_strain_LT = (tau_xy / shear_mod)

    print("\n------------------------------------- Strain Results -----------------------------------")
    print("\nValue of strain in 'X' direction:  ", strain_L)
    print("Value of strain in 'Y' direction:  ", strain_T)
    print("Value of shear strain in 'XY' direction:  ", shear_strain_LT)
    print("\n --------------------------------------------------------------------------------------")

else:
    strain_X = (sigma_x/long_mod_x) - ( poisson_YX*(sigma_y/trans_mod_y) ) - ( coff_Mx*(tau_xy/long_elastic_mod) )
    strain_Y = (sigma_y/trans_mod_y) - ( poisson_XY*(sigma_x/long_mod_x) ) - ( coff_My*(tau_xy/long_elastic_mod) )
    shear_strain_XY = (tau_xy/shear_mod_XY) - (coff_Mx*(sigma_x/long_elastic_mod)) - (coff_My*(sigma_y/long_elastic_mod))

    print("\n------------------------------------- Strain Results -----------------------------------")
    print("\nValue of strain in 'X' direction:  ", strain_X)
    print("Value of strain in 'Y' direction:  ", strain_Y)
    print("Value of shear strain in 'XY' direction:  ", shear_strain_XY)
    print("\n --------------------------------------------------------------------------------------")


"""------------------------- CALCULATION OF STIFFNESS MATRIX [Q] FOR L-T PLANE---------------------------"""

long_elastic_mod = float(input("Enter longitudinal elastic modulus:  "))
trans_elastic_mod = float(input("Enter transverse elastic modulus:  "))
shear_mod = float(input("Enter shear modulus:  "))
poisson_LT = float(input("Enter major poisson ratio:  "))
poisson_TL = float(input("Enter minor poisson ratio:  "))

Q11 = long_elastic_mod / (1 - (poisson_LT * poisson_TL))
Q22 = trans_elastic_mod / (1 - (poisson_LT * poisson_TL))
Q12 = (poisson_LT * trans_elastic_mod) / (1 - (poisson_LT * poisson_TL))
Q66 = shear_mod

Q_matrix = np.array([
    [ Q11 , Q12 , 0],
    [ Q12 , Q22 , 0],
    [ 0 , 0 , Q66] ])

print("\n [Q] Matrix:")
print("-----------------------------------------------")
print(Q_matrix)
print("-----------------------------------------------\n")

"""-------------------------- CALCULATION OF COMPLIANCE MATRIX [S] IN L-T PLANE--------------------"""

S11 = 1 / long_elastic_mod
S22 = 1 / trans_elastic_mod
S12 = (-1) * (poisson_LT / long_elastic_mod)
S66 = 1 / shear_mod

S_matrix = np.array([
    [ S11 , S12 , 0],
    [ S12 , S22 , 0],
    [ 0 , 0 , S66] ])

print("\n [S] Matrix:")
print("-----------------------------------------------")
print(S_matrix)
print("-----------------------------------------------\n")

s = np.linalg.inv(Q_matrix)
print(s)
#
# """------------------------- CALCULATION OF STIFFNESS MATRIX [Q BAR] FOR X-Y PLANE----------------------"""

theta = float(input("\nNote: Angle should be measure in anticlockwise direction from 'x' axis.\n"
                    "Enter angle of orientation of fiber:  "))
sine = (math.sin(math.radians(theta)))
cosine = (math.cos(math.radians(theta)))

print(cosine)


Q11_bar = Q11*(cosine**4) + Q22*(sine**4) + 2*(Q12+(2*Q66))*(sine**2)*(cosine**2)
Q22_bar = Q11*(sine**4) + Q22*(cosine**4) + 2*(Q12+(2*Q66))*(sine**2)*(cosine**2)
Q12_bar = (Q11+Q22-(4*Q66))*(sine**2)*(cosine**2) + Q12*((cosine**4)+(sine**4))
Q66_bar = (Q11+Q22-(2*Q12)-(2*Q66))*(sine**2)*(cosine**2) + Q66*((cosine**4)+(sine**4))
Q16_bar = (Q11+Q12-(2*Q66))*(cosine**3)*(sine) - (Q22-Q12-(2*Q66))*(sine**3)*(cosine)
Q26_bar = (Q11+Q12-(2*Q66))*(cosine)*(sine**3) - (Q22-Q12-(2*Q66))*(sine)*(cosine**3)

Q_bar_matrix = np.array([
    [Q11_bar , Q12_bar , Q16_bar],
    [Q12_bar , Q22_bar , Q26_bar],
    [Q16_bar , Q26_bar , Q66_bar] ])
#
print("\n [Q bar]  Matrix:")
print("-----------------------------------------------")
print(Q_bar_matrix)
print("-----------------------------------------------\n")



"""-------------------------- CALCULATION OF COMPLIANCE MATRIX [S bar] IN x-y PLANE--------------------"""

S11_bar = S11*(cosine**4) + S22*(sine**4) + ((2*S12)+S66)*(sine**2)*(cosine**2)
S22_bar = S11*(sine**4) + S22*(cosine**4) + ((2*S12)+S66)*(sine**2)*(cosine**2)
S12_bar = (S11+S22-S66)*(sine**2)*(cosine**2) + S12*((cosine**4)+(sine**4))
S66_bar = 2*((2*S11)+(2*S22)-(4*S12)-S66)*(sine**2)*(cosine**2) + S66*((cosine**4)+(sine**4))
S16_bar = ((2*S11)-(2*S12)-S66)*(cosine**3)*(sine) - ((2*S22-(2*S12)-S66))*(sine**3)*(cosine)
S26_bar = ((2*S11)-(2*S12)-S66)*(cosine)*(sine**3) - ((2*S22-(2*S12)-S66))*(sine)*(cosine**3)

S_bar_matrix = np.array([
    [S11_bar , S12_bar , S16_bar],
    [S12_bar , S22_bar , S26_bar],
    [S16_bar , S26_bar , S66_bar] ])

print("\n [S bar]  Matrix:")
print("-----------------------------------------------")
print(S_bar_matrix)
print("-----------------------------------------------\n")

S_bar_matrix = np.linalg.inv(Q_bar_matrix)

print("\n [S bar]  Matrix:")
print("-----------------------------------------------")
print(S_bar_matrix)
print("-----------------------------------------------\n")

"""------------------------------Calculation for stress and strains----------------------------"""

sigma_x = float(input("Enter normal stress value in 'X' direction:   "))
sigma_y = float(input("Enter normal stress value in 'Y' direction:   "))
tau_xy = float(input("Enter shear stress value in 'XY' direction:   "))

stress_matrix = np.array([
    [sigma_x],[sigma_y],[tau_xy] ])
print(stress_matrix)
# print(stress_matrix.shape)

strain_matrix = S_bar_matrix .dot(stress_matrix)
print(strain_matrix)
print("\n Strain Matrix :")
print("-----------------------------------------------")
print(strain_matrix)
print("-----------------------------------------------\n")



