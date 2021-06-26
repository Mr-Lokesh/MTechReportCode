
import math
import numpy as np
import pandas as pd
# -----------------------------------Input Parameters --------------------------------------------
vol_frac_fiber = float(input("Enter volume fraction of fiber:  "))
young_mod_fiber = float(input("Enter Young's Modulus (in Mpa) of fiber:  "))
young_mod_matrix = float(input("Enter Young's Modulus (in Mpa) of Matrix:  "))
shear_mod_fiber = float(input("Enter Shear Modulus (in Mpa) of fiber:  "))
shear_mod_matrix = float(input("Enter Shear Modulus (in Mpa) of Matrix:  "))
poisson_ratio_fiber = float(input("Enter poisson's ratio of fiber:  "))
poisson_ratio_matrix = float(input("Enter poisson's ratio of matrix:  "))
# ----------------------------------------------------------------------------------------------



#  Longitudinal Modulus
uni_longi_mod = young_mod_matrix * (vol_frac_fiber*((young_mod_fiber/young_mod_matrix)-1)+1)
print("\nLongitudinal modulus of composite: ",uni_longi_mod)

# Transverse Modulus
ratio = (young_mod_fiber / young_mod_matrix)
eta = (ratio - 1) / (ratio + 2)
trans_mod = young_mod_matrix * (1 + (2*eta*vol_frac_fiber)) / (1 - (eta*vol_frac_fiber))
print("Transverse Modulus value is: ",trans_mod)

# Shear Modulus
shear_mod_ratio = (shear_mod_fiber / shear_mod_matrix)
eta = (shear_mod_ratio - 1) / (shear_mod_ratio + 1)
shear_mod = shear_mod_matrix * (1 + (1 * eta * vol_frac_fiber)) / (1 - (eta * vol_frac_fiber))
print("Shear Modulus value (in Mpa) is: ", shear_mod)

# Poisson's Ratio
poisson_LT = (poisson_ratio_fiber * vol_frac_fiber)+(poisson_ratio_matrix*(1-vol_frac_fiber))
print("Major Poisson's Ratio is:   ", poisson_LT)

uni_longi_mod = young_mod_matrix * (vol_frac_fiber * ((young_mod_fiber / young_mod_matrix) - 1) + 1)
poisson_LT = (poisson_ratio_fiber * vol_frac_fiber) + (poisson_ratio_matrix * (1 - vol_frac_fiber))
poisson_TL = (young_mod_matrix / uni_longi_mod) * poisson_LT
print("Minor Poisson's Ratio is:   ",poisson_TL)


# Laminates Stress Distribution and properties

A11,A12,A22,A16,A26,A66,B11,B12,B16,B22,B26,B66,D11,D12,D16,D22,D26,D66 = 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
sum = 0
var = 0

orientation = []
thickness = []
Z_elements_from_top = []
list_stress_values_X = []
list_stress_values_Y = []
list_stress_values_XY = []
list_stress_values_L = []
list_stress_values_T = []
list_stress_values_LT = []
list_location = []

list_stress_values_top = []
list_stress_values_bottom = []

#-----------------------------------------------------Input Layers ----------------------------------------
print("\n")
layers = int(input("Enter number of layers in laminate:   "))
n = layers

for i in range(layers):
    height = float(input("Enter thickness of layer"+str([i+1])+" from top:   "))
    thickness.append(height)
    angle = float(input("Enter orientation angle of layer"+str([i+1])+" in degrees:   "))
    orientation.append(angle)
# ----------------------------------------------------------------------------------------------------------
# Calculation of Z values for each layer according to thickness
for j in range(0, len(thickness)):
    sum = sum + thickness[j]
median = sum / 2

for j in range(0, len(thickness)):
    if j == 0:
        Z_elements_from_top.append(-median)

    var = var + thickness[j]
    z_value = var - median
    Z_elements_from_top.append(z_value)

print("\nOrientation angle from top to bottom: ",orientation)
print("Thickness of each layer from top to bottom: ",thickness)
#print("Z elements required for [B] and [D] matrix:  ",Z_elements_from_top)

for i in range(layers):
    # uni_longi_mod = float(input("\nEnter longitudinal elastic modulus of layer"+str([i+1])+":  "))
    # trans_mod = float(input("Enter transverse elastic modulus of layer"+str([i+1])+":  "))
    # shear_mod = float(input("Enter shear modulus of layer"+str([i+1])+":  "))
    # poisson_LT = float(input("Enter major poisson ratio of layer"+str([i+1])+":  "))
    # poisson_TL = float(input("Enter minor poisson ratio of layer"+str([i+1])+":  "))

    theta = orientation[i]
    sine = (math.sin(math.radians(theta)))
    cosine = (math.cos(math.radians(theta)))

    Q11 = uni_longi_mod / (1 - (poisson_LT * poisson_TL))
    Q22 = trans_mod / (1 - (poisson_LT * poisson_TL))
    Q12 = (poisson_LT * trans_mod) / (1 - (poisson_LT * poisson_TL))
    Q66 = shear_mod

    Q11_bar = Q11 * (cosine ** 4) + Q22 * (sine ** 4) + 2 * (Q12 + (2 * Q66)) * (sine ** 2) * (cosine ** 2)
    Q22_bar = Q11 * (sine ** 4) + Q22 * (cosine ** 4) + 2 * (Q12 + (2 * Q66)) * (sine ** 2) * (cosine ** 2)
    Q12_bar = (Q11 + Q22 - (4 * Q66)) * (sine ** 2) * (cosine ** 2) + Q12 * ((cosine ** 4) + (sine ** 4))
    Q66_bar = (Q11 + Q22 - (2 * Q12) - (2 * Q66)) * (sine ** 2) * (cosine ** 2) + Q66 * ((cosine ** 4) + (sine ** 4))
    Q16_bar = (Q11 + Q12 - (2 * Q66)) * (cosine ** 3) * (sine) - (Q22 - Q12 - (2 * Q66)) * (sine ** 3) * (cosine)
    Q26_bar = (Q11 + Q12 - (2 * Q66)) * (cosine) * (sine ** 3) - (Q22 - Q12 - (2 * Q66)) * (sine) * (cosine ** 3)

    A11 = A11 + thickness[i] * float(Q11_bar)
    A22 = A22 + thickness[i] * float(Q22_bar)
    A12 = A12 + thickness[i] * float(Q12_bar)
    A66 = A66 + thickness[i] * float(Q66_bar)
    A16 = A16 + thickness[i] * float(Q16_bar)
    A26 = A26 + thickness[i] * float(Q26_bar)

    B11 = B11 + (((Z_elements_from_top[i + 1] ** 2) - (Z_elements_from_top[i] ** 2)) / 2) * float(Q11_bar)
    B22 = B22 + (((Z_elements_from_top[i + 1] ** 2) - (Z_elements_from_top[i] ** 2)) / 2) * float(Q22_bar)
    B12 = B12 + (((Z_elements_from_top[i + 1] ** 2) - (Z_elements_from_top[i] ** 2)) / 2) * float(Q12_bar)
    B66 = B66 + (((Z_elements_from_top[i + 1] ** 2) - (Z_elements_from_top[i] ** 2)) / 2) * float(Q66_bar)
    B16 = B16 + (((Z_elements_from_top[i + 1] ** 2) - (Z_elements_from_top[i] ** 2)) / 2) * float(Q16_bar)
    B26 = B26 + (((Z_elements_from_top[i + 1] ** 2) - (Z_elements_from_top[i] ** 2)) / 2) * float(Q26_bar)

    D11 = D11 + (((Z_elements_from_top[i + 1] ** 3) - (Z_elements_from_top[i] ** 3)) / 3) * float(Q11_bar)
    D22 = D22 + (((Z_elements_from_top[i + 1] ** 3) - (Z_elements_from_top[i] ** 3)) / 3) * float(Q22_bar)
    D12 = D12 + (((Z_elements_from_top[i + 1] ** 3) - (Z_elements_from_top[i] ** 3)) / 3) * float(Q12_bar)
    D66 = D66 + (((Z_elements_from_top[i + 1] ** 3) - (Z_elements_from_top[i] ** 3)) / 3) * float(Q66_bar)
    D16 = D16 + (((Z_elements_from_top[i + 1] ** 3) - (Z_elements_from_top[i] ** 3)) / 3) * float(Q16_bar)
    D26 = D26 + (((Z_elements_from_top[i + 1] ** 3) - (Z_elements_from_top[i] ** 3)) / 3) * float(Q26_bar)

A_matrix = np.matrix([
    [A11, A12, A16],
    [A12, A22, A26],
    [A16, A26, A66]])

B_matrix = np.matrix([
    [B11, B12, B16],
    [B12, B22, B26],
    [B16, B26, B66]])

D_matrix = np.matrix([
    [D11, D12, D16],
    [D12, D22, D26],
    [D16, D26, D66]])

print("\n [A]  Matrix:")
print("-------------------------------------------------------------------")
print(A_matrix)
print("-----------------------------------------------------------------\n")

print("\n [B]  Matrix:")
print("-------------------------------------------------------------------")
print(B_matrix)
print("-----------------------------------------------------------------\n")

print("\n [D]  Matrix:")
print("-------------------------------------------------------------------")
print(D_matrix)
print("-----------------------------------------------------------------\n")


Stiffness_matrix = np.matrix([
        [A11 , A12 , A16 , B11 , B12 , B16],
        [A12 , A22 , A26 , B12 , B22 , B26],
        [A16 , A26 , A66 , B16 , B26 , B66],
        [B11 , B12 , B16 , D11 , D12 , D16],
        [B12 , B22 , B26 , D12 , D22 , D26],
        [B16 , B26 , B66 , D16 , D26 , D66]
    ])
Inverse_stiffness_matrix = Stiffness_matrix.getI()
Inverse_D_matrix = D_matrix.getI()
Inverse_A_matrix = A_matrix.getI()
#print(Stiffness_matrix)
tensile_modulus_E11 = 1/ (sum * Inverse_A_matrix[0,0])
shear_modulus_G12 = 1/ (sum * Inverse_stiffness_matrix[2,2])
flexural_modulus = 12 / ((sum**3) * Inverse_D_matrix[0,0])
print("Tensile Modulus (E11):  ", tensile_modulus_E11)
print("Shear Modulus (G12):  ", shear_modulus_G12)
print("Flexural Modulus:  ", flexural_modulus)


# Stress Distribution Calculation inside each layer of laminate

N_x = float(input("\nEnter stress in 'X' direction:   "))
N_y = float(input("Enter stress in 'Y' direction:   "))
N_xy = float(input("Enter shear stress in 'XY' direction:   "))
M_x = float(input("Enter moment in 'X' direction:   "))
M_y = float(input("Enter moment in 'Y' direction:   "))
M_xy = float(input("Enter moment in 'XY' direction:   "))

NM_matrix = np.matrix([
        [N_x],[N_y],[N_xy],[M_x],[M_y],[M_xy]
    ])

Strain_curvature_matrix = Inverse_stiffness_matrix * NM_matrix

laminate_strains = Strain_curvature_matrix[0:3,0:]
laminate_curvature = Strain_curvature_matrix[3:,0:]

for i in range(0, len(Z_elements_from_top)):
    #print("\nStrains in " + str([i + 1]) + " layer from top:")
    lamina_strain = laminate_strains + Z_elements_from_top[i] * laminate_curvature

for i in range(layers):
    # long_elastic_mod = float(input("\nEnter longitudinal elastic modulus of layer"+str([i+1])+":  "))
    # trans_elastic_mod = float(input("Enter transverse elastic modulus of layer"+str([i+1])+":  "))
    # shear_mod = float(input("Enter shear modulus of layer"+str([i+1])+":  "))
    # poisson_LT = float(input("Enter major poisson ratio of layer"+str([i+1])+":  "))
    # poisson_TL = float(input("Enter minor poisson ratio of layer"+str([i+1])+":  "))
    theta = orientation[i]
    sine = (math.sin(math.radians(theta)))
    cosine = (math.cos(math.radians(theta)))

    Q11 = uni_longi_mod / (1 - (poisson_LT * poisson_TL))
    Q22 = trans_mod / (1 - (poisson_LT * poisson_TL))
    Q12 = (poisson_LT * trans_mod) / (1 - (poisson_LT * poisson_TL))
    Q66 = shear_mod

    Q11_bar = Q11 * (cosine ** 4) + Q22 * (sine ** 4) + 2 * (Q12 + (2 * Q66)) * (sine ** 2) * (cosine ** 2)
    Q22_bar = Q11 * (sine ** 4) + Q22 * (cosine ** 4) + 2 * (Q12 + (2 * Q66)) * (sine ** 2) * (cosine ** 2)
    Q12_bar = (Q11 + Q22 - (4 * Q66)) * (sine ** 2) * (cosine ** 2) + Q12 * ((cosine ** 4) + (sine ** 4))
    Q66_bar = (Q11 + Q22 - (2 * Q12) - (2 * Q66)) * (sine ** 2) * (cosine ** 2) + Q66 * ((cosine ** 4) + (sine ** 4))
    Q16_bar = (Q11 + Q12 - (2 * Q66)) * (cosine ** 3) * (sine) - (Q22 - Q12 - (2 * Q66)) * (sine ** 3) * (cosine)
    Q26_bar = (Q11 + Q12 - (2 * Q66)) * (cosine) * (sine ** 3) - (Q22 - Q12 - (2 * Q66)) * (sine) * (cosine ** 3)

    deformation_xy_top = laminate_strains + (Z_elements_from_top[i] * laminate_curvature)
    deformation_xy_bottom = laminate_strains + (Z_elements_from_top[i + 1] * laminate_curvature)
    # print("\n For z =  " + str([i]) + "deformations are: ", deformation)

    # Conversion of deformation matrix from x-y to L_T plane.
    transformation_matrix = np.matrix([
        [(cosine ** 2), (sine ** 2), (sine * cosine)],
        [(sine ** 2), (cosine ** 2), ((-1) * sine * cosine)],
        [((-2) * sine * cosine), (2 * sine * cosine), ((cosine ** 2) - (sine ** 2))]
    ])
    deformation_LT_top = transformation_matrix * deformation_xy_top
    deformation_LT_bottom = transformation_matrix * deformation_xy_bottom
    Q_bar_matrix = np.matrix([
        [Q11_bar, Q12_bar, Q16_bar],
        [Q12_bar, Q22_bar, Q26_bar],
        [Q16_bar, Q26_bar, Q66_bar]])

    stress_values_top_XY_cor = Q_bar_matrix * deformation_xy_top
    stress_values_bottom_XY_cor = Q_bar_matrix * deformation_xy_bottom

    stress_values_top_LT_cor = Q_bar_matrix * deformation_LT_top
    stress_values_bottom_LT_cor = Q_bar_matrix * deformation_LT_bottom
    # print("\n For layer " + str([i+1]) + " stress values on top side are: \n", stress_values_top)
    # print("\n For layer " + str([i+1]) + " stress values on bottom side are: \n", stress_values_bottom)

    location_top = "layer" + str([i + 1]) + " top:  "
    location_bottom = "layer" + str([i + 1]) + " bottom:  "

    list_location.append(location_top)
    list_location.append(location_bottom)

    # Converting stress
    # _values matrix elements into list.

    list_stress_values_X.append(stress_values_top_XY_cor.item(0, 0))
    list_stress_values_X.append(stress_values_bottom_XY_cor.item(0, 0))
    list_stress_values_Y.append(stress_values_top_XY_cor.item(1, 0))
    list_stress_values_Y.append(stress_values_bottom_XY_cor.item(1, 0))
    list_stress_values_XY.append(stress_values_top_XY_cor.item(2, 0))
    list_stress_values_XY.append(stress_values_bottom_XY_cor.item(2, 0))

    list_stress_values_L.append(stress_values_top_LT_cor.item(0, 0))
    list_stress_values_L.append(stress_values_bottom_LT_cor.item(0, 0))
    list_stress_values_T.append(stress_values_top_LT_cor.item(1, 0))
    list_stress_values_T.append(stress_values_bottom_LT_cor.item(1, 0))
    list_stress_values_LT.append(stress_values_top_LT_cor.item(2, 0))
    list_stress_values_LT.append(stress_values_bottom_LT_cor.item(2, 0))

    list_stress_values_top.append(stress_values_top_LT_cor.item(0, 0))
    list_stress_values_top.append(stress_values_top_LT_cor.item(1, 0))
    list_stress_values_top.append(stress_values_top_LT_cor.item(2, 0))

    list_stress_values_bottom.append(stress_values_bottom_LT_cor.item(0, 0))
    list_stress_values_bottom.append(stress_values_bottom_LT_cor.item(1, 0))
    list_stress_values_bottom.append(stress_values_bottom_LT_cor.item(2, 0))

stress_state_in_XY_plane = pd.DataFrame({
    'Location': list_location,
    'Sigma_X': list_stress_values_X,
    'Sigma_Y': list_stress_values_Y,
    'Tau_XY': list_stress_values_XY
})
stress_state_in_XY_plane.set_index('Location', inplace=True)
print("\n----------------------------- Stress State in XY plane ---------------------------------\n")
print(stress_state_in_XY_plane)
print("--------------------------------------------------------------------------------------------")

stress_state_in_LT_plane = pd.DataFrame({
    'Location': list_location,
    'Sigma_L': list_stress_values_L,
    'Sigma_T': list_stress_values_T,
    'Tau_LT': list_stress_values_LT
})
stress_state_in_LT_plane.set_index('Location', inplace=True)
print("\n------------------------------ Stress State in LT plane---------------------------------\n")
print(stress_state_in_LT_plane)
print("------------------------------------------------------------------------------------------")