import math

print("\t\t--------------- MECHANICAL PROPERTIES OF UNIDIRECTIONAL COMPOSITES ---------------\n")

print("1. Density\n2. Volume fraction of air\n3. Longitudinal tensile Modulus\n4. Load sharing\n"
      "5. Failure Analysis\n6. Critical volume fraction")
print("7. Transverse Modulus\n8. Transverse Strength\n9. Shear Modulus\n10.Poisson's Ratio\n"
      "11. Longitudinal Compressive strength\n12. Thermal Expansion Coefficient\n13. Moisture Expansion Coefficient")
num = int(input("\nEnter the option:  "))
if (num == 1):
    print("Assumption:\n1. There is no air voids in the composite.\n")
    m1 = float(input("Enter mass of fiber:   "))
    rho1 = float(input("Enter density of fiber:   "))
    m2 = float(input("Enter mass of matrix:   "))
    rho2 = float(input("Enter density of matrix:   "))

    m = m1 + m2
    rho_theo = (rho1*rho2*m) / ((m1*rho2)+(m2*rho1))
    print("\nDensity of Composite: ",rho_theo)

if (num == 2):
    rho_exp = float(input("Enter experimental density of material:   "))
    m1 = float(input("Enter mass of fiber:   "))
    rho1 = float(input("Enter density of fiber:   "))
    m2 = float(input("Enter mass of matrix:   "))
    rho2 = float(input("Enter density of matrix:   "))

    m = m1 + m2
    rho_theo = (rho1 * rho2 * m) / ((m1 * rho2) + (m2 * rho1))
    vol_frac_air = ((rho_theo - rho_exp)/ (rho_theo)) * 100

    print("\nVolume fraction of air: ",vol_frac_air,"%")

if (num == 3):
    print("Assumption:\n1. Fibers are perfectly bonded to the matrix material. ")
    print("2. Fibers leads the same amount of longitudinal strains in composite.")
    print("3. There is no air voids in the composite.\n")
    vol_frac_fiber = float(input("Enter volume fraction of fiber:  "))
    young_mod_fiber = float(input("Enter Young's Modulus (in Mpa) of fiber:  "))
    young_mod_matrix = float(input("Enter Young's Modulus (in Mpa) of Matrix:  "))
    uni_longi_mod = young_mod_matrix * (vol_frac_fiber*((young_mod_fiber/young_mod_matrix)-1)+1)

    print("\nLongitudinal modulus of composite: ",uni_longi_mod,"Mpa")

if (num == 4):
    print("Assumption:\n1. Perfectly bonded composite leads to equal amount of strains on loading.\n ")
    vol_frac_fiber = float(input("Enter volume fraction of fiber:  "))
    young_mod_fiber = float(input("Enter Young's Modulus (in Mpa) of fiber:  "))
    young_mod_matrix = float(input("Enter Young's Modulus (in Mpa) of Matrix:  "))

    ratio = (young_mod_fiber / young_mod_matrix )
    load_share = ratio / ((ratio)+((1/vol_frac_fiber)-1))
    print("\n Fiber to matrix load sharing ratio: ",load_share)

if (num == 5):
    uts_matrix = float(input("Enter ultimate tensile stress of matrix(in Mpa):  "))
    uts_fiber = float(input("Enter ultimate tensile stress of fiber(in Mpa):  "))
    vol_frac_fiber = float(input("Enter volume fraction of fiber:  "))
    stress_matrix_fiber_break = float(input("Enter stress in matrix when fiber breaks(in Mpa)  "))

    V_min = (uts_matrix - stress_matrix_fiber_break)/((uts_fiber+uts_matrix)-stress_matrix_fiber_break)
    print("\nTransition volume fraction of fiber (V_min) : ",V_min)

    if (vol_frac_fiber < V_min):
        uts_composite = (1-vol_frac_fiber) * uts_matrix
        print("The maximum theoretical tensile stress in composite:  ",uts_composite,"Mpa")
    else:
        uts_composite = (uts_fiber * vol_frac_fiber) + (stress_matrix_fiber_break * (1-vol_frac_fiber))
        print("The maximum theoretical tensile stress in composite:  ",uts_composite,"Mpa")

if (num == 6):
    uts_matrix = float(input("Enter ultimate tensile stress of matrix(in Mpa):  "))
    uts_fiber = float(input("Enter ultimate tensile stress of fiber(in Mpa):  "))
    stress_matrix_fiber_break = float(input("Enter stress in matrix when fiber breaks(in Mpa)  "))

    V_crit = (uts_matrix - stress_matrix_fiber_break) / (uts_fiber - stress_matrix_fiber_break)
    print("\n Critical Volume fraction is: ",V_crit)

if (num == 7):
    print("NOTE: HALPIN TSAI method is used here to calculate transverse modulus.")
    vol_frac_fiber = float(input("Enter volume fraction of fiber:  "))
    young_mod_fiber = float(input("Enter Young's Modulus (in Mpa) of fiber:  "))
    young_mod_matrix = float(input("Enter Young's Modulus (in Mpa) of Matrix:  "))
    print("\nChoose the fiber cross section area(CSA):\n1. Circular or squared\n2. Ractangular\n")
    choice = int(input("Enter the CSA:   "))

    if (choice == 1):
        ratio = (young_mod_fiber / young_mod_matrix)
        eta = (ratio - 1) / (ratio + 2)
        trans_mod = young_mod_matrix * (1 + (2*eta*vol_frac_fiber)) / (1 - (eta*vol_frac_fiber))
        print("\n Transverse Modulus value is: ",trans_mod)

    if (choice == 2):
        length = int(input("Enter the CSA side value parallel to the loading axis:  "))
        bredth = int(input("Enter the CSA side value perpendicular to the loading axis:  "))
        ratio = (young_mod_fiber / young_mod_matrix)
        shai = 2 * (length / bredth)
        eta = (ratio - 1) / (ratio + shai)
        trans_mod = young_mod_matrix * (1 + (shai * eta * vol_frac_fiber)) / (1 - (eta * vol_frac_fiber))
        print("\n Transverse Modulus value (in Mpa) is: ", trans_mod)

if (num == 8):
    uts_matrix = float(input("Enter ultimate tensile stress of matrix(in Mpa):  "))
    vol_frac_fiber = float(input("Enter volume fraction of fiber:  "))
    young_mod_fiber = float(input("Enter Young's Modulus (in Mpa) of fiber:  "))
    young_mod_matrix = float(input("Enter Young's Modulus (in Mpa) of Matrix:  "))

    mod_ratio = young_mod_matrix / young_mod_fiber
    var = math.sqrt((4 * vol_frac_fiber) / math.pi)
    SCF = (1 - vol_frac_fiber * (1 - mod_ratio)) / (1 - var * (1 - mod_ratio))
    trans_stress = uts_matrix / SCF
    print("\n Tensile strength in transverse direction(in Mpa) is:  ",trans_stress)


if (num == 9):
    print("NOTE:  HALPIN Tsai Technique is used here to calculate shear modulus.")
    vol_frac_fiber = float(input("Enter volume fraction of fiber:  "))
    shear_mod_fiber = float(input("Enter Shear Modulus (in Mpa) of fiber:  "))
    shear_mod_matrix = float(input("Enter Shear Modulus (in Mpa) of Matrix:  "))
    print("\nChoose the fiber cross section area(CSA):\n1. Circular or squared\n2. Ractangular\n")
    choice = int(input("Enter the CSA:   "))

    if (choice == 1):
        shear_mod_ratio = (shear_mod_fiber / shear_mod_matrix)
        eta = (shear_mod_ratio - 1) / (shear_mod_ratio + 1)
        shear_mod = shear_mod_matrix * (1 + (1 * eta * vol_frac_fiber)) / (1 - (eta * vol_frac_fiber))
        print("\n Shear Modulus value (in Mpa) is: ", shear_mod)

    if (choice == 2):
        length = int(input("Enter the length (bigger side):  "))
        bredth = int(input("Enter the width (smaller side):   "))
        shear_mod_ratio = (shear_mod_fiber / shear_mod_matrix)
        shear_shai = (length / bredth)
        eta = (shear_mod_ratio - 1) / (shear_mod_ratio + shear_shai)
        shear_mod = shear_mod_matrix * (1 + (shear_shai * eta * vol_frac_fiber)) / (1 - (eta * vol_frac_fiber))
        print("\n Shear Modulus value (in Mpa) is: ", shear_mod)

if (num == 10):
    print("Assumption:\n1. strain is considered same in fiber, matrix and composite until failure.\n")
    vol_frac_fiber = float(input("Enter volume fraction of fiber:  "))
    young_mod_fiber = float(input("Enter Young's Modulus (in Mpa) of fiber:  "))
    young_mod_matrix = float(input("Enter Young's Modulus (in Mpa) of Matrix:  "))
    poisson_ratio_fiber = float(input("Enter poisson's ratio of fiber:  "))
    poisson_ratio_matrix = float(input("Enter poisson's ratio of matrix:  "))
    print("\nMajor Poisson's Ratio:  composite loading in 'L' direction, poisson's effect is shown in 'T' direction.")
    print("Minor Poisson's Ratio:  composite loading in 'T' direction, poisson's effect is shown in 'L' direction.")
    print("\n1. Major Poisson's Ratio\n2. Minor Poisson's ratio\n")
    choice = int(input("Choose the option:   "))

    if (choice == 1):
        poisson_LT = (poisson_ratio_fiber * vol_frac_fiber)+(poisson_ratio_matrix*(1-vol_frac_fiber))
        print("\n Mojor Poisson's Ratio is:   ", poisson_LT)

    if (choice == 2):
        uni_longi_mod = young_mod_matrix * (vol_frac_fiber * ((young_mod_fiber / young_mod_matrix) - 1) + 1)
        poisson_LT = (poisson_ratio_fiber * vol_frac_fiber) + (poisson_ratio_matrix * (1 - vol_frac_fiber))
        poisson_TL = (young_mod_matrix / uni_longi_mod) * poisson_LT
        print("\n Minor Poisson's Ratio is:   ",poisson_TL)


if (num == 11):
    vol_frac_fiber = float(input("Enter volume fraction of fiber:  "))
    young_mod_fiber = float(input("Enter Young's Modulus (in Mpa) of fiber:  "))
    young_mod_matrix = float(input("Enter Young's Modulus (in Mpa) of Matrix:  "))
    poisson_ratio_fiber = float(input("Enter poisson's ratio of fiber:  "))
    poisson_ratio_matrix = float(input("Enter poisson's ratio of matrix:  "))
    UTstrain = float(input("Enter ultimate tensile strain of martix:  "))

    uni_longi_mod = young_mod_matrix * (vol_frac_fiber * ((young_mod_fiber / young_mod_matrix) - 1) + 1)
    poisson_LT = (poisson_ratio_fiber * vol_frac_fiber) + (poisson_ratio_matrix * (1 - vol_frac_fiber))
    longi_comp_stress = (UTstrain * uni_longi_mod) / (poisson_LT)
    print("\n Longitudinal Compressive Strength (in Mpa) is:   ",longi_comp_stress)


if (num == 12):
    vol_frac_fiber = float(input("Enter volume fraction of fiber:  "))
    young_mod_fiber = float(input("Enter Young's Modulus (in Mpa) of fiber:  "))
    young_mod_matrix = float(input("Enter Young's Modulus (in Mpa) of Matrix:  "))
    poisson_ratio_fiber = float(input("Enter poisson's ratio of fiber:  "))
    poisson_ratio_matrix = float(input("Enter poisson's ratio of matrix:  "))
    thermal_coff_fiber = float(input("Enter thermal expansion cofficient of fiber:   "))
    thermal_coff_matrix = float(input("Enter thermal expansion cofficient of matrix:   "))
    vol_frac_matrix = float(input("Enter volume fraction of matrix:  "))

    uni_longi_mod = young_mod_matrix * (vol_frac_fiber * ((young_mod_fiber / young_mod_matrix) - 1) + 1)
    poisson_LT = (poisson_ratio_fiber * vol_frac_fiber) + (poisson_ratio_matrix * (1 - vol_frac_fiber))

    print("What do you want to find:\n1. Thermal Expansion Coefficient in Longitudinal Direction\n"
          "2. Thermal Expansion Coefficient in Transverse Direction ")
    choice = int(input("Select the option:   "))

    if (choice == 1):
        var1 = thermal_coff_fiber * young_mod_fiber * vol_frac_fiber
        var2 = thermal_coff_matrix * young_mod_matrix * vol_frac_matrix
        longi_exp_coff = (var1 + var2) / uni_longi_mod
        print("\n Thermal Expansion Coefficient in Longitudinal Direction is:   ", longi_exp_coff)

    if (choice == 2):
        print("Select:\n1. If volume fraction of fiber is less than 0.2\n"
              "2. If volume fraction of fiber is greater than 0.2")
        sub_choice = int(input("\nEnter the option:   "))

        if (sub_choice == 1):
            var1 = thermal_coff_fiber * young_mod_fiber * vol_frac_fiber
            var2 = thermal_coff_matrix * young_mod_matrix * vol_frac_matrix
            longi_exp_coff = (var1 + var2) / uni_longi_mod
            var3 = longi_exp_coff * poisson_LT
            var4 = (1 + poisson_ratio_matrix) * vol_frac_matrix * thermal_coff_matrix
            var5 = (1 + poisson_ratio_fiber) * vol_frac_fiber * thermal_coff_fiber
            trans_exp_coff = var5 + var4 - var3
            print("\n Thermal Expansion Coefficient in Transverse Direction is:   ", trans_exp_coff)

        if (sub_choice == 2):
            var4 = (1 + poisson_ratio_matrix) * vol_frac_matrix * thermal_coff_matrix
            var6 = thermal_coff_fiber * vol_frac_fiber
            trans_exp_coff = var6 + var4
            print("\n Thermal Expansion Coefficient in Transverse Direction is:   ", trans_exp_coff)


if (num == 13):
    rho2 = float(input("Enter density of matrix:   "))
    rho_water = float(input("Enter density of water at same thermodynamic conditions:   "))
    rho_m = float(input("Enter density of matrix after absorbing moisture:   "))
    poisson_ratio_matrix = float(input("Enter poisson's ratio of matrix:  "))
    rho_compo = float(input("Enter total density of composite:   "))

    mois_coff = rho2 / (3 * rho_water)
    print("\nGeneral Matrix Expansion Coefficient is:   ",mois_coff)
    trans_mois = (rho_compo*(1+poisson_ratio_matrix)*mois_coff)/rho_m
    print("Moisture Expansion Coefficient in transverse direction of composite is:   ",trans_mois)
    print("\nNOTE: Longitudinal direction coefficient is considered as zero.\n"
          "Reason:  Fibers are very stiff in that direction,So they don't let matrix to expand.")







