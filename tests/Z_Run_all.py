import subprocess
import os
import sys
import shutil
import textwrap

# To run a range of tests pass a text string via the command line.
# i.e., "tests_to_run=[*range(14,16)]"

# This python script executes the regression test suite.
# In order to add your own test, copy one of the test blocks
# and modify the logic to what you would like tested.

# General guidance: Do not write a regression test that checks
# for number of iterations. Rather make it check for answers
# since something as simple as a preconditioner can change
# iteration count but still produce the same answer

kscript_path = os.path.dirname(os.path.abspath(__file__))
kchi_src_pth = kscript_path + '/../'
kpath_to_exe = kchi_src_pth + '/bin/ChiTech'
tests_to_run = []
print_only   = False

# Get the correct mpiexec executable for the local machine
# This can be overwritten
mpiexec = shutil.which("mpiexec")

print("")
if len(sys.argv) >= 2:
    exec(sys.argv[1])
print("************* ChiTech Regression Test *************")
print("")
test_number = 0
num_failed = 0

# Determine if we are on TACC
# (each test will require a separate job)
hostname = subprocess.check_output(['hostname']).decode('utf-8')
tacc = False
if "tacc.utexas.edu" in hostname:
    tacc = True

print("Using mpiexec at: ", mpiexec)

def format3(number):
    return "{:3d}".format(number)


def format_filename(filename):
    return "{:38s}".format(filename[:38])

# Numerical comparison:
#search[0] = "NumCompare"
#search[1] = string used to identify a line
#search[2] = Which word
#search[3] = numerical format ("float" or "int")
#search[4] = value it should be
#search[5] = tolerance (floats only)

# String comparison:
#search[0] = "StrCompare"
#search[1] = string used to identify a line
#search[2] = Which word
#search[3] = value it should be

def parse_output(out, search_strings_vals_tols):
    global num_failed
    test_passed = True
    for search in search_strings_vals_tols:
        if search[0] == "NumCompare" or search[0] == "StrCompare":
            continue
        find_str = search[0]

        # start of the string to find (<0 if not found)
        test_str_start = out.find(find_str)
        # end of the string to find
        test_str_end = test_str_start + len(find_str)
        # end of the line at which string was found
        test_str_line_end = out.find("\n", test_str_start)

        # test_passed = True
        if test_str_start >= 0:
            if len(search) > 1:
                true_val = search[1]
                tolerance = search[2]
                # convert value to number
                test_val = float(out[test_str_end:test_str_line_end])
                if not abs(test_val - true_val) < tolerance:
                    test_passed = False
                    print("\nTest failed:\nLine:" +
                          out[test_str_start:test_str_line_end] + "\n" +
                          "Test:", search)
                    break
        else:
            test_passed = False
            print("\nTest failed: identifying string not found",
                  search[0])
            break

    for search in search_strings_vals_tols:
        if search[0] != "NumCompare" and search[0] != "StrCompare":
            continue
        identifying_string = search[1]
        # start of the string to find (<0 if not found)
        id_str_start = out.find(identifying_string)
        # end of the line at which string was found
        id_str_line_end = out.find("\n", id_str_start)

        # If we didnt find this string we should quit
        if id_str_start < 0:
            test_passed = False
            print("\nTest failed: identifying string not found\n",
                  identifying_string)
            break

        line = out[id_str_start:id_str_line_end]
        words = line.split()

        if len(search) == 2:
            continue

        word_number = search[2]
        if word_number >= len(words):
            test_passed = False
            print("\nTest failed Comparison work not found:\nLine:" +
                  line + "\n" +
                  "Test:", search)
            break

        if search[0] == "NumCompare":
            numerical_format   = search[3]
            value_it_should_be = search[4]

            trial_value_str = words[word_number]

            if numerical_format == "float":
                trial_value = float(trial_value_str)
                tolerance = search[5]

                if abs(trial_value-value_it_should_be) > tolerance:
                    test_passed = False
                    print("\nTest failed:\nLine:" +
                          line + "\n" +
                          "Test:", search)
                    break

            if numerical_format == "int":
                trial_value = int(trial_value_str)

                if trial_value != value_it_should_be:
                    test_passed = False
                    print("\nTest failed:\nLine:" +
                          line + "\n" +
                          "Test:", search)
                    break

        if search[0] == "StrCompare":
            value_it_should_be = search[3]

            trial_word = words[word_number]

            if trial_word != value_it_should_be:
                test_passed = False
                print("\nTest failed:\nLine:" +
                      line + "\n" +
                      "Test:", search)
                break


    if test_passed:
        print(" - Passed")
    else:
        print(" - FAILED!")
        num_failed += 1
        # print(out)

    return test_passed

def run_test_tacc(file_name, comment, num_procs, search_strings_vals_tols):
    test_name = format_filename(file_name) + " " + comment + " " \
                + str(num_procs) + " MPI Processes"
    print("Running Test " + format3(test_number) + " " + test_name, end='', flush=True)
    if print_only: print(""); return
    with open(f"tests/{file_name}.job", 'w') as job_file:
        job_file.write(textwrap.dedent(f"""
            #!/usr/bin/bash
            #
            #SBATCH -J {file_name} # Job name
            #SBATCH -o tests/{file_name}.o # output file
            #SBATCH -e tests/{file_name}.e # error file
            #SBATCH -p skx-normal # Queue (partition) name
            #SBATCH -N {num_procs // 48 + 1} # Total # of nodes
            #SBATCH -n {num_procs} # Total # of mpi tasks
            #SBATCH -t 00:05:00 # Runtime (hh:mm:ss)
            #SBATCH -A Massively-Parallel-R # Allocation name (req'd if you have more than 1)

            export I_MPI_SHM=disable

            ibrun {kpath_to_exe} tests/{file_name}.lua master_export=false
            """
        ).strip())
    os.system(f"sbatch -W tests/{file_name}.job > /dev/null")  # -W means wait for job to finish
    with open(f"tests/{file_name}.o", 'r') as outfile:
        out = outfile.read()

    passed = parse_output(out, search_strings_vals_tols)

    # Cleanup
    if passed:
        os.system(f"rm tests/{file_name}.job tests/{file_name}.o tests/{file_name}.e")


def run_test_local(file_name, comment, num_procs, search_strings_vals_tols):
    test_name = format_filename(file_name) \
                + " - " + format_filename(comment) + " - " \
                + str(num_procs) + " MPI Processes"
    print("Running Test " + format3(test_number) + " " + test_name, end='', flush=True)
    if print_only: print(""); return

    process = subprocess.Popen([mpiexec, "-np", str(num_procs), kpath_to_exe,
                                "tests/" + file_name + ".lua", "master_export=false"],
                               cwd=kchi_src_pth,
                               stdout=subprocess.PIPE,
                               universal_newlines=True)
    process.wait()
    out, err = process.communicate()

    parse_output(out, search_strings_vals_tols)


def run_test(file_name, comment, num_procs, search_strings_vals_tols):
    global test_number
    test_number += 1
    if ((tests_to_run) and (test_number in tests_to_run)) or (not tests_to_run):
        if tacc:
            run_test_tacc(file_name, comment, num_procs, search_strings_vals_tols)
        else:
            run_test_local(file_name, comment, num_procs, search_strings_vals_tols)


# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Diffusion tests
#
### CFEM diffusion tests
run_test(
    file_name="CFEM_Diffusion/Diffusion_2D_1a_linear",
    comment="2D Diffusion with linear solution",
    num_procs=1,
    search_strings_vals_tols=[["[0]  Max-value=", 2.666667, 1.0e-10]])

run_test(
    file_name="CFEM_Diffusion/Diffusion_2D_2a_DirBCs",
    comment="2D Diffusion with Dirichlet BC",
    num_procs=1,
    search_strings_vals_tols=[["[0]  Avg-value=", 0.295902, 1.0e-10]])

run_test(
    file_name="CFEM_Diffusion/Diffusion_2D_2b_RobinBCs",
    comment="2D Diffusion with Robin BC",
    num_procs=1,
    search_strings_vals_tols=[["[0]  Avg-value=", 0.241751, 1.0e-10]])

run_test(
    file_name="CFEM_Diffusion/Diffusion_2D_3a_analytical_coef",
    comment="2D Diffusion with Analytical Coefficients",
    num_procs=1,
    search_strings_vals_tols=[["[0]  Max-value=", 0.021921, 1.0e-10]])

run_test(
    file_name="CFEM_Diffusion/Diffusion_2D_3b_analytical_coef2",
    comment="2D Diffusion with Manufactured Solution",
    num_procs=1,
    search_strings_vals_tols=[["[0]  Max-value=", 1.000244, 1.0e-10]])

### DFEM diffusion tests
run_test(
    file_name="DFEM_Diffusion/Diffusion_2D_1a_linear",
    comment="2D Diffusion with linear solution",
    num_procs=1,
    search_strings_vals_tols=[["[0]  Max-value=", 2.666667, 1.0e-10]])

run_test(
    file_name="DFEM_Diffusion/Diffusion_2D_2a_DirBCs",
    comment="2D Diffusion with Dirichlet BC",
    num_procs=1,
    search_strings_vals_tols=[["[0]  Avg-value=", 0.295892, 1.0e-10]])

run_test(
    file_name="DFEM_Diffusion/Diffusion_2D_2b_RobinBCs",
    comment="2D Diffusion with Robin BC",
    num_procs=1,
    search_strings_vals_tols=[["[0]  Avg-value=", 0.241757, 1.0e-10]])

run_test(
    file_name="DFEM_Diffusion/Diffusion_2D_3a_analytical_coef",
    comment="2D Diffusion with Analytical Coefficients",
    num_procs=1,
    search_strings_vals_tols=[["[0]  Max-value=", 0.021923, 1.0e-10]])

run_test(
    file_name="DFEM_Diffusion/Diffusion_2D_3b_analytical_coef2",
    comment="2D Diffusion with Manufactured Solution",
    num_procs=1,
    search_strings_vals_tols=[["[0]  Max-value=", 1.000586, 1.0e-10]])

#1
run_test(
    file_name="Diffusion1D",
    comment="1D Diffusion Test - CFEM",
    num_procs=1,
    search_strings_vals_tols=[["[0]  Max-value=", 2.5, 1.0e-10]])

#2
run_test(
    file_name="Diffusion1D_KBA",
    comment="1D Diffusion Test KBA partitioning - CFEM",
    num_procs=2,
    search_strings_vals_tols=[["[0]  Max-value=", 2.5, 1.0e-10]])

#3
run_test(
    file_name="Diffusion1D_IP",
    comment="1D Diffusion Test - DFEM",
    num_procs=2,
    search_strings_vals_tols=[["[0]  Max-value=", 0.5006523128, 1.0e-4]])

#4
run_test(
    file_name="Diffusion2D_1Poly",
    comment="2D Diffusion Test - CFEM",
    num_procs=1,
    search_strings_vals_tols=[["[0]  Max-value=", 0.29480, 1.0e-4]])

#5
run_test(
    file_name="Diffusion2D_1Poly_IP",
    comment="2D Diffusion Test - DFEM",
    num_procs=4,
    search_strings_vals_tols=[["[0]  Max-value=", 2.5, 1.0e-4]])

#6
run_test(
    file_name="Diffusion2D_2Unstructured",
    comment="2D Diffusion Test Unstr. Mesh - CFEM",
    num_procs=4,
    search_strings_vals_tols=[["[0]  Max-value=", 0.30384, 1.0e-4]])

#7
run_test(
    file_name="Diffusion2D_2Unstructured_IP",
    comment="2D Diffusion Test Unstr. Mesh - DFEM",
    num_procs=4,
    search_strings_vals_tols=[["[0]  Max-value=", 0.29685, 1.0e-4]])

#8
run_test(
    file_name="Diffusion3D_1Poly",
    comment="3D Diffusion Test - CFEM",
    num_procs=1,
    search_strings_vals_tols=[["[0]  Max-value=", 0.29480, 1.0e-4]])

#9
run_test(
    file_name="Diffusion3D_1Poly_IP",
    comment="3D Diffusion Test - DFEM",
    num_procs=4,
    search_strings_vals_tols=[["[0]  Max-value=", 0.29492, 1.0e-4]])

#10
run_test(
    file_name="Diffusion3D_2Ortho",
    comment="3D Diffusion Test Ortho Mesh - CFEM",
    num_procs=1,
    search_strings_vals_tols=[["[0]  Max-value=", 0.29480, 1.0e-4]])

#11
run_test(
    file_name="Diffusion3D_3Unstructured",
    comment="3D Diffusion Test Unstr. Mesh - CFEM",
    num_procs=4,
    search_strings_vals_tols=[["[0]  Max-value=", 0.29499, 1.0e-4]])

#12
run_test(
    file_name="Diffusion3D_3Unstructured_IP",
    comment="3D Diffusion Test Unstr. Mesh - DFEM",
    num_procs=4,
    search_strings_vals_tols=[["[0]  Max-value=", 0.29632, 1.0e-4]])

# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Transport cases
#13
run_test(
    file_name="Transport1D_1",
    comment="1D LinearBSolver Test - PWLD",
    num_procs=3,
    search_strings_vals_tols=[["[0]  Max-value1=", 0.49903, 1.0e-4],
                              ["[0]  Max-value2=", 7.18243e-4, 1.0e-4]])

#14
run_test(
    file_name="Transport1D_3a_DSA_ortho",
    comment="1D LinearBSolver test of a block of graphite with an air cavity. DSA and TG",
    num_procs=4,
    search_strings_vals_tols=[["StrCompare", "WGS groups [0-62] Iteration    28", 7, "CONVERGED"],
                              ["StrCompare", "WGS groups [63-167] Iteration    55", 7, "CONVERGED"],
                              ["NumCompare", "WGS groups [0-62] Iteration    28", 6, "float", 6.7433e-07, 1.0e-9],
                              ["NumCompare", "WGS groups [63-167] Iteration    55", 6, "float", 5.67431e-07, 1.0e-9]])

#15
run_test(
    file_name="Transport2D_1Poly",
    comment="2D LinearBSolver Test - PWLD",
    num_procs=4,
    search_strings_vals_tols=[["[0]  Max-value1=", 0.50758, 1.0e-4],
                              ["[0]  Max-value2=", 2.52527e-04, 1.0e-4]])

#16
run_test(
    file_name="Transport2D_2Unstructured",
    comment="2D LinearBSolver Test Unstructured grid - PWLD",
    num_procs=4,
    search_strings_vals_tols=[["[0]  Max-value1=", 0.51187, 1.0e-4],
                              ["[0]  Max-value2=", 1.42458e-03, 1.0e-4]])

#17
run_test(
    file_name="Transport2D_3Poly_quad_mod",
    comment="2D LinearBSolver Test Polar-Optimized quadrature - PWLD",
    num_procs=4,
    search_strings_vals_tols=[["[0]  Max-value1=", 0.50758, 1.0e-4],
                              ["[0]  Max-value2=", 2.52527e-04, 1.0e-4]])

#18
run_test(
    file_name="Transport2D_4a_DSA_ortho",
    comment="2D LinearBSolver test of a block of graphite with an air cavity. DSA and TG",
    num_procs=4,
    search_strings_vals_tols=[["StrCompare", "WGS groups [0-62] Iteration    53", 7, "CONVERGED"],
                              ["StrCompare", "WGS groups [63-167] Iteration    59", 7, "CONVERGED"],
                              ["NumCompare", "WGS groups [0-62] Iteration    53", 6, "float", 6.01304e-07, 1.0e-9],
                              ["NumCompare", "WGS groups [63-167] Iteration    59", 6, "float", 6.21411e-07, 1.0e-9]])

#19
run_test(
    file_name="Transport2D_4b_DSA_ortho",
    comment="2D LinearBSolver test of a block of graphite with an air cavity. DSA and TG",
    num_procs=4,
    search_strings_vals_tols=[["StrCompare", "WGS groups [0-62] Iteration    54", 7, "CONVERGED"],
                              ["StrCompare", "WGS groups [63-167] Iteration    57", 7, "CONVERGED"],
                              ["NumCompare", "WGS groups [0-62] Iteration    54", 6, "float", 4.97136e-07, 1.0e-9],
                              ["NumCompare", "WGS groups [63-167] Iteration    57", 6, "float", 6.88134e-07, 1.0e-9]])

#20
run_test(
    file_name="Transport3D_1a_Extruder",
    comment="3D LinearBSolver Test - PWLD",
    num_procs=4,
    search_strings_vals_tols=[["[0]  Max-value1=", 5.27450e-01, 1.0e-4],
                              ["[0]  Max-value2=", 3.76339e-04, 1.0e-4]])

#21
run_test(
    file_name="Transport3D_1b_Ortho",
    comment="3D LinearBSolver Test - PWLD Reflecting BC",
    num_procs=4,
    search_strings_vals_tols=[["[0]  Max-value1=", 5.28310e-01, 1.0e-4],
                              ["[0]  Max-value2=", 8.04576e-04, 1.0e-4]])

#22
run_test(
    file_name="Transport3D_1Poly_parmetis",
    comment="3D LinearBSolver Test Ortho Grid Parmetis - PWLD",
    num_procs=4,
    search_strings_vals_tols=[["[0]  Max-value1=", 5.27450e-01, 1.0e-4],
                              ["[0]  Max-value2=", 3.76339e-04, 1.0e-4]])

#----------------------------------------------------
#23
run_test(
    file_name="Transport3D_1Poly_qmom_part1",
    comment="3D LinearBSolver Test Source moment writing - PWLD",
    num_procs=4,
    search_strings_vals_tols=[["[0]  Max-value1=", 1.08320e-01, 1.0e-6],
                              ["[0]  Max-value2=", 0.000000000, 1.0e-10]])

#24
run_test(
    file_name="Transport3D_1Poly_qmom_part2",
    comment="3D LinearBSolver Test Source moment reading - PWLD",
    num_procs=4,
    search_strings_vals_tols=[["[0]  Max-value1=", 1.01701e-04, 1.0e-6],
                              ["[0]  Max-value2=", 9.14681e-06, 1.0e-10]])

#----------------------------------------------------
#25
run_test(
    file_name="Transport3D_2Unstructured",
    comment="3D LinearBSolver Test Extruded Unstructured - PWLD",
    num_procs=4,
    search_strings_vals_tols=[["[0]  Max-value1=", 5.41465e-01, 1.0e-4],
                              ["[0]  Max-value2=", 3.78243e-04, 1.0e-4]])

#26
run_test(
    file_name="Transport3D_3a_DSA_ortho",
    comment="3D LinearBSolver test of a block of graphite with an air cavity. DSA and TG",
    num_procs=4,
    search_strings_vals_tols=[["StrCompare", "WGS groups [0-62] Iteration    54", 7, "CONVERGED"],
                              ["StrCompare", "WGS groups [63-167] Iteration    69", 7, "CONVERGED"],
                              ["NumCompare", "WGS groups [0-62] Iteration    54", 6, "float", 7.88852e-07, 1.0e-9],
                              ["NumCompare", "WGS groups [63-167] Iteration    69", 6, "float", 9.78723e-07, 1.0e-9]])

#27
run_test(
    file_name="Transport3D_4Cycles1",
    comment="3D LinearBSolver Test Extruded-Unstructured Mesh - PWLD",
    num_procs=4,
    search_strings_vals_tols=[["[0]  Max-value1=", 5.55349e-01, 1.0e-4],
                              ["[0]  Max-value2=", 3.74343e-04, 1.0e-4]])

#28
run_test(
    file_name="Transport3D_5Cycles2",
    comment="3D LinearBSolver Test STAR-CCM+ mesh - PWLD",
    num_procs=4,
    search_strings_vals_tols=[["[0]  Max-value1=", 6.55396, 1.0e-4],
                              ["[0]  Max-value2=", 1.02943, 1.0e-4]])

#29
run_test(
    file_name="KEigenvalueTransport1D_1G",
    comment="1D KSolver LinearBSolver Test - PWLD",
    num_procs=4,
    search_strings_vals_tols=[["[0]          Final k-eigenvalue    :", 0.99954, 1.0e-5]])

#30
run_test(
    file_name="Transport2DCyl_1Monoenergetic",
    comment="2D LinearBSolver Cylindrical Test mono-energetic - PWLD",
    num_procs=4,
    search_strings_vals_tols=[["[0]  Max-value=", 1.00000, 1.0e-09]])

#31
run_test(
    file_name="Transport2DCyl_2Multigroup",
    comment="2D LinearBSolver Cylindrical Test multi-group - PWLD",
    num_procs=4,
    search_strings_vals_tols=[["[0]  Max-valueG1=", 1.00000, 1.0e-09],
                              ["[0]  Max-valueG2=", 0.25000, 1.0e-09]])

#------------------------------------------------ Adjoints
#32
run_test(
    file_name="Adjoint2D_1a_forward",
    comment="2D Transport test with localized material source FWD",
    num_procs=4,
    search_strings_vals_tols=[["QOI-value=", 1.38397e-05, 1.0e-08]])

#33
run_test(
    file_name="Adjoint2D_1b_adjoint",
    comment="2D Transport test with localized material source Adjoint generation",
    num_procs=4,
    search_strings_vals_tols=[])

#34
run_test(
    file_name="Adjoint2D_1c_response",
    comment="2D Transport test with localized material source Adjoint inner product",
    num_procs=4,
    search_strings_vals_tols=[["Inner-product=", 1.38405e-05, 1.0e-08]])


#35
run_test(
    file_name="Adjoint2D_2a_forward",
    comment="2D Transport test with point source FWD",
    num_procs=4,
    search_strings_vals_tols=[["QOI-value=", 2.90386e-05 , 1.0e-08]])

#36
run_test(
    file_name="Adjoint2D_2b_adjoint",
    comment="2D Transport test with point source Adjoint generation",
    num_procs=4,
    search_strings_vals_tols=[])

#37
run_test(
    file_name="Adjoint2D_2c_response",
    comment="2D Transport test with point source Adjoint response",
    num_procs=4,
    search_strings_vals_tols=[["Inner-product=", 2.90543e-05, 1.0e-08]])


#38
run_test(
    file_name="Adjoint2D_3a_forward",
    comment="2D Transport test with point source Multigroup FWD",
    num_procs=4,
    search_strings_vals_tols=[["QOI-value[0]=", 1.12687e-06, 1.0e-09],
                              ["QOI-value[1]=", 2.95934e-06, 1.0e-09],
                              ["QOI-value[2]=", 3.92975e-06, 1.0e-09],
                              ["QOI-value[3]=", 4.18474e-06, 1.0e-09],
                              ["QOI-value[4]=", 3.89649e-06, 1.0e-09],
                              ["QOI-value[5]=", 3.30482e-06, 1.0e-09],
                              ["QOI-value[6]=", 1.54506e-06, 1.0e-09],
                              ["QOI-value[7]=", 6.74868e-07, 1.0e-09],
                              ["QOI-value[8]=", 3.06178e-07, 1.0e-09],
                              ["QOI-value[9]=", 2.07284e-07, 1.0e-09],
                              ["QOI-value[sum]=", 2.21354e-05, 1.0e-08]])

#39
run_test(
    file_name="Adjoint2D_3b_adjoint",
    comment="2D Transport test with point source Multigroup Adjoint generation",
    num_procs=4,
    search_strings_vals_tols=[])

#40
run_test(
    file_name="Adjoint2D_3c_response",
    comment="2D Transport test with point source Multigroup Adjoint Response",
    num_procs=4,
    search_strings_vals_tols=[["Inner-product=", 3.30607e-06, 1.0e-09]])

# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ END OF TESTS
print("")
if num_failed == 0:
    print("All regression tests passed!")
else:
    print("ERROR: Not all regression tests passed!")
    print("Number of Tests failed = " + str(num_failed))
print("")
print("************* End of Regression Test *************")
print("")
if num_failed == 0:
    sys.exit(0)
else:
    sys.exit(1)
