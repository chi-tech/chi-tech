import subprocess
import os
import sys
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

def format3(number):
    return "{:3d}".format(number)


def format_filename(filename):
    return "{:35s}".format(filename)

def parse_output(out, search_strings_vals_tols):
    global num_failed
    test_passed = True
    for search in search_strings_vals_tols:
        find_str = search[0]
        true_val = search[1]
        tolerance = search[2]

        # start of the string to find (<0 if not found)
        test_str_start = out.find(find_str)
        # end of the string to find
        test_str_end = test_str_start + len(find_str)
        # end of the line at which string was found
        test_str_line_end = out.find("\n", test_str_start)

        test_passed = True
        if test_str_start >= 0:
            # convert value to number
            test_val = float(out[test_str_end:test_str_line_end])
            if not abs(test_val - true_val) < tolerance:
                test_passed = False
        else:
            test_passed = False

    if test_passed:
        print(" - Passed")
    else:
        print(" - FAILED!")
        num_failed += 1
        # print(out)

    return test_passed

def run_test_tacc(file_name, comment, num_procs,
		search_strings_vals_tols):
    test_name = format_filename(file_name) + " " + comment + " " + str(num_procs) + " MPI Processes"
    print("Running Test " + format3(test_number) + " " + test_name, end='', flush=True)
    if print_only: print(""); return
    with open(f"ChiTest/{file_name}.job", 'w') as job_file:
        job_file.write(textwrap.dedent(f"""
            #!/usr/bin/bash
            #
            #SBATCH -J {file_name} # Job name
            #SBATCH -o ChiTest/{file_name}.o # output file
            #SBATCH -e ChiTest/{file_name}.e # error file
            #SBATCH -p skx-normal # Queue (partition) name
            #SBATCH -N {num_procs // 48 + 1} # Total # of nodes
            #SBATCH -n {num_procs} # Total # of mpi tasks
            #SBATCH -t 00:05:00 # Runtime (hh:mm:ss)
            #SBATCH -A Massively-Parallel-R # Allocation name (req'd if you have more than 1)

            export I_MPI_SHM=disable

            ibrun {kpath_to_exe} ChiTest/{file_name}.lua master_export=false
            """
        ).strip())
    os.system(f"sbatch -W ChiTest/{file_name}.job > /dev/null")  # -W means wait for job to finish
    with open(f"ChiTest/{file_name}.o", 'r') as outfile:
        out = outfile.read()

    passed = parse_output(out, search_strings_vals_tols)

    # Cleanup
    if passed:
        os.system(f"rm ChiTest/{file_name}.job ChiTest/{file_name}.o ChiTest/{file_name}.e")

def run_test_local(file_name, comment, num_procs,
		search_strings_vals_tols):
    test_name = format_filename(file_name) + " " + comment + " " + str(num_procs) + " MPI Processes"
    print("Running Test " + format3(test_number) + " " + test_name, end='', flush=True)
    if print_only: print(""); return
    process = subprocess.Popen(["mpiexec", "-np", str(num_procs), kpath_to_exe,
                                "ChiTest/" + file_name + ".lua", "master_export=false"],
                               cwd=kchi_src_pth,
                               stdout=subprocess.PIPE,
                               universal_newlines=True)
    process.wait()
    out, err = process.communicate()

    parse_output(out, search_strings_vals_tols)

def run_test(file_name, comment, num_procs,
		search_strings_vals_tols):
    global test_number
    test_number += 1
    if ((tests_to_run) and (test_number in tests_to_run)) or \
       (not tests_to_run):
        if tacc:
            run_test_tacc(file_name, comment, num_procs,
                search_strings_vals_tols)
        else:
            run_test_local(file_name, comment, num_procs,
                search_strings_vals_tols)

# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Diffusion tests
run_test(
    file_name="Diffusion1D",
    comment="1D Diffusion Test - CFEM",
    num_procs=1,
    search_strings_vals_tols=[["[0]  Max-value=", 2.5, 1.0e-10]])

run_test(
    file_name="Diffusion1D_KBA",
    comment="1D Diffusion Test KBA partitioning - CFEM",
    num_procs=2,
    search_strings_vals_tols=[["[0]  Max-value=", 2.5, 1.0e-10]])

run_test(
    file_name="Diffusion1D_IP",
    comment="1D Diffusion Test - DFEM",
    num_procs=2,
    search_strings_vals_tols=[["[0]  Max-value=", 0.5006523128, 1.0e-4]])

run_test(
    file_name="Diffusion2D_1Poly",
    comment="2D Diffusion Test - CFEM",
    num_procs=1,
    search_strings_vals_tols=[["[0]  Max-value=", 0.29480, 1.0e-4]])

run_test(
    file_name="Diffusion2D_1Poly_IP",
    comment="2D Diffusion Test - DFEM",
    num_procs=4,
    search_strings_vals_tols=[["[0]  Max-value=", 2.5, 1.0e-4]])

run_test(
    file_name="Diffusion2D_2Unstructured",
    comment="2D Diffusion Test Unstr. Mesh - CFEM",
    num_procs=4,
    search_strings_vals_tols=[["[0]  Max-value=", 0.30384, 1.0e-4]])

run_test(
    file_name="Diffusion2D_2Unstructured_IP",
    comment="2D Diffusion Test Unstr. Mesh - DFEM",
    num_procs=4,
    search_strings_vals_tols=[["[0]  Max-value=", 0.29685, 1.0e-4]])

run_test(
    file_name="Diffusion3D_1Poly",
    comment="3D Diffusion Test - CFEM",
    num_procs=1,
    search_strings_vals_tols=[["[0]  Max-value=", 0.29480, 1.0e-4]])

run_test(
    file_name="Diffusion3D_1Poly_IP",
    comment="3D Diffusion Test - DFEM",
    num_procs=4,
    search_strings_vals_tols=[["[0]  Max-value=", 0.29492, 1.0e-4]])

run_test(
    file_name="Diffusion3D_2Ortho",
    comment="3D Diffusion Test Ortho Mesh - CFEM",
    num_procs=1,
    search_strings_vals_tols=[["[0]  Max-value=", 0.29480, 1.0e-4]])

run_test(
    file_name="Diffusion3D_3Unstructured",
    comment="3D Diffusion Test Unstr. Mesh - CFEM",
    num_procs=4,
    search_strings_vals_tols=[["[0]  Max-value=", 0.29499, 1.0e-4]])

run_test(
    file_name="Diffusion3D_3Unstructured_IP",
    comment="3D Diffusion Test Unstr. Mesh - DFEM",
    num_procs=4,
    search_strings_vals_tols=[["[0]  Max-value=", 0.29632, 1.0e-4]])

# run_test(
#     file_name="Diffusion3D_4VTU",
#     comment="3D Diffusion Test VTU Mesh - CFEM",
#     num_procs=4,
#     search_strings_vals_tols=[["[0]  Max-value=", 0.07373, 1.0e-6]])

# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Transport cases
run_test(
    file_name="Transport1D_1",
    comment="1D LinearBSolver Test - PWLD",
    num_procs=3,
    search_strings_vals_tols=[["[0]  Max-value1=", 0.49903, 1.0e-4],
                              ["[0]  Max-value2=", 7.18243e-4, 1.0e-4]])

run_test(
    file_name="Transport2D_1Poly",
    comment="2D LinearBSolver Test - PWLD",
    num_procs=4,
    search_strings_vals_tols=[["[0]  Max-value1=", 0.50758, 1.0e-4],
                              ["[0]  Max-value2=", 2.52527e-04, 1.0e-4]])

run_test(
    file_name="Transport2D_2Unstructured",
    comment="2D LinearBSolver Test Unstructured grid - PWLD",
    num_procs=4,
    search_strings_vals_tols=[["[0]  Max-value1=", 0.51187, 1.0e-4],
                              ["[0]  Max-value2=", 1.42458e-03, 1.0e-4]])

run_test(
    file_name="Transport3D_1a_Extruder",
    comment="3D LinearBSolver Test - PWLD",
    num_procs=4,
    search_strings_vals_tols=[["[0]  Max-value1=", 5.27450e-01, 1.0e-4],
                              ["[0]  Max-value2=", 3.76339e-04, 1.0e-4]])

run_test(
    file_name="Transport3D_1b_Ortho",
    comment="3D LinearBSolver Test - PWLD Reflecting BC",
    num_procs=4,
    search_strings_vals_tols=[["[0]  Max-value1=", 5.28310e-01, 1.0e-4],
                              ["[0]  Max-value2=", 8.04576e-04, 1.0e-4]])

run_test(
    file_name="Transport3D_1Poly_parmetis",
    comment="3D LinearBSolver Test Ortho Grid Parmetis - PWLD",
    num_procs=4,
    search_strings_vals_tols=[["[0]  Max-value1=", 5.27450e-01, 1.0e-4],
                              ["[0]  Max-value2=", 3.76339e-04, 1.0e-4]])

#----------------------------------------------------

run_test(
    file_name="Transport3D_1Poly_qmom_part1",
    comment="3D LinearBSolver Test Source moment writing - PWLD",
    num_procs=4,
    search_strings_vals_tols=[["[0]  Max-value1=", 1.08320e-01, 1.0e-6],
                              ["[0]  Max-value2=", 0.000000000, 1.0e-10]])

run_test(
    file_name="Transport3D_1Poly_qmom_part2",
    comment="3D LinearBSolver Test Source moment reading - PWLD",
    num_procs=4,
    search_strings_vals_tols=[["[0]  Max-value1=", 1.08197e-01, 1.0e-6],
                              ["[0]  Max-value2=", 9.14681e-06, 1.0e-10]])

#----------------------------------------------------

run_test(
    file_name="Transport3D_2Unstructured",
    comment="3D LinearBSolver Test Extruded Unstructured - PWLD",
    num_procs=4,
    search_strings_vals_tols=[["[0]  Max-value1=", 5.41465e-01, 1.0e-4],
                              ["[0]  Max-value2=", 3.78243e-04, 1.0e-4]])

run_test(
    file_name="Transport3D_3BlockPoly_DSA",
    comment="3D LinearBSolver Test WGDSA+TGDSA test - PWLD",
    num_procs=4,
    search_strings_vals_tols=[])

run_test(
    file_name="Transport3D_4Cycles1",
    comment="3D LinearBSolver Test Extruded-Unstructured Mesh - PWLD",
    num_procs=4,
    search_strings_vals_tols=[["[0]  Max-value1=", 5.55349e-01, 1.0e-4],
                              ["[0]  Max-value2=", 3.74343e-04, 1.0e-4]])

run_test(
    file_name="Transport3D_5Cycles2",
    comment="3D LinearBSolver Test STAR-CCM+ mesh - PWLD",
    num_procs=4,
    search_strings_vals_tols=[["[0]  Max-value1=", 6.55396, 1.0e-4],
                              ["[0]  Max-value2=", 1.02943, 1.0e-4]])

run_test(
    file_name="KEigenvalueTransport1D_1G",
    comment="1D KSolver LinearBSolver Test - PWLD",
    num_procs=4,
    search_strings_vals_tols=[["[0]          Final k-eigenvalue    :", 0.99954, 1.0e-5]])

run_test(
    file_name="Transport2DCyl_1Monoenergetic",
    comment="2D LinearBSolver Cylindrical Test mono-energetic - PWLD",
    num_procs=4,
    search_strings_vals_tols=[["[0]  Max-value=", 1.00000, 1.0e-09]])

run_test(
    file_name="Transport2DCyl_2Multigroup",
    comment="2D LinearBSolver Cylindrical Test multi-group - PWLD",
    num_procs=4,
    search_strings_vals_tols=[["[0]  Max-valueG1=", 1.00000, 1.0e-09],
                              ["[0]  Max-valueG2=", 0.25000, 1.0e-09]])

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
    sys.exit(0);
else:
    sys.exit(1); 
