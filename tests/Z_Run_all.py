import subprocess
import os
import sys
import shutil
import textwrap
import argparse


class CustomFormatter(argparse.RawTextHelpFormatter,
                      argparse.MetavarTypeHelpFormatter,
                      argparse.ArgumentDefaultsHelpFormatter):
    pass

arguments_help = textwrap.dedent('''\
Run the regression suite, or optionally, a series of individual
tests or a range of tests. 

To run all tests, use
python test/Z_Run_all.py

To run a list of individual tests, use 
python test/Z_Run_all.py --test-list test1 test2 ... testN

To run a range of tests, use
python test/Z_Run_all.py --test-range test_start test_end

For verbose fail messages, use
python test/Z_Run_all.py --verbose_fail

In each of the above examples, all arguments are separated by spaces 
and each of the arguments are integers.
''')

print(arguments_help)

parser = argparse.ArgumentParser(
    description="A script to run the regression test suite.",
    formatter_class=CustomFormatter,
    epilog=arguments_help
)

parser.add_argument(
    '--test-list',
    nargs='*',
    type=int,
    required=False,
    help="A list of test IDs to run."
)

parser.add_argument(
    '--test-range',
    nargs=2,
    default=None,
    type=int,
    required=False,
    help="The first and last test ID of a range of tests to run."
)

parser.add_argument(
    '--verbose_fail',
    action="store_true",
    help="If set, prints the error stream if tests failed."
)

argv = parser.parse_args()

# To run a particular set of tests, pass "--test-list" followed by space
# separated integers corresponding to the test numbers to run.

# To run a particular range of tests, pass "--test-range" followed by
# the first test desired and the last test desired.

# This python script executes the regression test suite.
# In order to add your own test, copy one of the test blocks
# and modify the logic to what you would like tested.

# General guidance: Do not write a regression test that checks
# for number of iterations. Rather make it check for answers
# since something as simple as a preconditioner can change
# iteration count but still produce the same answer

kscript_path = os.path.dirname(os.path.abspath(__file__))
kchi_src_pth = os.path.join(kscript_path, "..")
kpath_to_exe = os.path.join(kchi_src_pth, "bin/ChiTech")
print_only = False

# Get the correct mpiexec executable for the local machine
# This can be overwritten
mpiexec = shutil.which("mpiexec")

# Create the list of tests to run
if argv.test_list and argv.test_range:
    raise ValueError(
        "Either a list of individual tests or the bounds "
        "of the range of tests can be provided, not both."
    )

if argv.test_list:
    tests_to_run = argv.test_list
elif argv.test_range:
    first, last = argv.test_range
    tests_to_run = list(range(first, last + 1))
else:
    tests_to_run = list()

print()
print("************* ChiTech Regression Test *************")
print()

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
    return f"{number:3d}"


def format_filename(filename):
    return f"{filename[:50]:50s}"


# Numerical comparison:
# search[0] = "NumCompare"
# search[1] = string used to identify a line
# search[2] = Which word
# search[3] = numerical format ("float" or "int")
# search[4] = value it should be
# search[5] = tolerance (floats only)

# String comparison:
# search[0] = "StrCompare"
# search[1] = string used to identify a line
# search[2] = Which word
# search[3] = value it should be

def parse_output(out, err, search_strings_vals_tols):
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
                    if argv.verbose_fail:
                        print("\nTest failed:\nLine:" +
                              out[test_str_start:test_str_line_end] + "\n" +
                              "Test:", search)
                    break
        else:
            test_passed = False
            if argv.verbose_fail:
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
            if argv.verbose_fail:
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
            if argv.verbose_fail:
                print("\nTest failed Comparison word not found:\nLine:" +
                      line + "\n" +
                      "Test:", search)
            break

        if search[0] == "NumCompare":
            numerical_format = search[3]
            value_it_should_be = search[4]

            trial_value_str = words[word_number]

            if numerical_format == "float":
                trial_value = float(trial_value_str)
                tolerance = search[5]

                if abs(trial_value - value_it_should_be) > tolerance:
                    test_passed = False
                    if argv.verbose_fail:
                        print("\nTest failed:\nLine:" +
                              line + "\n" +
                              "Test:", search)
                    break

            if numerical_format == "int":
                trial_value = int(trial_value_str)

                if trial_value != value_it_should_be:
                    test_passed = False
                    if argv.verbose_fail:
                        print("\nTest failed:\nLine:" +
                              line + "\n" +
                              "Test:", search)
                    break

        if search[0] == "StrCompare":
            value_it_should_be = search[3]

            trial_word = words[word_number]

            if trial_word != value_it_should_be:
                test_passed = False
                if argv.verbose_fail:
                    print("\nTest failed:\nLine:" +
                          line + "\n" +
                          "Test:", search)
                break

    if test_passed:
        print(" - \033[32mPassed\033[39m")
    else:
        print(" - \033[31mFAILED!\033[39m")
        num_failed += 1
        if argv.verbose_fail:
            print(err)

    return test_passed


def run_test_tacc(file_name, comment, num_procs, search_strings_vals_tols):
    test_name = f"{format_filename(file_name)} {comment} " \
                f"{num_procs} MPI Processes"

    msg = f"Running Test {format3(test_number)} {test_name}"
    print(msg, end='', flush=True)

    if print_only:
        print()
        return

    with open(f"tests/{file_name}.job", 'w') as job_file:
        job_file.write(textwrap.dedent(f"""
            #!/usr/bin/bash
            #
            #SBATCH -J {file_name} # Job name
            #SBATCH -o tests/{file_name}.o # output file
            #SBATCH -e tests/{file_name}.e # error file
            #SBATCH -p skx-normal_ # Queue (partition) name
            #SBATCH -N {num_procs // 48 + 1} # Total # of nodes
            #SBATCH -n {num_procs} # Total # of mpi tasks
            #SBATCH -t 00:05:00 # Runtime (hh:mm:ss)
            #SBATCH -A Massively-Parallel-R # Allocation name (req'd if you have more than 1)

            export I_MPI_SHM=disable

            ibrun {kpath_to_exe} tests/{file_name}.lua master_export=false
            """
                                       ).strip())

    # -W means wait for job to finish
    os.system(f"sbatch -W tests/{file_name}.job > /dev/null")
    with open(f"tests/{file_name}.o", 'r') as outfile:
        out = outfile.read()

    passed = parse_output(out, search_strings_vals_tols)

    # Cleanup
    if passed:
        os.system(f"rm tests/{file_name}.job "
                  f"tests/{file_name}.o tests/{file_name}.e")


def run_test_local(file_name, comment, num_procs, search_strings_vals_tols):
    test_name = f"{format_filename(file_name)} - {comment} - " \
                f"{num_procs} MPI Processes"

    msg = f"Running Test {format3(test_number)} {test_name}"
    print(msg, end='', flush=True)
    if print_only:
        print()
        return

    cmd = f"{mpiexec} -np {num_procs} {kpath_to_exe} " \
          f"tests/{file_name}.lua master_export=false"
    process = subprocess.Popen(cmd.split(),
                               cwd=kchi_src_pth,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               universal_newlines=True)
    process.wait()
    out, err = process.communicate()

    parse_output(out, err, search_strings_vals_tols)



def run_test(file_name, comment, num_procs, search_strings_vals_tols):
    global test_number

    test_number += 1
    if (tests_to_run and test_number in tests_to_run) or (not tests_to_run):
        if tacc:
            run_test_tacc(file_name, comment, num_procs,
                          search_strings_vals_tols)
        else:
            run_test_local(file_name, comment, num_procs,
                           search_strings_vals_tols)


def run_unit_test(file_name, comment, num_procs):
    global test_number
    test_number += 1

    global num_failed
    test_name = f"{format_filename(file_name)} - {comment} - " \
                f"{num_procs} MPI Processes"

    msg = f"Running Test {format3(test_number)} {test_name}"
    print(msg, end='', flush=True)
    if print_only:
        print()
        return

    cmd = f"./tests/unit_test_inputs/ZUnitTestRunner.sh {mpiexec} " \
          f"{num_procs} {kpath_to_exe} tests/{file_name}"
    process = subprocess.Popen(cmd.split(),
                               cwd=kchi_src_pth,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               universal_newlines=True)
    process.wait()
    out, err = process.communicate()

    if process.returncode == 0:
        print(" - \033[32mPassed\033[39m")
    else:
        print(" - \033[31mFAILED!\033[39m")
        num_failed += 1
        if argv.verbose_fail:
            print(err)

# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#                             Unit Tests
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
run_unit_test(file_name="unit_test_inputs/test_00_chi_math.lua",
              comment="Testing objects in chi_math",
              num_procs=1)
run_unit_test(file_name="unit_test_inputs/test_01_chi_misc_utils.lua",
              comment="Testing objects in chi_misc_utils",
              num_procs=1)
run_unit_test(file_name="unit_test_inputs/test_02_chi_data_types.lua",
              comment="Testing chi_data_types",
              num_procs=2)
run_unit_test(file_name="unit_test_inputs/test_02a_wdd_ijk_sweep.lua",
              comment="Testing utilities for orthogonal meshes",
              num_procs=1)
run_unit_test(file_name="unit_test_inputs/test_02b_paramblock.lua",
              comment="Testing parameter blocks",
              num_procs=1)

# SimTests
run_test(
    file_name="unit_test_inputs/SimTest01_FV",
    comment="Finite Volume developer's example",
    num_procs=1,
    search_strings_vals_tols=[
        ["NumCompare", "iteration   10", 4, "float", 2.2349712e-07, 1.0e-10]])

run_test(
    file_name="unit_test_inputs/SimTest02_FV",
    comment="Finite Volume gradient developer's example",
    num_procs=1,
    search_strings_vals_tols=[
        ["NumCompare", "iteration   10", 4, "float", 2.2349712e-07, 1.0e-10]])

run_test(
    file_name="unit_test_inputs/SimTest03_PWLC",
    comment="CFEM PWLC developer's example",
    num_procs=4,
    search_strings_vals_tols=[
        ["NumCompare", "iteration   10", 4, "float", 1.3725535e-07, 1.0e-11]])

run_test(
    file_name="unit_test_inputs/SimTest04_PWLC",
    comment="CFEM PWLC with MMS developer's example",
    num_procs=4,
    search_strings_vals_tols=[
        ["NumCompare", "iteration    6", 4, "float", 1.1613008e-07, 1.0e-10],
        ["NumCompare", "Error:"        , 1, "float", 7.032369e-02, 1.0e-7]])

run_test(
    file_name="unit_test_inputs/SimTest06_WDD",
    comment="FV Transport WDD developer's example",
    num_procs=1,
    search_strings_vals_tols=[
        ["NumCompare", "Iteration     8", 2, "float", 5.483e-07, 1.0e-10]])

run_test(
    file_name="unit_test_inputs/SimTest91_PWLD",
    comment="DFEM Transport PWLD developer's example",
    num_procs=1,
    search_strings_vals_tols=[
        ["NumCompare", "Iteration     8", 2, "float", 5.872e-07, 1.0e-10]])

run_test(
    file_name="unit_test_inputs/SimTest92_DSA",
    comment="MIP Diffusion solver tests",
    num_procs=1,
    search_strings_vals_tols=[
        ["NumCompare", "SimTest92_DSA iteration    8", 5, "float", 1.9731053e-11, 1.0e-14]])

run_test(
    file_name="unit_test_inputs/SimTest93_RayTracing",
    comment="Raytracing Transport developer's example",
    num_procs=1,
    search_strings_vals_tols=[])

# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#                             Mesh IO tests
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Tests 14 - 14
run_test(
    file_name="MeshIO/ReadWavefrontObj1",
    comment="Mesh reading 2D Wavefront.obj",
    num_procs=4,
    search_strings_vals_tols=[
        ["StrCompare", "VolumeMesherPredefinedUnpartitioned: Cells created = 3242"]])

# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#                             CFEM Diffusion Tests
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Tests 15 - 19

run_test(
    file_name="CFEM_Diffusion/cDiffusion_2D_1a_linear",
    comment="2D Diffusion with linear solution",
    num_procs=1,
    search_strings_vals_tols=[["[0]  Max-value=", 2.666667, 1.0e-10]])

run_test(
    file_name="CFEM_Diffusion/cDiffusion_2D_2a_DirBCs",
    comment="2D Diffusion with Dirichlet BC",
    num_procs=1,
    search_strings_vals_tols=[["[0]  Avg-value=", 0.295902, 1.0e-10]])

run_test(
    file_name="CFEM_Diffusion/cDiffusion_2D_2b_RobinBCs",
    comment="2D Diffusion with Robin BC",
    num_procs=1,
    search_strings_vals_tols=[["[0]  Avg-value=", 0.241751, 1.0e-10]])

run_test(
    file_name="CFEM_Diffusion/cDiffusion_2D_3a_analytical_coef",
    comment="2D Diffusion with Analytical Coefficients",
    num_procs=1,
    search_strings_vals_tols=[["[0]  Max-value=", 0.021921, 1.0e-10]])

run_test(
    file_name="CFEM_Diffusion/cDiffusion_2D_3b_analytical_coef2",
    comment="2D Diffusion with Manufactured Solution",
    num_procs=1,
    search_strings_vals_tols=[["[0]  Max-value=", 1.000244, 1.0e-10]])

# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#                             CFEM Diffusion Tests
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Tests 20 - 24

run_test(
    file_name="DFEM_Diffusion/dDiffusion_2D_1a_linear",
    comment="2D Diffusion with linear solution",
    num_procs=1,
    search_strings_vals_tols=[["[0]  Max-value=", 2.666667, 1.0e-10]])

run_test(
    file_name="DFEM_Diffusion/dDiffusion_2D_2a_DirBCs",
    comment="2D Diffusion with Dirichlet BC",
    num_procs=1,
    search_strings_vals_tols=[["[0]  Avg-value=", 0.295892, 1.0e-10]])

run_test(
    file_name="DFEM_Diffusion/dDiffusion_2D_2b_RobinBCs",
    comment="2D Diffusion with Robin BC",
    num_procs=1,
    search_strings_vals_tols=[["[0]  Avg-value=", 0.241757, 1.0e-10]])

run_test(
    file_name="DFEM_Diffusion/dDiffusion_2D_3a_analytical_coef",
    comment="2D Diffusion with Analytical Coefficients",
    num_procs=1,
    search_strings_vals_tols=[["[0]  Max-value=", 0.021923, 1.0e-10]])

run_test(
    file_name="DFEM_Diffusion/dDiffusion_2D_3b_analytical_coef2",
    comment="2D Diffusion with Manufactured Solution",
    num_procs=1,
    search_strings_vals_tols=[["[0]  Max-value=", 1.000586, 1.0e-10]])

# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#                                  Diffusion Tests
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Tests 25 - 36


run_test(
    file_name="Diffusion/Diffusion1D",
    comment="1D Diffusion Test - CFEM",
    num_procs=1,
    search_strings_vals_tols=[["[0]  Max-value=", 2.5, 1.0e-10]])

run_test(
    file_name="Diffusion/Diffusion1D_KBA",
    comment="1D Diffusion Test KBA partitioning - CFEM",
    num_procs=2,
    search_strings_vals_tols=[["[0]  Max-value=", 2.5, 1.0e-10]])

run_test(
    file_name="Diffusion/Diffusion1D_IP",
    comment="1D Diffusion Test - DFEM",
    num_procs=2,
    search_strings_vals_tols=[["[0]  Max-value=", 0.5006523128, 1.0e-4]])

run_test(
    file_name="Diffusion/Diffusion2D_1Poly",
    comment="2D Diffusion Test - CFEM",
    num_procs=1,
    search_strings_vals_tols=[["[0]  Max-value=", 0.29480, 1.0e-4]])

run_test(
    file_name="Diffusion/Diffusion2D_1Poly_IP",
    comment="2D Diffusion Test - DFEM",
    num_procs=4,
    search_strings_vals_tols=[["[0]  Max-value=", 2.5, 1.0e-4]])

run_test(
    file_name="Diffusion/Diffusion2D_2Unstructured",
    comment="2D Diffusion Test Unstr. Mesh - CFEM",
    num_procs=4,
    search_strings_vals_tols=[["[0]  Max-value=", 0.30384, 1.0e-4]])

run_test(
    file_name="Diffusion/Diffusion2D_2Unstructured_IP",
    comment="2D Diffusion Test Unstr. Mesh - DFEM",
    num_procs=4,
    search_strings_vals_tols=[["[0]  Max-value=", 0.29685, 1.0e-4]])

run_test(
    file_name="Diffusion/Diffusion3D_1Poly",
    comment="3D Diffusion Test - CFEM",
    num_procs=1,
    search_strings_vals_tols=[["[0]  Max-value=", 0.29480, 1.0e-4]])

run_test(
    file_name="Diffusion/Diffusion3D_1Poly_IP",
    comment="3D Diffusion Test - DFEM",
    num_procs=4,
    search_strings_vals_tols=[["[0]  Max-value=", 0.29492, 1.0e-4]])

run_test(
    file_name="Diffusion/Diffusion3D_2Ortho",
    comment="3D Diffusion Test Ortho Mesh - CFEM",
    num_procs=1,
    search_strings_vals_tols=[["[0]  Max-value=", 0.29480, 1.0e-4]])

run_test(
    file_name="Diffusion/Diffusion3D_3Unstructured",
    comment="3D Diffusion Test Unstr. Mesh - CFEM",
    num_procs=4,
    search_strings_vals_tols=[["[0]  Max-value=", 0.29499, 1.0e-4]])

run_test(
    file_name="Diffusion/Diffusion3D_3Unstructured_IP",
    comment="3D Diffusion Test Unstr. Mesh - DFEM",
    num_procs=4,
    search_strings_vals_tols=[["[0]  Max-value=", 0.29632, 1.0e-4]])

# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#                     Steady State Transport Tests
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Tests 37 - 53

run_test(
    file_name="Transport_Steady/Transport1D_1",
    comment="1D LinearBSolver Test - PWLD",
    num_procs=3,
    search_strings_vals_tols=[["[0]  Max-value1=", 0.49903, 1.0e-4],
                              ["[0]  Max-value2=", 7.18243e-4, 1.0e-4]])

run_test(
    file_name="Transport_Steady/Transport1D_3a_DSA_ortho",
    comment="1D LinearBSolver test of a block of graphite "
            "with an air cavity. DSA and TG",
    num_procs=4,
    search_strings_vals_tols=[
        ["StrCompare", "WGS groups [0-62] Iteration    28", 7, "CONVERGED"],
        ["StrCompare", "WGS groups [63-167] Iteration    39", 7, "CONVERGED"],
        ["NumCompare", "WGS groups [0-62] Iteration    28", 6, "float",
         6.74299e-07, 1.0e-9],
        ["NumCompare", "WGS groups [63-167] Iteration    39", 6, "float",
         8.73816e-07, 1.0e-9]])

run_test(
    file_name="Transport_Steady/Transport2D_1Poly",
    comment="2D LinearBSolver Test - PWLD",
    num_procs=4,
    search_strings_vals_tols=[["[0]  Max-value1=", 0.50758, 1.0e-4],
                              ["[0]  Max-value2=", 2.52527e-04, 1.0e-4]])

run_test(
    file_name="Transport_Steady/Transport2D_2Unstructured",
    comment="2D LinearBSolver Test Unstructured grid - PWLD",
    num_procs=4,
    search_strings_vals_tols=[["[0]  Max-value1=", 0.51187, 1.0e-4],
                              ["[0]  Max-value2=", 1.42458e-03, 1.0e-4]])

run_test(
    file_name="Transport_Steady/Transport2D_3Poly_quad_mod",
    comment="2D LinearBSolver Test Polar-Optimized quadrature - PWLD",
    num_procs=4,
    search_strings_vals_tols=[["[0]  Max-value1=", 0.50758, 1.0e-4],
                              ["[0]  Max-value2=", 2.52527e-04, 1.0e-4]])

run_test(
    file_name="Transport_Steady/Transport2D_4a_DSA_ortho",
    comment="2D LinearBSolver test of a block of graphite "
            "with an air cavity. DSA and TG",
    num_procs=4,
    search_strings_vals_tols=[
        ["StrCompare", "WGS groups [0-62] Iteration    53", 7, "CONVERGED"],
        ["StrCompare", "WGS groups [63-167] Iteration    59", 7, "CONVERGED"],
        ["NumCompare", "WGS groups [0-62] Iteration    53", 6, "float",
         5.96018e-07, 1.0e-9],
        ["NumCompare", "WGS groups [63-167] Iteration    59", 6, "float",
         5.96296e-07, 1.0e-9]])

run_test(
    file_name="Transport_Steady/Transport2D_4b_DSA_ortho",
    comment="2D LinearBSolver test of a block of graphite "
            "with an air cavity. DSA and TG",
    num_procs=4,
    search_strings_vals_tols=[
        ["StrCompare", "WGS groups [0-62] Iteration    54", 7, "CONVERGED"],
        ["StrCompare", "WGS groups [63-167] Iteration    56", 7, "CONVERGED"],
        ["NumCompare", "WGS groups [0-62] Iteration    54", 6, "float",
         5.00021e-07, 1.0e-9],
        ["NumCompare", "WGS groups [63-167] Iteration    56", 6, "float",
         9.73954e-07, 1.0e-9]])

run_test(
    file_name="Transport_Steady/Transport2D_5PolyA_AniHeteroBndry",
    comment="2D LinearBSolver Test Anisotropic Hetero BC - PWLD",
    num_procs=4,
    search_strings_vals_tols=[["[0]  Max-value1=", 3.18785, 1.0e-5]])

run_test(
    file_name="Transport_Steady/Transport3D_1a_Extruder",
    comment="3D LinearBSolver Test - PWLD",
    num_procs=4,
    search_strings_vals_tols=[["[0]  Max-value1=", 5.27450e-01, 1.0e-4],
                              ["[0]  Max-value2=", 3.76339e-04, 1.0e-4]])

run_test(
    file_name="Transport_Steady/Transport3D_1b_Ortho",
    comment="3D LinearBSolver Test - PWLD Reflecting BC",
    num_procs=4,
    search_strings_vals_tols=[["[0]  Max-value1=", 5.28310e-01, 1.0e-4],
                              ["[0]  Max-value2=", 8.04576e-04, 1.0e-4]])

run_test(
    file_name="Transport_Steady/Transport3D_1Poly_parmetis",
    comment="3D LinearBSolver Test Ortho Grid Parmetis - PWLD",
    num_procs=4,
    search_strings_vals_tols=[["[0]  Max-value1=", 5.27450e-01, 1.0e-4],
                              ["[0]  Max-value2=", 3.76339e-04, 1.0e-4]])

run_test(
    file_name="Transport_Steady/Transport3D_1Poly_qmom_part1",
    comment="3D LinearBSolver Test Source moment writing - PWLD",
    num_procs=4,
    search_strings_vals_tols=[["[0]  Max-value1=", 1.08320e-01, 1.0e-6],
                              ["[0]  Max-value2=", 0.000000000, 1.0e-10]])

run_test(
    file_name="Transport_Steady/Transport3D_1Poly_qmom_part2",
    comment="3D LinearBSolver Test Source moment reading - PWLD",
    num_procs=4,
    search_strings_vals_tols=[["[0]  Max-value1=", 1.01701e-04, 1.0e-6],
                              ["[0]  Max-value2=", 9.14681e-06, 1.0e-10]])

run_test(
    file_name="Transport_Steady/Transport3D_2Unstructured",
    comment="3D LinearBSolver Test Extruded Unstructured - PWLD",
    num_procs=4,
    search_strings_vals_tols=[["[0]  Max-value1=", 5.41465e-01, 1.0e-4],
                              ["[0]  Max-value2=", 3.78243e-04, 1.0e-4]])

run_test(
    file_name="Transport_Steady/Transport3D_3a_DSA_ortho",
    comment="3D LinearBSolver test of a block of graphite "
            "with an air cavity. DSA and TG",
    num_procs=4,
    search_strings_vals_tols=[
        ["StrCompare", "WGS groups [0-62] Iteration    54", 7, "CONVERGED"],
        ["StrCompare", "WGS groups [63-167] Iteration    63", 7, "CONVERGED"],
        ["NumCompare", "WGS groups [0-62] Iteration    54", 6, "float",
         7.92062e-07, 1.0e-9],
        ["NumCompare", "WGS groups [63-167] Iteration    63", 6, "float",
         8.1975e-07, 1.0e-9]])

run_test(
    file_name="Transport_Steady/Transport3D_4Cycles1",
    comment="3D LinearBSolver Test Extruded-Unstructured Mesh - PWLD",
    num_procs=4,
    search_strings_vals_tols=[["[0]  Max-value1=", 5.55349e-01, 1.0e-4],
                              ["[0]  Max-value2=", 3.74343e-04, 1.0e-4]])

run_test(
    file_name="Transport_Steady/Transport3D_5Cycles2",
    comment="3D LinearBSolver Test STAR-CCM+ mesh - PWLD",
    num_procs=4,
    search_strings_vals_tols=[["[0]  Max-value1=", 6.55396, 1.0e-4],
                              ["[0]  Max-value2=", 1.02943, 1.0e-4]])

# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#                     k-Eigenvalue Transport Tests
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Tests 54 - 56

run_test(
    file_name="Transport_Keigen/KEigenvalueTransport1D_1G",
    comment="1D KSolver LinearBSolver Test - PWLD",
    num_procs=4,
    search_strings_vals_tols=[
        # ["[0]          Final k-eigenvalue    :", 0.99954, 1.0e-5],
        ["NumCompare", "Final k-eigenvalue", 3, "float", 0.99954, 1.0e-5]])

run_test(
    file_name="Transport_Keigen/KEigenvalueTransport2D_1a_QBlock",
    comment="2D 2G KEigenvalue::Solver test using Power Iteration",
    num_procs=4,
    search_strings_vals_tols=[
        ["StrCompare", "Iteration    21", 9, "CONVERGED"],
        ["NumCompare", "Final k-eigenvalue", 3, "float", 0.5969127, 1.0e-7] ])

run_test(
    file_name="Transport_Keigen/KEigenvalueTransport2D_1b_QBlock",
    comment="2D 2G KEigenvalue::Solver test using NonLinearK",
    num_procs=4,
    search_strings_vals_tols=[
        ["NumCompare", "Iteration     3", 5, "float", 0.5969127, 1.0e-7],
        ["NumCompare", "Final k-eigenvalue", 3, "float", 0.5969127, 1.0e-7] ])

# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#         Steady-State Cylindrical Transport Tests
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Tests 57 - 58

run_test(
    file_name="Transport_Steady_Cyl/Transport2DCyl_1Monoenergetic",
    comment="2D LinearBSolver Cylindrical Test mono-energetic - PWLD",
    num_procs=4,
    search_strings_vals_tols=[["[0]  Max-value=", 1.00000, 1.0e-09]])

run_test(
    file_name="Transport_Steady_Cyl/Transport2DCyl_2Multigroup",
    comment="2D LinearBSolver Cylindrical Test multi-group - PWLD",
    num_procs=4,
    search_strings_vals_tols=[["[0]  Max-valueG1=", 1.00000, 1.0e-09],
                              ["[0]  Max-valueG2=", 0.25000, 1.0e-09]])

# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#             Steady-State Adjoint Transport Tests
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Tests 59 - 67

run_test(
    file_name="Transport_Adjoint/Adjoint2D_1a_forward",
    comment="2D Transport test with localized material source FWD",
    num_procs=4,
    search_strings_vals_tols=[["QOI-value=", 1.38397e-05, 1.0e-08]])

run_test(
    file_name="Transport_Adjoint/Adjoint2D_1b_adjoint",
    comment="2D Transport test with localized material source "
            "Adjoint generation",
    num_procs=4,
    search_strings_vals_tols=[])

run_test(
    file_name="Transport_Adjoint/Adjoint2D_1c_response",
    comment="2D Transport test with localized material source "
            "Adjoint inner product",
    num_procs=4,
    search_strings_vals_tols=[["Inner-product=", 1.38405e-05, 1.0e-08]])

run_test(
    file_name="Transport_Adjoint/Adjoint2D_2a_forward",
    comment="2D Transport test with point source FWD",
    num_procs=4,
    search_strings_vals_tols=[["QOI-value=", 2.90386e-05, 1.0e-08]])

run_test(
    file_name="Transport_Adjoint/Adjoint2D_2b_adjoint",
    comment="2D Transport test with point source Adjoint generation",
    num_procs=4,
    search_strings_vals_tols=[])

run_test(
    file_name="Transport_Adjoint/Adjoint2D_2c_response",
    comment="2D Transport test with point source Adjoint response",
    num_procs=4,
    search_strings_vals_tols=[["Inner-product=", 2.90543e-05, 1.0e-08]])

run_test(
    file_name="Transport_Adjoint/Adjoint2D_3a_forward",
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

run_test(
    file_name="Transport_Adjoint/Adjoint2D_3b_adjoint",
    comment="2D Transport test with point source Multigroup Adjoint generation",
    num_procs=4,
    search_strings_vals_tols=[])

run_test(
    file_name="Transport_Adjoint/Adjoint2D_3c_response",
    comment="2D Transport test with point source Multigroup Adjoint Response",
    num_procs=4,
    search_strings_vals_tols=[["Inner-product=", 3.30607e-06, 1.0e-09]])

# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#                        Transient Transport Tests
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ END OF TESTS
os.system("rm *.pvtu *.vtu *.e *.data")

print()
if num_failed == 0:
    print("All regression tests passed!")
else:
    print("\033[31mERROR: Not all regression tests passed!\033[39m")
    print("Number of Tests failed = " + str(num_failed))
print("")
print("************* End of Regression Test Suite *************")
print("")
if num_failed == 0:
    sys.exit(0)
else:
    sys.exit(1)
