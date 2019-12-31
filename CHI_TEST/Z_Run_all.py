
import subprocess
import os

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

print("")
print("************* ChiTech Regression Test *************")
print("")
test_number = 0
num_failed  = 0


#=========================================== Test
test_number += 1
test_name = "1D Diffusion Test - CFEM 1 MPI Process"
print("Running Test " + str(test_number) + " " + test_name,end='',flush=True)
process = subprocess.Popen([kpath_to_exe,
                            "CHI_TEST/Diffusion1D.lua", "master_export=false"],
                           cwd=kchi_src_pth,
                           stdout=subprocess.PIPE,
                           universal_newlines=True)
process.wait()
out,err = process.communicate()

#string to find in output
find_str          = "[0]  Max-value="
#start of the string (<0 if not found)
test_str_start    = out.find(find_str)
#end of the string to find
test_str_end      = test_str_start + len(find_str)
#end of the line at which string was found
test_str_line_end = out.find("\n",test_str_start)

test_passed = False
if (test_str_start >= 0):
    #convert value to number
    test_val = float(out[test_str_end:test_str_line_end])
    if (abs(test_val-2.5) < 1.0e-10):
        test_passed = True
else:
    test_passed = False

if (test_passed):
    print(" - Passed")
else:
    print(" - FAILED!")
    num_failed += 1

#=========================================== Test
test_number += 1
test_name = "1D Diffusion Test - MIP 2 MPI Processes"
print("Running Test " + str(test_number) + " " + test_name,end='',flush=True)
process = subprocess.Popen(["mpiexec","-np","2",kpath_to_exe,
                            "CHI_TEST/Diffusion1D_IP.lua", "master_export=false"],
                           cwd=kchi_src_pth,
                           stdout=subprocess.PIPE,
                           universal_newlines=True)
process.wait()
out,err = process.communicate()

#string to find in output
find_str          = "[0]  Max-value="
#start of the string (<0 if not found)
test_str_start    = out.find(find_str)
#end of the string to find
test_str_end      = test_str_start + len(find_str)
#end of the line at which string was found
test_str_line_end = out.find("\n",test_str_start)

test_passed = False
if (test_str_start >= 0):
    #convert value to number
    test_val = float(out[test_str_end:test_str_line_end])
    if (abs(test_val-0.5006523128) < 1.0e-4):
        test_passed = True
else:
    test_passed = False

if (test_passed):
    print(" - Passed")
else:
    print(" - FAILED!")
    num_failed += 1

#=========================================== Test
test_number += 1
test_name = "2D Diffusion Test - CFEM 1 MPI Processes"
print("Running Test " + str(test_number) + " " + test_name,end='',flush=True)
process = subprocess.Popen(["mpiexec","-np","1",kpath_to_exe,
                            "CHI_TEST/Diffusion2D_1Poly.lua", "master_export=false"],
                           cwd=kchi_src_pth,
                           stdout=subprocess.PIPE,
                           universal_newlines=True)
process.wait()
out,err = process.communicate()

#string to find in output
find_str          = "[0]  Max-value="
#start of the string (<0 if not found)
test_str_start    = out.find(find_str)
#end of the string to find
test_str_end      = test_str_start + len(find_str)
#end of the line at which string was found
test_str_line_end = out.find("\n",test_str_start)

test_passed = False
if (test_str_start >= 0):
    #convert value to number
    test_val = float(out[test_str_end:test_str_line_end])
    if (abs(test_val-0.29480) < 1.0e-4):
        test_passed = True
else:
    test_passed = False

if (test_passed):
    print(" - Passed")
else:
    print(" - FAILED!")
    num_failed += 1

#=========================================== Test
test_number += 1
test_name = "2D Diffusion Test - MIP 4 MPI Processes"
print("Running Test " + str(test_number) + " " + test_name,end='',flush=True)
process = subprocess.Popen(["mpiexec","-np","4",kpath_to_exe,
                            "CHI_TEST/Diffusion2D_1Poly_IP.lua", "master_export=false"],
                           cwd=kchi_src_pth,
                           stdout=subprocess.PIPE,
                           universal_newlines=True)
process.wait()
out,err = process.communicate()

#string to find in output
find_str          = "[0]  Max-value="
#start of the string (<0 if not found)
test_str_start    = out.find(find_str)
#end of the string to find
test_str_end      = test_str_start + len(find_str)
#end of the line at which string was found
test_str_line_end = out.find("\n",test_str_start)

test_passed = False
if (test_str_start >= 0):
    #convert value to number
    test_val = float(out[test_str_end:test_str_line_end])
    if (abs(test_val-2.5) < 1.0e-4):
        test_passed = True
else:
    test_passed = False

if (test_passed):
    print(" - Passed")
else:
    print(" - FAILED!")
    num_failed += 1

#=========================================== Test
test_number += 1
test_name = "3D Diffusion Test - CFEM 1 MPI Processes"
print("Running Test " + str(test_number) + " " + test_name,end='',flush=True)
process = subprocess.Popen(["mpiexec","-np","1",kpath_to_exe,
                            "CHI_TEST/Diffusion3D_1Poly.lua", "master_export=false"],
                           cwd=kchi_src_pth,
                           stdout=subprocess.PIPE,
                           universal_newlines=True)
process.wait()
out,err = process.communicate()

#string to find in output
find_str          = "[0]  Max-value="
#start of the string (<0 if not found)
test_str_start    = out.find(find_str)
#end of the string to find
test_str_end      = test_str_start + len(find_str)
#end of the line at which string was found
test_str_line_end = out.find("\n",test_str_start)

test_passed = False
if (test_str_start >= 0):
    #convert value to number
    test_val = float(out[test_str_end:test_str_line_end])
    if (abs(test_val-0.29480) < 1.0e-4):
        test_passed = True
else:
    test_passed = False

if (test_passed):
    print(" - Passed")
else:
    print(" - FAILED!")
    num_failed += 1

#=========================================== Test
test_number += 1
test_name = "3D Diffusion Test - MIP 4 MPI Processes"
print("Running Test " + str(test_number) + " " + test_name,end='',flush=True)
process = subprocess.Popen(["mpiexec","-np","4",kpath_to_exe,
                            "CHI_TEST/Diffusion3D_1Poly_IP.lua", "master_export=false"],
                           cwd=kchi_src_pth,
                           stdout=subprocess.PIPE,
                           universal_newlines=True)
process.wait()
out,err = process.communicate()

#string to find in output
find_str          = "[0]  Max-value="
#start of the string (<0 if not found)
test_str_start    = out.find(find_str)
#end of the string to find
test_str_end      = test_str_start + len(find_str)
#end of the line at which string was found
test_str_line_end = out.find("\n",test_str_start)

test_passed = False
if (test_str_start >= 0):
    #convert value to number
    test_val = float(out[test_str_end:test_str_line_end])
    if (abs(test_val-0.29492) < 1.0e-4):
        test_passed = True
else:
    test_passed = False

if (test_passed):
    print(" - Passed")
else:
    print(" - FAILED!")
    num_failed += 1

#=========================================== Test
test_number += 1
test_name = "1D LinearBSolver Test - PWLD 3 MPI Processes"
print("Running Test " + str(test_number) + " " + test_name,end='',flush=True)
process = subprocess.Popen(["mpiexec","-np","3",kpath_to_exe,
                            "CHI_TEST/Transport1D_1.lua", "master_export=false"],
                           cwd=kchi_src_pth,
                           stdout=subprocess.PIPE,
                           universal_newlines=True)
process.wait()
out,err = process.communicate()

#string to find in output
find_str          = "[0]  Max-value1="
#start of the string (<0 if not found)
test_str_start    = out.find(find_str)
#end of the string to find
test_str_end      = test_str_start + len(find_str)
#end of the line at which string was found
test_str_line_end = out.find("\n",test_str_start)

test_passed = True
if (test_str_start >= 0):
    #convert value to number
    test_val = float(out[test_str_end:test_str_line_end])
    if (not abs(test_val-0.49903) < 1.0e-4):
        test_passed = False
else:
    test_passed = False

#string to find in output
find_str          = "[0]  Max-value2="
#start of the string (<0 if not found)
test_str_start    = out.find(find_str)
#end of the string to find
test_str_end      = test_str_start + len(find_str)
#end of the line at which string was found
test_str_line_end = out.find("\n",test_str_start)

# test_passed = True
if (test_str_start >= 0):
    #convert value to number
    test_val = float(out[test_str_end:test_str_line_end])
    if (not abs(test_val-7.18243e-4) < 1.0e-4):
        test_passed = False
else:
    test_passed = False

if (test_passed):
    print(" - Passed")
else:
    print(" - FAILED!")
    num_failed += 1

#=========================================== Test
test_number += 1
test_name = "2D LinearBSolver Test - PWLD 4 MPI Processes"
print("Running Test " + str(test_number) + " " + test_name,end='',flush=True)
process = subprocess.Popen(["mpiexec","-np","4",kpath_to_exe,
                            "CHI_TEST/Transport2D_1Poly.lua", "master_export=false"],
                           cwd=kchi_src_pth,
                           stdout=subprocess.PIPE,
                           universal_newlines=True)
process.wait()
out,err = process.communicate()

#string to find in output
find_str          = "[0]  Max-value1="
#start of the string (<0 if not found)
test_str_start    = out.find(find_str)
#end of the string to find
test_str_end      = test_str_start + len(find_str)
#end of the line at which string was found
test_str_line_end = out.find("\n",test_str_start)

test_passed = True
if (test_str_start >= 0):
    #convert value to number
    test_val = float(out[test_str_end:test_str_line_end])
    if (not abs(test_val-0.50758) < 1.0e-4):
        test_passed = False
else:
    test_passed = False

#string to find in output
find_str          = "[0]  Max-value2="
#start of the string (<0 if not found)
test_str_start    = out.find(find_str)
#end of the string to find
test_str_end      = test_str_start + len(find_str)
#end of the line at which string was found
test_str_line_end = out.find("\n",test_str_start)

# test_passed = True
if (test_str_start >= 0):
    #convert value to number
    test_val = float(out[test_str_end:test_str_line_end])
    if (not abs(test_val-2.52527e-04) < 1.0e-4):
        test_passed = False
else:
    test_passed = False

if (test_passed):
    print(" - Passed")
else:
    print(" - FAILED!")
    num_failed += 1

#=========================================== Test
test_number += 1
test_name = "3D LinearBSolver Test - PWLD 4 MPI Processes"
print("Running Test " + str(test_number) + " " + test_name,end='',flush=True)
process = subprocess.Popen(["mpiexec","-np","4",kpath_to_exe,
                            "CHI_TEST/Transport3D_1Poly.lua", "master_export=false"],
                           cwd=kchi_src_pth,
                           stdout=subprocess.PIPE,
                           universal_newlines=True)
process.wait()
out,err = process.communicate()

#string to find in output
find_str          = "[0]  Max-value1="
#start of the string (<0 if not found)
test_str_start    = out.find(find_str)
#end of the string to find
test_str_end      = test_str_start + len(find_str)
#end of the line at which string was found
test_str_line_end = out.find("\n",test_str_start)

test_passed = True
if (test_str_start >= 0):
    #convert value to number
    test_val = float(out[test_str_end:test_str_line_end])
    if (not abs(test_val-5.27450e-01) < 1.0e-4):
        test_passed = False
else:
    test_passed = False

#string to find in output
find_str          = "[0]  Max-value2="
#start of the string (<0 if not found)
test_str_start    = out.find(find_str)
#end of the string to find
test_str_end      = test_str_start + len(find_str)
#end of the line at which string was found
test_str_line_end = out.find("\n",test_str_start)

# test_passed = True
if (test_str_start >= 0):
    #convert value to number
    test_val = float(out[test_str_end:test_str_line_end])
    if (not abs(test_val-3.76339e-04) < 1.0e-4):
        test_passed = False
else:
    test_passed = False

if (test_passed):
    print(" - Passed")
else:
    print(" - FAILED!")
    num_failed += 1

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ END OF TESTS
print("")
if (num_failed == 0):
    print("All regression tests passed!")
else:
    print("ERROR: Not all regression tests passed!")
    print("Number of Tests failed = " + str(num_failed))
print("")
print("************* End of Regression Test *************")
print("")