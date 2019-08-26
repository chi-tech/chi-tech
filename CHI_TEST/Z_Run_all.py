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
print("Running Test " + str(test_number) + " " + test_name,end='')
process = subprocess.Popen([kpath_to_exe,
                            "CHI_TEST/Diffusion1D.lua"],
                           cwd=kchi_src_pth,
                           stdout=subprocess.PIPE,
                           universal_newlines=True)
process.wait()
out,err = process.communicate()

#string to find in output
find_str          = "Diffusion Solver: Number of iterations ="
#start of the string (<0 if not found)
test_str_start    = out.find(find_str)
#end of the string to find
test_str_end      = test_str_start + len(find_str)
#end of the line at which string was found
test_str_line_end = out.find("\n",test_str_start)

test_passed = False
if (test_str_start >= 0):
    #convert value to number
    test_val = int(out[test_str_end:test_str_line_end])
    if (test_val == 5):
        test_passed = True

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