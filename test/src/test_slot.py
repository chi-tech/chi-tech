import os
import subprocess


class TestSlot:
    """Data structure to hold information regarding the parallel execution
       of a test"""
    def __init__(self, test, argv):
        self.process = None
        self.test = test
        self.passed = False
        self.argv = argv

        self._Run()

    def _Run(self):
        """Protected method to actually run the test"""
        test = self.test
        self.test.submitted = True

        cmd = "mpiexec "
        cmd += "-np " + str(test.num_procs) + " "
        cmd += self.argv.exe + " "
        cmd += test.filename + " "
        cmd += "--suppress_color "
        cmd += "--supress_beg_end_timelog "
        cmd += "master_export=false "
        # cmd += f"> out/{test.filename}.out "
        # cmd += "2>&1 "
        self.process = subprocess.Popen(cmd,
                                        cwd=test.file_dir,
                                        shell=True,
                                        stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE,
                                        universal_newlines=True)
        test_path = os.path.relpath(test.file_dir + test.filename)
        # print("Submitting test " + test_path)

    def Probe(self):
        """Probes the test to see if it finished"""
        running = True
        test = self.test

        if test.ran:
            return False

        if self.process.poll() is not None:

            out, err = self.process.communicate()

            file = open(test.file_dir + f"out/{test.filename}.out", "w")
            file.write(out + "\n")
            file.write(err + "\n")
            file.close()

            if not self.test.ran:
                self.PerformChecks()

            self.test.ran = True
            # print("done with " + test.GetTestPath())

            running = False

        return running

    def PerformChecks(self):
        """Applies to check-suite for the test"""
        test = self.test
        passed = True
        output_filename = f"{test.file_dir}out/{test.filename}.out"

        error_code = self.process.returncode
        for check in self.test.checks:
            verbose = self.argv.verbose
            check_passed = check.PerformCheck(output_filename,
                                              error_code, verbose)
            passed = passed and check_passed

            check_annotations = check.GetAnnotations()
            for ann in check_annotations:
                test.annotations.append(ann)

        test_path = os.path.relpath(test.file_dir + test.filename)

        if not os.path.isfile(test_path):
            test.annotations.append("lua file missing")

        pad = 0
        if passed:
            self.passed = True
            message = "\033[32mPassed\033[0m"
            pad += 5 + 4
        else:
            self.passed = False
            message = "\033[31mFailed\033[0m"
            pad += 5 + 4

        prefix = "\033[33m[{:2d}]\033[0m".format(test.num_procs)
        pad += 5 + 4

        for annotation in test.annotations:
            message = f"\033[36m[{annotation}]\033[0m" + message
            pad += 5 + 4

        width = 120 - len(prefix + test_path) + pad
        message = message.rjust(width, ".")

        print(prefix + test_path + message)
