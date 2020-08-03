import subprocess

class SubOutput:
  def __init__(self,success,output,error):
    self.success = success
    self.output = output
    self.error = error

#######################################
# Runs a subprocess
def ExecSub(command,log,env_vars=None):
  success = True
  output = ""
  error = b"No Error"

  result = subprocess.Popen(command,
                            stdout=log,
                            stderr=subprocess.PIPE,
                            shell=True,
                            env=env_vars)
  output, error = result.communicate()

  if (result.returncode != 0):
    success = False

  return SubOutput(success,error.decode("utf8"),output)


#######################################
INCLUDES="-I../ -I../ChiConsole -I../ChiLua -I../ChiMath -I../ChiPhysics"
INCLUDES=INCLUDES+" -I../ChiTimer -I../ChiMesh -I../ChiMPI -I../ChiLog"

log_file = open("log.txt","w+")
file = open("YFileList.txt","r")

ExecSub("rm Output/*",log_file)

i=-1
for line in file:
  i+=1
  src_file = line.rstrip()
  dot = src_file.find(".")
  outfile = src_file[0:dot]+".out"
  cmd = "g++ "+INCLUDES+" "+src_file+" -o test"

  # Compile
  print("Unit test {:3d} ".format(i)+src_file,end='')
  result = ExecSub(cmd,log_file)
  success = True
  if not result.success:
    print(" Compilation Failed!")
    success = False
  else:
    # Run
    result = ExecSub("./test >> Output/"+outfile,log_file)
    if not result.success:
      print(" Execution Failed!")
      success = False
    #Check output
    else:
      result = ExecSub("diff Output/"+outfile+" OutputRef/"+outfile,log_file)
      if len(result.output) != 0:
        print(" Failed!")
        success = False

  if success:
    print(" Success.")

ExecSub("rm test",log_file)
file.close()
log_file.close()
