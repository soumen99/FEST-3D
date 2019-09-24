import os
from subprocess import call

def run(case,f):
    os.chdir(case)
    call(['bash', "run.sh"], stdout=f, stderr=f)
    os.chdir("../")


if __name__=="__main__":
    print " "
    print "  ----- Integrated Tests Started -----  "
    print "Total two processes will be used with MPICH library"
    with open("make.log", 'w+') as f:
        call(['make', 'all'], stdout=f, stderr=f)
    with open("Report.txt", "w+") as f:
        print "Running Test number 1  --->  Subsonic flow over a smooth bump"
        run('SmoothBump',f)
        print "Running Test number 2  --->  Laminar flow over a flat plate"
        run('Lfp', f)
        print "Running Test number 3  --->  Turbulent flow over a flat plate"
        #print "This test might take few minutes ..."
        run('Tfp', f)
        print " ----- All tests completed -----"
        print " "

    with open("Report.txt", 'r') as f:
        data = f.read()
        print "Tests passed: ", str(data.count("Passed"))+" out of 3"
        print "Check test summary in 'Report.txt' file."
        print " "
