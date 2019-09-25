title: Test

Since it is not possible to perform unit test with current status of FEST-3D,
we have defined a few integrated test cases. Use following command to run all
tests in the `test/` directory.

```
$cd <FEST-3D/root/directory>/test/
$python Test.py
```
You will see following output on the screen:
```

  ----- Integrated Tests Started -----  
Total two processes will be used with MPICH library
Running Test number 1  --->  Subsonic flow over a smooth bump
Running Test number 2  --->  Laminar flow over a flat plate
Running Test number 3  --->  Turbulent flow over a flat plate
 ----- All tests completed -----
 
Tests passed:  3 out of 3
Check test summary in 'Report.txt' file.

```

The integrated tests use two processes with MPICH library. Three different test cases are defined:

1. Inviscid test case: [Subsonic flow over a 2D smooth bump](./05_tutorials/02_2dbump.html) (0.1% Tolerance).
2. Laminar test csae: [Laminar flow over a flat plate](./05_tutorials/04_LamFp.html) (1% Tolerance).
3. Turbulent test case: [Fully turbulent flow over a flat plate](./05_tutorials/05_TurbFp.html) (2% Tolerance. Higher tolerance due to coarser grid).

For more details about domain, boundary conditions and flow conditions of these test cases, check the Tutorial section. The solver setup (domain, grid, flow and boundary conditions, schemes) is same as listed in the separate tutorials.

Once the tests are complete, you can check the test summary in `Report.txt` file.

```
Ran Test number 1  --->  Subsonic flow over a smooth bump
  __________Report__________
 ---------- Inviscid Test case: Smooth Bump ----------
 Expected Change in entropy           : 0.000E+00
 Calculated relative change in entropy: 1.046E-06
 Difference                           : 1.046E-04 %
------------ >>> Test Passed  <<< --------------


Ran Test number 2  --->  Laminar flow over a flat plate
  __________Report__________
 ------ Laminar Test case: Flat plate ------
 Expected drag coeffcient    : 1.330E-03
 Calculated drag coefficient : 1.329E-03
 Difference                  : 4.638E-02 %
------------ >>> Test Passed  <<< --------------


Ran Test number 3  --->  Turbulent flow over a flat plate
  __________Report__________
 ------ Turbulent Test case: Flat plate ------
 Expected drag coeffcient    : 2.960E-03
 Calculated drag coefficient : 2.966E-03
 Difference                  : 1.883E-01 %
------------ >>> Test Passed  <<< --------------

```

