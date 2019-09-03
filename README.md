# JAMAJniGsl
JAMAJniGsl is a JAVA package providing a java interface for GSL library and using the classes defined by JAMA Package. It's a tool for numerical linear algebra computations and random number generators. JAMAJniG calls GSL libraries which require previous installation (Detailed installation instructions are listed below).


Build Instructions
------------------

* JAMAJniGsl requires the installation of GSL2.5 by using following commands.
 ```
wget ftp://ftp.gnu.org/gnu/gsl/gsl-2.5.tar.gz
cd gsl-2.5
mkdir /home/yourname/gsl
./configure --prefix=/home/yourname/gsl
make
make check
make install
 ```

* To compile the package, enter src directory and copy the makefile based on the operating system.

For Linux user, use the following commands:
```
cd src
cp Makefile.linux Makefile
Make
```

For Mac user, use the following commands:
```
cd src
cp Makefile.mac Makefile
Make
```

* Then execute `make`.

* To clean generated file, type `make clean` on the command line.  

Running the tests
-----------------
* For testing, enter test directory and execute `make` . If you want to clean testing results and all class files, type `make clean`. 

* There are four test files. The "JAMAJniGslTest.java" will test all the methods in JAMAJniGsl and report the errors. The "JAMAJniGslExamples.java" will provide specific examples for basic linear algebra operations. It can clearly show you how to use methods defined in JAMAJniGsl. If you are interested in how to use functions in GSL libraries to do matrix operations, the "JAMAJniGslExamplesBLAS.java" and "JAMAJniGslExamplesLAPACK.java"  and "GslRngtest.java" will give you specific examples. However, it is not necessary to go into blas and lapack if you just want to be a user of JAMAJniGsl.

Notes
---------
This package is intended for some basic problems we encounter when solving linear algebra problems. So we only include several most basic and widely used routines of the GSL library.


Source Repository
-----------------
JAMAJniGsl's source-code repository is hosted here on GitHub.


Authors
---------

| Name   | Email       |              |
|:------ |:----------- | :----------- |
| Yuzheng Dun (maintainer)| yuzhengdun@hust.edu.cn   | Visiting student, Department of Biostatistics  UCLA|
| Lu Zhang | lu.zhang@ucla.edu    | PhD student, Department of Biostatistics UCLA  |
| Sudipto Banerjee | sudipto@ucla.edu   | Professor, Department of Biostatistics  UCLA |
<!--- --->
                             


Licensing
---------
JAMAJniGsl is licensed under the Creative Commons Attribution License. 



