DEST = .
CC = gcc
CFLAGS = $(INCLUDE) $(LIB) 

JAVA_HOME = ..

JNI_LAPACK_SRC = $(JAVA_HOME)/src/jni_lapack
JNI_BLAS_SRC = $(JAVA_HOME)/src/jni_blas
JNI_RANDOM_SRC = $(JAVA_HOME)/src/jni_random
JNI_UTIL = $(JAVA_HOME)/src/util

INCLUDE = -I$(JAVA_HOME)/include -I/usr/include/gsl

LIB = -L/usr/local/lib #-L/usr/local/lib64 -L/usr/lib 

C_SRC_LAPACK1 = $(JNI_LAPACK_SRC)/c/LUDecomposition.c 
C_SRC_LAPACK2 = $(JNI_LAPACK_SRC)/c/CholeskyDecomposition.c
C_SRC_LAPACK3 = $(JNI_LAPACK_SRC)/c/QRDecomposition.c
C_SRC_LAPACK4 = $(JNI_LAPACK_SRC)/c/EigenvalueDecomposition.c
C_SRC_LAPACK5 = $(JNI_LAPACK_SRC)/c/SingularValueDecomposition.c
C_SRC_BLAS = $(JNI_BLAS_SRC)/c/Matrix.c
C_SRC_RANDOM = $(JNI_RANDOM_SRC)/c/JniGslRng.c
DYLIB_LAPACK1 = $(JAVA_HOME)/lib/libgsl_LUDecomposition.so
DYLIB_LAPACK2 = $(JAVA_HOME)/lib/libgsl_CholeskyDecomposition.so
DYLIB_LAPACK3 = $(JAVA_HOME)/lib/libgsl_QRDecomposition.so
DYLIB_LAPACK4 = $(JAVA_HOME)/lib/libgsl_EigenvalueDecomposition.so
DYLIB_LAPACK5 = $(JAVA_HOME)/lib/libgsl_SingularValueDecomposition.so
DYLIB_BLAS = $(JAVA_HOME)/lib/libgsl_Gsl.so
DYLIB_RANDOM = $(JAVA_HOME)/lib/libgsl_random.so

	

# For package JAMAJniGsl
JC = javac
JFLAGS = -d
JAVA_SRC = $(JNI_UTIL)/Maths.java  $(JNI_BLAS_SRC)/java/Matrix.java $(JNI_LAPACK_SRC)/java/LUDecomposition.java $(JNI_LAPACK_SRC)/java/CholeskyDecomposition.java $(JNI_LAPACK_SRC)/java/QRDecomposition.java $(JNI_LAPACK_SRC)/java/EigenvalueDecomposition.java $(JNI_LAPACK_SRC)/java/SingularValueDecomposition.java $(JNI_RANDOM_SRC)/java/JniGslRng.java
PACKAGE = $(JAVA_HOME)/class/JAMAJniGsl


#If you would like to link to openblas library, you can justly simply replace all the -llapack and -lblas below by the -lopenblas.
all: $(PACKAGE) $(DYLIB_LAPACK1) $(DYLIB_LAPACK2) $(DYLIB_LAPACK3) $(DYLIB_LAPACK4) $(DYLIB_LAPACK5) $(DYLIB_BLAS) $(DYLIB_RANDOM)

$(PACKAGE): $(JAVA_SRC)
	@mkdir -p $(@D)
	$(JC) $(JFLAGS) $(JAVA_HOME)/class $(JAVA_SRC)

$(DYLIB_LAPACK1): $(C_SRC_LAPACK1)
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -shared -fPIC $(C_SRC_LAPACK1) -lgsl -lcblas -lm -o $(DYLIB_LAPACK1)

$(DYLIB_LAPACK2): $(C_SRC_LAPACK2)
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -shared -fPIC $(C_SRC_LAPACK2) -lgsl -lcblas -lm -o $(DYLIB_LAPACK2)
$(DYLIB_LAPACK3): $(C_SRC_LAPACK3)
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -shared -fPIC $(C_SRC_LAPACK3) -lgsl -lcblas -lm -o $(DYLIB_LAPACK3)

$(DYLIB_LAPACK4): $(C_SRC_LAPACK4)
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -shared -fPIC $(C_SRC_LAPACK4) -lgsl -lcblas -lm -o $(DYLIB_LAPACK4)

$(DYLIB_LAPACK5): $(C_SRC_LAPACK5)
	@mkdir -p $(@D)
	 $(CC) $(CFLAGS) -shared -fPIC $(C_SRC_LAPACK5) -lgsl -lcblas -lm -o $(DYLIB_LAPACK5)

$(DYLIB_BLAS): $(C_SRC_BLAS) 
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -shared -fPIC $(C_SRC_BLAS) -lgsl -lcblas -lm -o $(DYLIB_BLAS)

$(DYLIB_RANDOM): $(C_SRC_RANDOM)
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -shared -fPIC $(C_SRC_RANDOM) -lgsl -lcblas -lm -o $(DYLIB_RANDOM)


clean: 
	rm -fr $(JAVA_HOME)/class
	rm -fr $(JAVA_HOME)/lib
