JAVA_HOME = ..

JAMAJNIGSL_CLASS = JAMAJniGslTest.class
JAMAJNIGSL_RUN = JAMAJniGsl_test_result.txt
BLAS_CLASS = JAMAJniGslExamplesBLAS.class
BLAS_RUN = BLAS_examples_result.txt
LAPACK_CLASS = JAMAJniGslExamplesLAPACK.class
LAPACK_RUN = LAPACK_examples_result.txt
EXAMPLE_CLASS = JAMAJniGslExamples.class
EXAMPLE_RUN = JAMAJniGsl_examples_result.txt
RNG_CLASS = GslRngtest.class
RNG_RUN = GslRng_examples_result.txt

all:  $(JAMAJNIGSL_CLASS) $(JAMAJNIGSL_RUN) $(BLAS_CLASS) $(BLAS_RUN) $(LAPACK_CLASS) $(LAPACK_RUN) $(EXAMPLE_CLASS) $(EXAMPLE_RUN) $(RNG_CLASS) $(RNG_RUN)

$(JAMAJNIGSL_CLASS):
	javac -classpath $(JAVA_HOME)/class:. JAMAJniGslTest.java 

$(JAMAJNIGSL_RUN):
	java -Djava.library.path=$(JAVA_HOME)/lib -classpath $(JAVA_HOME)/class:. JAMAJniGslTest > JAMAJniGsl_test_result.txt

$(BLAS_CLASS):
	javac -classpath $(JAVA_HOME)/class:. JAMAJniGslExamplesBLAS.java 

$(BLAS_RUN):
	java -Djava.library.path=$(JAVA_HOME)/lib -classpath $(JAVA_HOME)/class:. JAMAJniGslExamplesBLAS > BLAS_examples_result.txt

$(LAPACK_CLASS):
	javac -classpath $(JAVA_HOME)/class:. JAMAJniGslExamplesLAPACK.java 

$(LAPACK_RUN):
	java -Djava.library.path=$(JAVA_HOME)/lib -classpath $(JAVA_HOME)/class:. JAMAJniGslExamplesLAPACK > LAPACK_examples_result.txt

$(EXAMPLE_CLASS):
	javac -classpath $(JAVA_HOME)/class:. JAMAJniGslExamples.java 

$(EXAMPLE_RUN):
	java -Djava.library.path=$(JAVA_HOME)/lib -classpath $(JAVA_HOME)/class:. JAMAJniGslExamples > JAMAJniGslexamples_result.txt

$(RNG_CLASS):
	javac -classpath $(JAVA_HOME)/class:. GslRngtest.java

$(RNG_RUN):
	java -Djava.library.path=$(JAVA_HOME)/lib -classpath $(JAVA_HOME)/class:. GslRngtest > GslRng_examples_result.txt


clean: 
	rm -f *.txt *.class *.out *.serial
