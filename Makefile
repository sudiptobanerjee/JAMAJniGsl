JAVA_SRC = src
PACKAGE = class/JAMAJniGsl
SUBPACKAGES_CLASS_PATH = class
JC = javac
JFLAGS = -d

all: 
	$(MAKE) -C src
	$(MAKE) -C test

#package:

#test:

clean:
	$(MAKE) clean -C src
	$(MAKE) clean -C test

