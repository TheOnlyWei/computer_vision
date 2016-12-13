########################################
##
## Makefile
## LINUX compilation
##
##############################################


#FLAGS
C++FLAG = -g -std=c++11

MATH_LIBS = -lm

EXEC_DIR=.


.cc.o:
	g++ $(C++FLAG) $(INCLUDES)  -c $< -o $@


#Including
INCLUDES=  -I.

#-->All libraries (without LEDA)
LIBS_ALL =  -L/usr/lib -L/usr/local/lib

Cpp_OBJ1=image.o disjoint_sets.o s1.o
Cpp_OBJ2=image.o disjoint_sets.o s2.o
Cpp_OBJ3=image.o disjoint_sets.o s3.o
Cpp_OBJ4=image.o disjoint_sets.o s4.o
#Cpp_OBJ5=image.o disjoint_sets.o test.o
#Cpp_OBJ6=image.o disjoint_sets.o hough_section.o

PROGRAM1=s1
PROGRAM2=s2
PROGRAM3=s3
PROGRAM4=s4
#PROGRAM5=test
#PROGRAM6=hough_section

$(PROGRAM1): $(Cpp_OBJ1)
	g++ $(C++FLAG) -o $(EXEC_DIR)/$@ $(Cpp_OBJ1) $(INCLUDES) $(LIBS_ALL)
$(PROGRAM2): $(Cpp_OBJ2)
	g++ $(C++FLAG) -o $(EXEC_DIR)/$@ $(Cpp_OBJ2) $(INCLUDES) $(LIBS_ALL)
$(PROGRAM3): $(Cpp_OBJ3)
	g++ $(C++FLAG) -o $(EXEC_DIR)/$@ $(Cpp_OBJ3) $(INCLUDES) $(LIBS_ALL)
$(PROGRAM4): $(Cpp_OBJ4)
	g++ $(C++FLAG) -o $(EXEC_DIR)/$@ $(Cpp_OBJ4) $(INCLUDES) $(LIBS_ALL)
#$(PROGRAM5): $(Cpp_OBJ5)
	#g++ $(C++FLAG) -o $(EXEC_DIR)/$@ $(Cpp_OBJ5) $(INCLUDES) $(LIBS_ALL)
#$(PROGRAM6): $(Cpp_OBJ6)
	#g++ $(C++FLAG) -o $(EXEC_DIR)/$@ $(Cpp_OBJ6) $(INCLUDES) $(LIBS_ALL)





all:
	make $(PROGRAM1) # make s1 creates only s1
	make $(PROGRAM2) # make s2 creates only s2
	make $(PROGRAM3) # make s3 creates only s3
	make $(PROGRAM4) # make s4 creates only s4

clean:
	(rm -f *.o;)

(:
