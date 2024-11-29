CC = g++
CFLAGS = -Wall -Werror -g -std=c++17
LDLIBS = -lm

BUILDDIR = ../build

# ------------------------------------------------------------------------------------
# Executable
# ------------------------------------------------------------------------------------
$(BUILDDIR)/geom3d: $(BUILDDIR)/lex.yy.o $(BUILDDIR)/geom3d_grammar.tab.o $(BUILDDIR)/geom3d_object.o $(BUILDDIR)/main.o \
                    $(BUILDDIR)/exponentiation.o $(BUILDDIR)/variable.o $(BUILDDIR)/term.o $(BUILDDIR)/polinom.o \
                    $(BUILDDIR)/point.o $(BUILDDIR)/line.o $(BUILDDIR)/solids.o $(BUILDDIR)/statement.o $(BUILDDIR)/sphere.o \
					$(BUILDDIR)/number.o $(BUILDDIR)/plane.o ndgConditions.o
	$(CC) $(CFLAGS) -o $@ $^

# ------------------------------------------------------------------------------------
# Scanner and Grammar compile
# ------------------------------------------------------------------------------------
$(BUILDDIR)/lex.yy.o: $(BUILDDIR)/lex.yy.c $(BUILDDIR)/geom3d_grammar.tab.hpp
	$(CC) $(CFLAGS) -c -o $@ $<

$(BUILDDIR)/geom3d_grammar.tab.o: $(BUILDDIR)/geom3d_grammar.tab.cpp
	$(CC) $(CFLAGS) -c -o $@ $<


$(BUILDDIR)/geom3d_grammar.tab.cpp $(BUILDDIR)/geom3d_grammar.tab.hpp: geom3d_grammar.ypp
	bison --defines=$(BUILDDIR)/geom3d_grammar.tab.hpp --output=$(BUILDDIR)/geom3d_grammar.tab.cpp -v geom3d_grammar.ypp 

$(BUILDDIR)/lex.yy.c: geom3d_scanner.l
	flex -o $(BUILDDIR)/lex.yy.c geom3d_scanner.l

$(BUILDDIR)/main.o: main.cpp main.h
	$(CC) $(CFLAGS) -c -o $@ $<	

# ------------------------------------------------------------------------------------
# Compile header files for geometry objects and statements
# ------------------------------------------------------------------------------------
$(BUILDDIR)/ndgConditions.o: ndgConditions.cpp ndgConditions.h
	$(CC) $(CFLAGS) -c -o $@ $<		

$(BUILDDIR)/statement.o: statement.cpp statement.h
	$(CC) $(CFLAGS) -c -o $@ $<	

$(BUILDDIR)/solids.o: solids.cpp solids.h
	$(CC) $(CFLAGS) -c -o $@ $<		

$(BUILDDIR)/line.o: line.cpp line.h
	$(CC) $(CFLAGS) -c -o $@ $<	

$(BUILDDIR)/plane.o: plane.cpp plane.h polinom.h
	$(CC) $(CFLAGS) -c -o $@ $<		

$(BUILDDIR)/point.o: point.cpp point.h
	$(CC) $(CFLAGS) -c -o $@ $<	

$(BUILDDIR)/sphere.o: sphere.cpp sphere.h
	$(CC) $(CFLAGS) -c -o $@ $<	

$(BUILDDIR)/number.o: number.cpp number.h
	$(CC) $(CFLAGS) -c -o $@ $<	

$(BUILDDIR)/geom3d_object.o: geom3d_object.cpp geom3d_object.h
	$(CC) $(CFLAGS) -c -o $@ $<	


# ------------------------------------------------------------------------------------
# Compile Variable and Polinome
# ------------------------------------------------------------------------------------
$(BUILDDIR)/polinom.o: polinom.cpp polinom.h
	$(CC) $(CFLAGS) -c -o $@ $<	

$(BUILDDIR)/term.o: term.cpp term.h
	$(CC) $(CFLAGS) -c -o $@ $<	

$(BUILDDIR)/exponentiation.o: exponentiation.cpp exponentiation.h
	$(CC) $(CFLAGS) -c -o $@ $<		

$(BUILDDIR)/variable.o: variable.cpp variable.h
	$(CC) $(CFLAGS) -c -o $@ $<	



.PHONY: clean


clean: 
	rm -fr $(BUILDDIR)/*