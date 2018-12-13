VERSION 	= 0
PATHLEVEL	= 1
SUBLEVEL	= 0
EXTRAVERSION 	= 
NAME		= N-dim Integration - Stefano Mandelli

ESEGUIBILE	:= Integ.x
OBJ		:= main.o Integrate.o funct.o
CC		:= g++
FLAGS		:= -O3
LIBS		:=
INC		:=

%.o : %.cpp
	$(CC) ${LIBS} ${FLAGS} -c $< -o $@

${ESEGUIBILE}:${OBJ}
	$(CC) ${LIBS} ${FLAGS} -o $@ $^ 
	mv *.o obj/.


.PHONY: clean

clean:
	rm -rf *.out *.x *.o


