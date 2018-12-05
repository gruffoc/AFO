VERSION 	= 0
PATHLEVEL	= 1
SUBLEVEL	= 0
EXTRAVERSION 	= 
NAME		= AFO (Atmosphere Fluctuations For CMB Observations) - Stefano Mandelli

ESEGUIBILE	:= Signal.x
OBJ		:= Signal.o xorshift.o
CC		:= g++
FLAGS		:= -O3
LIBS		:=
INC		:=

%.o : %.cpp
	$(CC) ${LIBS} ${FLAGS} -c $< -o $@

${ESEGUIBILE}: ${OBJ}
	$(CC) ${LIBS} ${FLAGS} -o $@ $^	


.PHONY: clean

clean:
	rm -rf *.out *.x *.o


