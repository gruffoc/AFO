VERSION 	= 0
PATHLEVEL	= 1
SUBLEVEL	= 0
EXTRAVERSION 	= 
NAME		= AUTO-Correlation Coefficient - Stefano Mandelli
## Switch #######################################
DEVICE		:= true
DEBUGMODE	:= true
#################################################


#Makefile for the Device ########################


OBJ		:= corr_coef.o funct.o
ESEGUIBILE	:= corr_coef.x
CC		:= g++

FLAGS 		:= -O3  -g 
LIBS		:= -lgsl 
INC		:= 

%.o : %.cpp
	$(CC) ${LIBS} ${FLAGS} -c $< -o $@

${ESEGUIBILE}: ${OBJ}
	$(CC) ${LIBS} ${FLAGS} -o $@ $^	


.PHONY: clean

clean:
	rm -rf *.out *.x *.o


