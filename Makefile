CC 		= g++
CFLAGS 	= -c
STDFLAGS = -std=c++11
EXNAME 	= Disk_Jamming

OBJS 	= main.o MD.o VerletList.o Energy.o

all: ${OBJS}
	${CC} ${STDFLAGS} -o ${EXNAME} ${OBJS}

main.o: main.cpp main.h
	${CC} ${STDFLAGS} ${CFLAGS} main.cpp
	
MD.o: MD.cpp MD.h
	${CC} ${STDFLAGS} ${CFLAGS} MD.cpp

VerletList.o: VerletList.cpp VerletList.h
	${CC} ${STDFLAGS} ${CFLAGS} VerletList.cpp

Energy.o: Energy.cpp Energy.h
	${CC} ${STDFLAGS} ${CFLAGS} Energy.cpp
	
clean:
	rm -f ${OBJS} !${EXNAME}
	@echo "all objects cleaned up!"