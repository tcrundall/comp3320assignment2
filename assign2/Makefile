CC	= gcc
LDFLAGS	= -lm 

#link for dependencies on Apple OS X
UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
GLIBS  = -framework OpenGL -framework GLUT
else
#link for dependencies for Linux
GLIBS   = -lGL -lGLU -lX11 -lglut -lXi -lXmu
endif

OBJS = cloth_code.o

all: 	opengl_main kernel_main

opengl_main:	opengl_main.o ${OBJS}
	$(CC) $(CFLAGS) -o $@ opengl_main.o ${OBJS} $(LDFLAGS) $(GLIBS)

kernel_main:	kernel_main.o ${OBJS}
	$(CC) $(CFLAGS) -o $@ kernel_main.o ${OBJS} $(LDFLAGS)

.SUFFIXES:	.o .c

%.o : %.c	cloth_code.h cloth_param.h
	$(CC) $(CFLAGS) -c $<


clean:
	/bin/rm a.out *.o *~ *_main




































