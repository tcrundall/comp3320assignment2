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
OBJS_OPT = cloth_code_opt.o
OBJS_SSE = cloth_code_sse.o

all: 	opengl_main kernel_main kernel_main_opt kernel_main_sse

opengl_main:	opengl_main.o ${OBJS}
	$(CC) $(CFLAGS) -o $@ opengl_main.o ${OBJS} $(LDFLAGS) $(GLIBS)

kernel_main:	kernel_main.o ${OBJS}
	$(CC) $(CFLAGS) -o $@ kernel_main.o ${OBJS} $(LDFLAGS)

kernel_main_opt:	kernel_main_opt.o ${OBJS_OPT}
	$(CC) $(CFLAGS) -o $@ kernel_main_opt.o ${OBJS_OPT} $(LDFLAGS)

kernel_main_sse:	kernel_main_sse.o ${OBJS_SSE}
	$(CC) $(CFLAGS) -msse2 -o $@ kernel_main_sse.o ${OBJS_SSE} $(LDFLAGS)

.SUFFIXES:	.o .c

%.o : %.c	cloth_code.h cloth_param.h
	$(CC) $(CFLAGS) -c $<


clean:
	/bin/rm *.o *~ *_main *_opt *_sse




































