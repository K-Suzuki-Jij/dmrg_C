PROGNAME := 1D_KLM_TVF_DMRG.out
SRCDIR   := src
MAINDIR  := main
EXPDIR   := expectations
ALLOCDIR := allocations
STATDIR  := status
OUTDIR   := build
INCDIR   := ../../../include
LIBDIR   := ../../../lib
TARGET   := $(OUTDIR)/$(PROGNAME)
SRCS     := $(wildcard $(SRCDIR)/$(MAINDIR)/*.c) $(wildcard $(SRCDIR)/$(EXPDIR)/*.c) $(wildcard $(SRCDIR)/$(ALLOCDIR)/*.c) $(wildcard $(SRCDIR)/$(STATDIR)/*.c)
OBJS     := $(addprefix $(OUTDIR)/,$(patsubst %.c,%.o,$(SRCS)))

#CC = icc
#CFLAGS = -Wall -fast -parallel -qopenmp -I $(INCDIR)
#LFLAGS = -mkl

CC = gcc
CFLAGS = -Ofast -Xpreprocessor -fopenmp -I $(INCDIR) -I include
LFLAGS = -llapack -lomp -lsml -ldmrg -lonsite_op -L $(LIBDIR)

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

$(OUTDIR)/%.o:%.c
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) $(CFLAGS) -o $@ -c $<

clean:
	rm -rf $(OUTDIR)
