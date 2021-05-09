
#CFLAGS=-I /usr/include/nauty/ -O3 -DWORDSIZE=32 -DMAXN=WORDSIZE -DOUTPROC=surgeproc -DZLIB -DPRUNE=surgeprune -DGENG_MAIN=geng_main

CFLAGS=-I /usr/include/nauty/ -O3 -DOUTPROC=surgeproc -DZLIB -DPRUNE=surgeprune -DGENG_MAIN=geng_main

LDFLAGS=-lnauty -lz

OBJECTS=surge.o geng.o

GENGC=/usr/share/nauty/geng.c

geng.o: $(GENGC)
	$(CC) $(CFLAGS) -c $< -o $@

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

surge: $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

.PHONY: clean

clean:
	$(RM) surge $(OBJECTS)
