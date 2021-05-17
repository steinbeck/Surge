CFLAGS=-I /usr/include/nauty/ -O3 -DWORDSIZE=32 -DMAXN=WORDSIZE -DOUTPROC=surgeproc -DZLIB -DPRUNE=surgeprune -DGENG_MAIN=geng_main

LDFLAGS=nautyW1.a -lz

OBJECTS=surge.o geng.o

surge: $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

.PHONY: clean

clean:
	$(RM) surge $(OBJECTS)
