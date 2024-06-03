CC=gcc
CFLAGS=-ansi -Wall -Wextra -Werror -pedantic-errors

all: symnmf

symnmf: symnmf.o
	$(CC) $(CFLAGS) -o symnmf symnmf.o -lm

symnmf.o: symnmf.c
	$(CC) $(CFLAGS) -c symnmf.c

clean:
	rm -f *.o symnmf