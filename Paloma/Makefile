all:
	gcc -Wall -g -pg -mavx2 -mtune=native -march=native -ffast-math -funroll-loops  -fopenmp Patterson.c
old:
	gcc -Wall -g -pg -O3 -mtune=native -march=native -fopenmp Patterson_old.c
p2:
	gcc -Wall -g -pg -O3 -mtune=native -march=native -fopenmp P2.c
clang:
	clang -Wall -g -pg -O3 -mavx2 -mtune=native -march=native -ffast-math -funroll-loops Patterson.c

clean:
	rm -f a.out

